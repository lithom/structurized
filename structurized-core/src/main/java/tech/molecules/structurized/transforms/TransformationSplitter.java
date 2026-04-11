package tech.molecules.structurized.transforms;

import java.util.*;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import tech.molecules.structurized.OpenChemLibUtil;
import tech.molecules.structurized.context.ContextLabelCodec;

/**
 * TransformationSplitter
 * ----------------------
 *
 * Implements the strict-MCS-based splitting of an A→B change into one or more independent
 * {@link TransformationGroup}s and builds canonical {@link TransformationSignature}s for each group.
 *
 * Assumptions:
 *  - A and B are pre-standardized (aromaticity/tautomer policy, salts removed, etc.).
 *  - A strict MCS mapping (atoms + bonds identical in the core) is provided by caller.
 *  - No bond-type mismatches inside the preserved core (strictness enforced upstream).
 *
 * Granularity is controlled by:
 *  - radiusR: context shell radius within the preserved core away from attachments (0..3 typical)
 *  - featureMask: reserved metadata field carried in the signature hash; current context encoding does
 *    not branch on this mask yet
 *  - expanded raw context radius: fixed larger shell used only as an auxiliary downstream-analysis
 *    representation; it is not part of the canonical signature hash
 *
 * Context shell encoding:
 *  - The stored context shell is a copied core-only fragment around the ordered attachment atoms.
 *  - Each copied shell atom receives a custom label that encodes parent-core annotations.
 *  - These annotations preserve parent-core ring/aromatic/pi/ring-size semantics even when the shell
 *    itself truncates a larger ring system.
 *  - Attachment atoms also encode pairwise relations to later attachments, which allows downstream
 *    recovery of multi-anchor distances and same-ring/system relationships from the IDCode alone.
 *
 * Usage (sketch):
 *  StereoMolecule A = ...; // normalized
 *  StereoMolecule B = ...; // normalized
 *  TransformationSplitter.MCSMap mcs = TransformationSplitter.MCSMap.fromAtoB(mapAtoB);
 *  List<TransformationGroup> groups =
 *      TransformationSplitter.splitIntoTransformations(A, B, mcs, 1, FeatureMask.DEFAULT);
 *
 *  for (var tg : groups) {
 *      TransformationSignature sig = tg.signature;
 *      System.out.println("sigId="+sig.sigId+" type="+tg.type);
 *  }
 */
public class TransformationSplitter {
    /**
     * Default radius of the optional expanded raw context shell.
     *
     * This larger shell is kept alongside the canonical context shell for downstream analysis that wants
     * to derive richer local descriptors directly from a broader neighborhood.
     */
    public static final int DEFAULT_EXPANDED_RAW_CONTEXT_RADIUS = 3;

    // ================================
    // ======= Public API ============
    // ================================

    /** Main entry: split A→B into independent transformation groups and build signatures. */
    public static List<TransformationGroup> splitIntoTransformations(
            StereoMolecule A,
            StereoMolecule B,
            MCSMap mcs,
            int radiusR,
            int featureMask
    ) {
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        Objects.requireNonNull(mcs, "mcs");
        validateMCSMap(A, B, mcs);
        if (radiusR < 0) radiusR = 0;

        // 1) Define Core sets (mapped atoms) on A and B
        BitSet coreA = new BitSet(A.getAllAtoms());
        BitSet coreB = new BitSet(B.getAllAtoms());
        for (int a = 0; a < A.getAllAtoms(); a++) {
            int b = mcs.mapAtoB[a];
            if (b >= 0) {
                coreA.set(a);
                coreB.set(b);
            }
        }

        // 2) Compute Outside components on each side
        List<BitSet> removedComps = connectedComponentsExcluding(A, coreA);
        List<BitSet> addedComps   = connectedComponentsExcluding(B, coreB);
        CoreGraphData coreGraph = buildCoreGraph(A, coreA);

        // 3) Attachment atoms on each side = core atoms with ≥1 neighbor outside
        BitSet attachA = coreBoundaryAttachments(A, coreA);
        BitSet attachB = coreBoundaryAttachments(B, coreB);

        // 4) Build tripartite graph R–T–A (Removed–Attachment–Added)
        TripartiteGraph G = new TripartiteGraph(removedComps.size(), countBits(coreA), addedComps.size());

        // We define unified T indices by the order of core atoms in A (mapped set). Map B core to same via MCS
        int[] coreIndexAtoT = new int[A.getAllAtoms()];
        Arrays.fill(coreIndexAtoT, -1);
        int Tcount = 0;
        for (int a = coreA.nextSetBit(0); a >= 0; a = coreA.nextSetBit(a+1)) {
            coreIndexAtoT[a] = Tcount++;
        }
        // Sanity: Tcount == number of core atoms

        // Edges R–T (A-side)
        for (int i = 0; i < removedComps.size(); i++) {
            BitSet comp = removedComps.get(i);
            BitSet neighCore = neighborsInSet(A, comp, coreA);
            for (int a = neighCore.nextSetBit(0); a >= 0; a = neighCore.nextSetBit(a+1)) {
                if (attachA.get(a)) {
                    int t = coreIndexAtoT[a];
                    if (t >= 0) G.addEdgeRT(i, t);
                }
            }
        }

        // Edges A–T (B-side)
        for (int j = 0; j < addedComps.size(); j++) {
            BitSet comp = addedComps.get(j);
            BitSet neighCoreB = neighborsInSet(B, comp, coreB);
            for (int b = neighCoreB.nextSetBit(0); b >= 0; b = neighCoreB.nextSetBit(b+1)) {
                if (attachB.get(b)) {
                    int a = mcs.mapBtoA[b];
                    if (a >= 0) {
                        int t = coreIndexAtoT[a];
                        if (t >= 0) G.addEdgeAT(j, t);
                    }
                }
            }
        }

        // 5) Connected components in tripartite graph → independent transformation groups
        List<TripartiteComponent> comps = G.connectedComponents();

        int[] TtoA = new int[Tcount];
        Arrays.fill(TtoA, -1);
        for (int a = coreA.nextSetBit(0); a >= 0; a = coreA.nextSetBit(a+1)) {
            int t = coreIndexAtoT[a];
            TtoA[t] = a;
        }

        // 6) Build signatures for each component
        List<TransformationGroup> result = new ArrayList<>();
        for (TripartiteComponent cc : comps) {
            if (cc.R.isEmpty() && cc.A.isEmpty()) {
                continue;
            }

            // Collect attachments as actual A atom indices
            List<Integer> orderedT = new ArrayList<>(cc.T);
            Collections.sort(orderedT);
            List<Integer> attachmentsA = new ArrayList<>(orderedT.size());
            for (int t : orderedT) {
                attachmentsA.add(TtoA[t]);
            }

            // Removed fragment atoms (A side)
            BitSet Ratoms = new BitSet(A.getAllAtoms());
            for (int ridx : cc.R) Ratoms.or(removedComps.get(ridx));
            // Added fragment atoms (B side)
            BitSet Aatoms = new BitSet(B.getAllAtoms());
            for (int aidx : cc.A) Aatoms.or(addedComps.get(aidx));

            // Extract fragments with labeled attachment dummies in canonical order
            FragmentWithAttachments removedFrag = extractFragmentWithAttachments(A, Ratoms, attachmentsA);

            // Map attachments to B indices using MCS
            List<Integer> attachmentsB = new ArrayList<>(attachmentsA.size());
            for (int a : attachmentsA) attachmentsB.add(mcs.mapAtoB[a]);
            FragmentWithAttachments addedFrag = extractFragmentWithAttachments(B, Aatoms, attachmentsB);

            // Context shell within core (radiusR) around attachments
            BitSet coreContextA = coreContextWithinRadius(A, coreA, new HashSet<>(attachmentsA), radiusR);
            // Canonicalize context as attachment-aware IDCode of the induced subgraph (core-only atoms in coreContextA)
            String contextId = canonizeContextShell(A, coreContextA, attachmentsA, coreGraph);
            int expandedRawContextRadius = Math.max(radiusR, DEFAULT_EXPANDED_RAW_CONTEXT_RADIUS);
            BitSet expandedCoreContextA =
                    coreContextWithinRadius(A, coreA, new HashSet<>(attachmentsA), expandedRawContextRadius);
            String expandedRawContextId = canonizeContextShell(A, expandedCoreContextA, attachmentsA, coreGraph);

            // Attachment pattern (ordered)
            String attachPattern = attachmentPattern(A, attachmentsA);

            // Reaction class (heuristic string from fragments)
            String rxnClass = classifyReaction(removedFrag.mol, addedFrag.mol);

            // Build signature (canonical; deterministic order of attachments)
            TransformationSignature sig = new TransformationSignature(
                    removedFrag.idcode,
                    addedFrag.idcode,
                    attachPattern,
                    contextId,
                    expandedRawContextId,
                    radiusR,
                    expandedRawContextRadius,
                    featureMask,
                    rxnClass
            );

            // Classify TG type
            TransformationType type = classifyTG(cc);

            // Build group
            TransformationGroup tg = new TransformationGroup(sig, attachmentsA, type);
            result.add(tg);
        }

        return result;
    }

    // ================================
    // ======= Data Structures ========
    // ================================

    /** Strict MCS mapping: array mapAtoB and its inverse. Unmapped = -1. */
    public static class MCSMap {
        public final int[] mapAtoB;
        public final int[] mapBtoA;
        public MCSMap(int[] mapAtoB, int[] mapBtoA) {
            this.mapAtoB = mapAtoB;
            this.mapBtoA = mapBtoA;
        }
        public static MCSMap fromAtoB(int[] mapAtoB) {
            int maxB = -1;
            for (int v : mapAtoB) maxB = Math.max(maxB, v);
            int[] mapBtoA = new int[maxB+1];
            Arrays.fill(mapBtoA, -1);
            for (int a=0; a<mapAtoB.length; a++) {
                int b = mapAtoB[a];
                if (b>=0) mapBtoA[b]=a;
            }
            return new MCSMap(mapAtoB, mapBtoA);
        }
    }

    private static void validateMCSMap(StereoMolecule A, StereoMolecule B, MCSMap mcs) {
        if (mcs.mapAtoB.length != A.getAtoms()) {
            throw new IllegalArgumentException("mapAtoB length does not match A atom count");
        }
        if (mcs.mapBtoA.length < B.getAtoms()) {
            throw new IllegalArgumentException("mapBtoA length does not cover all B atoms");
        }

        boolean[] seenB = new boolean[B.getAtoms()];
        for (int a = 0; a < mcs.mapAtoB.length; a++) {
            int b = mcs.mapAtoB[a];
            if (b < 0) {
                continue;
            }
            if (b >= B.getAtoms()) {
                throw new IllegalArgumentException("mapAtoB contains out-of-range B atom index: " + b);
            }
            if (seenB[b]) {
                throw new IllegalArgumentException("mapAtoB is not one-to-one for B atom index: " + b);
            }
            seenB[b] = true;
            if (mcs.mapBtoA[b] != a) {
                throw new IllegalArgumentException("mapBtoA is inconsistent with mapAtoB at A atom index: " + a);
            }
        }
    }

    /**
     * Reserved bitmask for future mask-aware context abstractions.
     *
     * Today these flags are carried as signature metadata and included in the stable hash,
     * but the raw context shell encoding remains the same regardless of mask choice.
     * Downstream algorithms should derive context features via
     * tech.molecules.structurized.context.ContextFeatureExtractor.
     */
    public static class FeatureMask {
        public static final int AROMATIC     = 1 << 0;
        public static final int RING         = 1 << 1;
        public static final int HYBRID       = 1 << 2;
        public static final int FORMAL_CHG   = 1 << 3;
        public static final int H_COUNT      = 1 << 4;
        public static final int HET_NEIGHBOR = 1 << 5;
        public static final int AROM_PI      = 1 << 6;
        public static final int DEFAULT = AROMATIC | RING | HYBRID;
    }

    // ================================
    // ======= Internal helpers =======
    // ================================

    /**
     * Return connected components of the subgraph induced by atoms NOT in keepSet.
     */
    private static List<BitSet> connectedComponentsExcluding(StereoMolecule mol, BitSet keepSet) {
        int n = mol.getAllAtoms();
        BitSet visited = new BitSet(n);
        List<BitSet> comps = new ArrayList<>();
        for (int v=0; v<n; v++) {
            if (keepSet.get(v) || visited.get(v)) continue;
            BitSet comp = new BitSet(n);
            ArrayDeque<Integer> dq = new ArrayDeque<>();
            dq.add(v); visited.set(v); comp.set(v);
            while (!dq.isEmpty()) {
                int a = dq.poll();
                for (int i=0;i<mol.getConnAtoms(a);i++) {
                    int nb = mol.getConnAtom(a,i);
                    if (keepSet.get(nb) || visited.get(nb)) continue;
                    visited.set(nb); comp.set(nb); dq.add(nb);
                }
            }
            comps.add(comp);
        }
        return comps;
    }

    /** Core boundary attachments: core atoms with ≥1 neighbor outside core. */
    private static BitSet coreBoundaryAttachments(StereoMolecule mol, BitSet core) {
        BitSet att = new BitSet(mol.getAllAtoms());
        for (int a = core.nextSetBit(0); a >= 0; a = core.nextSetBit(a+1)) {
            for (int i=0;i<mol.getConnAtoms(a);i++) {
                int nb = mol.getConnAtom(a,i);
                if (!core.get(nb)) { att.set(a); break; }
            }
        }
        return att;
    }

    /** For a set S, return the set of neighbors that lie in the given targetSet. */
    private static BitSet neighborsInSet(StereoMolecule mol, BitSet S, BitSet targetSet) {
        BitSet out = new BitSet(mol.getAllAtoms());
        for (int a = S.nextSetBit(0); a >= 0; a = S.nextSetBit(a+1)) {
            for (int i=0;i<mol.getConnAtoms(a);i++) {
                int nb = mol.getConnAtom(a,i);
                if (targetSet.get(nb)) out.set(nb);
            }
        }
        return out;
    }

    /** Extract fragment induced by fragAtoms, plus create attachment dummies for the given ordered attachment atoms. */
    private static FragmentWithAttachments extractFragmentWithAttachments(StereoMolecule mol, BitSet fragAtoms, List<Integer> attachmentsOrdered) {
        StereoMolecule frag = new StereoMolecule();
        boolean[] fragAtomsArr = OpenChemLibUtil.bitsetToBool(fragAtoms,mol.getAtoms());

        int[] mapOldToNew = new int[mol.getAllAtoms()];
        Arrays.fill(mapOldToNew, -1);


        mol.copyMoleculeByAtoms(frag,fragAtomsArr,true,mapOldToNew);
        frag.ensureHelperArrays(Molecule.cHelperCIP);
//        // Copy fragment atoms
//        for (int a = fragAtoms.nextSetBit(0); a >= 0; a = fragAtoms.nextSetBit(a+1)) {
////            int newIdx = frag.addAtom(mol.getAtomicNo(a));
////            frag.setAtomCharge(newIdx, mol.getAtomCharge(a));
////            frag.setAtomMass(newIdx, mol.getAtomMass(a));
////            frag.setAtomRadical(newIdx, mol.getAtomRadical(a));
////            frag.setAtomQueryFeature(newIdx, mol.getAtomQueryFeatures(a));
//            int newIdx = mol.copyAtom(frag,a,0,0);
//            mapOldToNew[a] = newIdx;
//        }
//        // Copy bonds internal to fragment
//        for (int a = fragAtoms.nextSetBit(0); a >= 0; a = fragAtoms.nextSetBit(a+1)) {
//            for (int i=0;i<mol.getConnAtoms(a);i++) {
//
//
//                int nb = mol.getConnAtom(a,i);
//                if (!fragAtoms.get(nb)) continue;
//                int bOld = mol.getConnBond(a,i);
//                int aNew = mapOldToNew[a];
//                int nbNew = mapOldToNew[nb];
//                if (aNew<0 || nbNew<0) continue;
//                if (aNew < nbNew) { // add each bond once
//                    int bNew = frag.addBond(aNew, nbNew, mol.getBondOrder(bOld));
//                    frag.setBondType(bNew, mol.getBondType(bOld));
//                    frag.setBondQueryFeature(bNew, mol.getBondQueryFeatures(bOld));
//                }
//            }
//        }
        // Create labeled attachment dummies in the given order
        int label = 0;
        for (Integer aCore : attachmentsOrdered) {
            // Find neighbors of aCore that are in fragAtoms (i.e., bonds crossing core↔outside). For each such neighbor in fragAtoms, add a dummy and bond it to that neighbor.
            for (int i=0;i<mol.getConnAtoms(aCore);i++) {
                int nb = mol.getConnAtom(aCore,i);
                if (!fragAtoms.get(nb)) continue;
                int nbNew = mapOldToNew[nb];
                int d = frag.addAtom(0); // 0 = wildcard/"*"
                frag.setAtomCustomLabel(d, "*" + label);
                frag.addBond(d, nbNew, Molecule.cBondTypeSingle);
            }
            label++;
        }
        frag.ensureHelperArrays(Molecule.cHelperRings);
        String idcode = new Canonizer(frag, Canonizer.ENCODE_ATOM_CUSTOM_LABELS).getIDCode();
        return new FragmentWithAttachments(frag, idcode);
    }

    private static class FragmentWithAttachments {
        final StereoMolecule mol; final String idcode;
        FragmentWithAttachments(StereoMolecule m, String id) { this.mol=m; this.idcode=id; }
    }

    /** Context within the preserved core at radius r from attachments. */
    private static BitSet coreContextWithinRadius(StereoMolecule A, BitSet coreA, Set<Integer> attachmentsA, int r) {
        BitSet ctx = new BitSet(A.getAllAtoms());
        ArrayDeque<Integer> q = new ArrayDeque<>();
        int[] dist = new int[A.getAllAtoms()]; Arrays.fill(dist, -1);
        for (int a : attachmentsA) { ctx.set(a); dist[a]=0; q.add(a); }
        while (!q.isEmpty()) {
            int v = q.poll();
            if (dist[v] == r) continue;
            for (int i=0;i<A.getConnAtoms(v);i++) {
                int nb = A.getConnAtom(v,i);
                if (!coreA.get(nb)) continue; // stay in core only
                if (dist[nb] != -1) continue;
                dist[nb] = dist[v] + 1; ctx.set(nb); q.add(nb);
            }
        }
        return ctx;
    }

    /**
     * Canonize the induced core subgraph specified by atoms in {@code keep} as an attachment-aware
     * IDCode.
     *
     * The copied shell fragment is structurally only the retained core subgraph. To preserve parent-core
     * semantics across truncation, every copied atom is annotated before canonization with a custom label
     * that encodes:
     *  - ordered attachment identity, if applicable
     *  - parent-core ring membership
     *  - parent-core aromaticity
     *  - parent-core atom pi value
     *  - parent-core smallest ring size
     *  - pairwise attachment relations for multi-anchor contexts
     *
     * The resulting IDCode is therefore a self-contained raw context object: it can be decoded later to
     * recover both the shell graph and these parent-derived annotations.
     */
    private static String canonizeContextShell(
            StereoMolecule mol,
            BitSet keep,
            List<Integer> attachmentsOrdered,
            CoreGraphData coreGraph
    ) {
        StereoMolecule sub = new StereoMolecule();
        boolean[] keepAtoms = OpenChemLibUtil.bitsetToBool(keep, mol.getAtoms());
        int[] map = new int[mol.getAllAtoms()];
        Arrays.fill(map, -1);

        mol.copyMoleculeByAtoms(sub, keepAtoms, true, map);

        Map<Integer, ContextLabelCodec.ContextAtomData> annotations =
                buildContextAnnotations(coreGraph, keep, attachmentsOrdered);

        for (int oldAtom = keep.nextSetBit(0); oldAtom >= 0; oldAtom = keep.nextSetBit(oldAtom + 1)) {
            int newAtom = map[oldAtom];
            if (newAtom < 0) {
                continue;
            }
            sub.setAtomCustomLabel(newAtom, ContextLabelCodec.encode(annotations.get(oldAtom)));
        }

        sub.ensureHelperArrays(Molecule.cHelperRings);
        return new Canonizer(sub, Canonizer.ENCODE_ATOM_CUSTOM_LABELS).getIDCode();
    }

    /**
     * Build parent-core annotations for all shell atoms and pairwise relation metadata for ordered
     * attachment atoms.
     *
     * Pairwise attachment relations are stored only on the lower-index attachment atom for each pair.
     * This keeps labels compact while remaining sufficient for downstream reconstruction.
     */
    private static Map<Integer, ContextLabelCodec.ContextAtomData> buildContextAnnotations(
            CoreGraphData coreGraph,
            BitSet keep,
            List<Integer> attachmentsOrdered
    ) {
        Map<Integer, List<ContextLabelCodec.PairRelation>> pairRelations = new HashMap<>();
        for (int attachmentIndex = 0; attachmentIndex < attachmentsOrdered.size(); attachmentIndex++) {
            pairRelations.put(attachmentsOrdered.get(attachmentIndex), new ArrayList<>());
        }

        for (int i = 0; i < attachmentsOrdered.size(); i++) {
            int oldAtomA = attachmentsOrdered.get(i);
            int coreAtomA = coreGraph.oldToCore()[oldAtomA];
            for (int j = i + 1; j < attachmentsOrdered.size(); j++) {
                int oldAtomB = attachmentsOrdered.get(j);
                int coreAtomB = coreGraph.oldToCore()[oldAtomB];
                pairRelations.get(oldAtomA).add(new ContextLabelCodec.PairRelation(
                        j,
                        sharesAny(coreGraph.ringMembership().get(coreAtomA), coreGraph.ringMembership().get(coreAtomB)),
                        sharesAny(coreGraph.aromaticSystemMembership().get(coreAtomA), coreGraph.aromaticSystemMembership().get(coreAtomB)),
                        sharesAny(coreGraph.ringSystemMembership().get(coreAtomA), coreGraph.ringSystemMembership().get(coreAtomB)),
                        coreGraph.coreMol().getPathLength(coreAtomA, coreAtomB)
                ));
            }
        }

        Map<Integer, Integer> attachmentIndexByAtom = new HashMap<>();
        for (int attachmentIndex = 0; attachmentIndex < attachmentsOrdered.size(); attachmentIndex++) {
            attachmentIndexByAtom.put(attachmentsOrdered.get(attachmentIndex), attachmentIndex);
        }

        Map<Integer, ContextLabelCodec.ContextAtomData> annotations = new HashMap<>();
        for (int oldAtom = keep.nextSetBit(0); oldAtom >= 0; oldAtom = keep.nextSetBit(oldAtom + 1)) {
            int coreAtom = coreGraph.oldToCore()[oldAtom];
            annotations.put(oldAtom, new ContextLabelCodec.ContextAtomData(
                    attachmentIndexByAtom.get(oldAtom),
                    coreGraph.coreMol().isRingAtom(coreAtom),
                    coreGraph.coreMol().isAromaticAtom(coreAtom),
                    coreGraph.coreMol().getAtomPi(coreAtom),
                    coreGraph.coreMol().getAtomRingSize(coreAtom),
                    pairRelations.getOrDefault(oldAtom, List.of())
            ));
        }
        return annotations;
    }

    /**
     * Copy the preserved core once into its own fragment representation and precompute ring / ring-system
     * membership on that complete core graph.
     *
     * These precomputed parent-core relationships are the source of truth for context-shell annotation.
     * They intentionally come from the full copied core rather than the later truncated shell fragment.
     */
    private static CoreGraphData buildCoreGraph(StereoMolecule mol, BitSet coreAtoms) {
        StereoMolecule coreMol = new StereoMolecule();
        boolean[] keepAtoms = OpenChemLibUtil.bitsetToBool(coreAtoms, mol.getAtoms());
        int[] oldToCore = new int[mol.getAllAtoms()];
        Arrays.fill(oldToCore, -1);
        mol.copyMoleculeByAtoms(coreMol, keepAtoms, true, oldToCore);
        coreMol.ensureHelperArrays(Molecule.cHelperRings);

        List<Set<Integer>> ringMembership = collectRingMembership(coreMol, false);
        List<Set<Integer>> ringSystemMembership = collectRingSystemMembership(coreMol, false);
        List<Set<Integer>> aromaticSystemMembership = collectRingSystemMembership(coreMol, true);
        return new CoreGraphData(coreMol, oldToCore, ringMembership, ringSystemMembership, aromaticSystemMembership);
    }

    private static List<Set<Integer>> collectRingMembership(StereoMolecule mol, boolean aromaticOnly) {
        RingCollection rings = mol.getRingSet();
        List<Set<Integer>> membership = new ArrayList<>(mol.getAllAtoms());
        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            membership.add(new HashSet<>());
        }

        for (int ring = 0; ring < rings.getSize(); ring++) {
            if (aromaticOnly && !rings.isAromatic(ring)) {
                continue;
            }
            int[] ringAtoms = rings.getRingAtoms(ring);
            for (int atom : ringAtoms) {
                membership.get(atom).add(ring);
            }
        }
        return membership;
    }

    private static List<Set<Integer>> collectRingSystemMembership(StereoMolecule mol, boolean aromaticOnly) {
        RingCollection rings = mol.getRingSet();
        int ringCount = rings.getSize();
        int[] component = new int[ringCount];
        Arrays.fill(component, -1);

        int componentIndex = 0;
        for (int ring = 0; ring < ringCount; ring++) {
            if (component[ring] != -1) {
                continue;
            }
            if (aromaticOnly && !rings.isAromatic(ring)) {
                continue;
            }

            ArrayDeque<Integer> queue = new ArrayDeque<>();
            queue.add(ring);
            component[ring] = componentIndex;

            while (!queue.isEmpty()) {
                int current = queue.removeFirst();
                for (int candidate = 0; candidate < ringCount; candidate++) {
                    if (component[candidate] != -1) {
                        continue;
                    }
                    if (aromaticOnly && !rings.isAromatic(candidate)) {
                        continue;
                    }
                    if (sharesBond(rings.getRingBonds(current), rings.getRingBonds(candidate))) {
                        component[candidate] = componentIndex;
                        queue.add(candidate);
                    }
                }
            }

            componentIndex++;
        }

        List<Set<Integer>> membership = new ArrayList<>(mol.getAllAtoms());
        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            membership.add(new HashSet<>());
        }

        for (int ring = 0; ring < ringCount; ring++) {
            if (component[ring] == -1) {
                continue;
            }
            for (int atom : rings.getRingAtoms(ring)) {
                membership.get(atom).add(component[ring]);
            }
        }

        return membership;
    }

    private static boolean sharesBond(int[] ringBondsA, int[] ringBondsB) {
        for (int bondA : ringBondsA) {
            for (int bondB : ringBondsB) {
                if (bondA == bondB) {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean sharesAny(Set<Integer> left, Set<Integer> right) {
        for (int value : left) {
            if (right.contains(value)) {
                return true;
            }
        }
        return false;
    }

    private record CoreGraphData(
            StereoMolecule coreMol,
            int[] oldToCore,
            List<Set<Integer>> ringMembership,
            List<Set<Integer>> ringSystemMembership,
            List<Set<Integer>> aromaticSystemMembership
    ) {}

    /** Simple attachment pattern string (ordered), can be extended. */
    private static String attachmentPattern(StereoMolecule mol, List<Integer> attachmentsA) {
        StringBuilder sb = new StringBuilder();
        for (int idx=0; idx<attachmentsA.size(); idx++) {
            int a = attachmentsA.get(idx);
            sb.append(idx).append(":");
            sb.append(mol.getAtomicNo(a)).append("/");
            sb.append(mol.isAromaticAtom(a)?"Ar":"sp"+mol.getAtomPi(a));
            sb.append("/");
            sb.append(mol.isRingAtom(a)?"R":"N");
            if (idx+1<attachmentsA.size()) sb.append(";");
        }
        return sb.toString();
    }

    /** Very light reaction class heuristic from fragments; replace with your own. */
    private static String classifyReaction(StereoMolecule removed, StereoMolecule added) {
        int atomsRem = removed.getAllAtoms();
        int atomsAdd = added.getAllAtoms();
        if (atomsRem==0 && atomsAdd>0) return "INSERT";
        if (atomsRem>0 && atomsAdd==0) return "DELETE";
        return "REPLACE";
    }

    /** CC→Type classification based on component sizes. */
    private static TransformationType classifyTG(TripartiteComponent cc) {
        int r = cc.R.size(), a = cc.A.size(), t = cc.T.size();
        if (r==1 && a==1) return TransformationType.REPLACEMENT;
        if (r==0 && a==1) return TransformationType.INSERTION;
        if (r==1 && a==0) return TransformationType.DELETION;
        if (r>=2 && a==1) return TransformationType.MERGE;
        if (r==1 && a>=2) return TransformationType.SPLIT;
        return TransformationType.MULTI_CENTER;
    }

    private static int countBits(BitSet bs) { return bs.cardinality(); }

    // ================================
    // ===== Tripartite Graph =========
    // ================================

    private static class TripartiteGraph {
        final int Rn, Tn, An;
        final List<Integer>[] RT; // adjacency from R to T
        final List<Integer>[] AT; // adjacency from A to T
        TripartiteGraph(int Rn, int Tn, int An) {
            this.Rn=Rn; this.Tn=Tn; this.An=An;
            //noinspection unchecked
            RT = new List[Rn];
            //noinspection unchecked
            AT = new List[An];
            for (int i=0;i<Rn;i++) RT[i]=new ArrayList<>();
            for (int j=0;j<An;j++) AT[j]=new ArrayList<>();
        }
        void addEdgeRT(int r, int t) { RT[r].add(t); }
        void addEdgeAT(int a, int t) { AT[a].add(t); }
        List<TripartiteComponent> connectedComponents() {
            // Build a combined graph over nodes indexed in three blocks: [0..Rn-1]=R, [Rn..Rn+Tn-1]=T, [Rn+Tn..Rn+Tn+An-1]=A
            int N = Rn+Tn+An;
            List<Integer>[] adj = new List[N];
            for (int i=0;i<N;i++) adj[i]=new ArrayList<>();
            // Undirected edges for CC purposes
            for (int r=0;r<Rn;r++) for (int t : RT[r]) { int R=r, T=Rn+t; adj[R].add(T); adj[T].add(R); }
            for (int a=0;a<An;a++) for (int t : AT[a]) { int A=Rn+Tn+a, T=Rn+t; adj[A].add(T); adj[T].add(A); }

            boolean[] vis = new boolean[N];
            List<TripartiteComponent> out = new ArrayList<>();
            for (int s=0;s<N;s++) if (!vis[s]) {
                ArrayDeque<Integer> dq=new ArrayDeque<>(); dq.add(s); vis[s]=true;
                TripartiteComponent cc = new TripartiteComponent();
                while(!dq.isEmpty()){
                    int v = dq.poll();
                    if (v < Rn) cc.R.add(v);
                    else if (v < Rn+Tn) cc.T.add(v - Rn);
                    else cc.A.add(v - (Rn+Tn));
                    for (int nb : adj[v]) if (!vis[nb]) { vis[nb]=true; dq.add(nb);}                }
                out.add(cc);
            }
            return out;
        }
    }

    private static class TripartiteComponent {
        final List<Integer> R = new ArrayList<>(); // indices into removedComps
        final List<Integer> T = new ArrayList<>(); // indices of core attachments (A-side ordering)
        final List<Integer> A = new ArrayList<>(); // indices into addedComps
    }
}


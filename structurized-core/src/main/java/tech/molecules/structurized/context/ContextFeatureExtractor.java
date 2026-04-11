package tech.molecules.structurized.context;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

/**
 * Derives structured context descriptors from the raw attachment-aware context shell IDCode.
 *
 * <p>The extractor first decodes the stored shell fragment and then reads any custom-label annotations
 * written by {@code TransformationSplitter}. When parent-aware annotations are present, they take
 * precedence over shell-local helper-array values for ring, aromatic, pi, ring-size, and attachment-pair
 * semantics. This allows a low-radius shell to retain chemically relevant information from the full
 * preserved core.</p>
 */
public final class ContextFeatureExtractor {
    private ContextFeatureExtractor() {}

    /**
     * Compute the full current context feature set from the canonical raw shell representation.
     */
    public static ContextFeatures computeAll(String contextShellIdcode) {
        DecodedContextShell decoded = decodeContextShell(contextShellIdcode);
        AttachmentPairFeatureSet pairFeatureSet = computeAttachmentPairFeatureSet(decoded);
        return new ContextFeatures(
                contextShellIdcode,
                computeShellSummary(decoded),
                computeAttachmentAtomFeatures(decoded),
                pairFeatureSet.pairFeatures(),
                pairFeatureSet.distanceMatrix(),
                pairFeatureSet.minDistance(),
                pairFeatureSet.maxDistance()
        );
    }

    /** Compute shell-global summary counts. */
    public static ContextFeatures.ShellSummaryFeatures computeShellSummary(String contextShellIdcode) {
        return computeShellSummary(decodeContextShell(contextShellIdcode));
    }

    /** Compute ordered per-attachment atom features. */
    public static List<ContextFeatures.AttachmentAtomFeatures> computeAttachmentAtomFeatures(String contextShellIdcode) {
        return computeAttachmentAtomFeatures(decodeContextShell(contextShellIdcode));
    }

    /** Compute pairwise relation features for all ordered attachment pairs. */
    public static List<ContextFeatures.AttachmentPairFeatures> computeAttachmentPairFeatures(String contextShellIdcode) {
        return computeAttachmentPairFeatureSet(decodeContextShell(contextShellIdcode)).pairFeatures();
    }

    private static ContextFeatures.ShellSummaryFeatures computeShellSummary(DecodedContextShell decoded) {
        StereoMolecule mol = decoded.mol();
        int heteroAtomCount = 0;
        int aromaticAtomCount = 0;
        int ringAtomCount = 0;
        int formalChargeSum = 0;
        int formalChargeAbsSum = 0;

        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            ContextLabelCodec.ContextAtomData annotation = decoded.annotationByAtom().get(atom);
            if (isHeteroAtom(mol.getAtomicNo(atom))) {
                heteroAtomCount++;
            }
            if (isAromatic(atom, annotation, mol)) {
                aromaticAtomCount++;
            }
            if (isRingAtom(atom, annotation, mol)) {
                ringAtomCount++;
            }
            formalChargeSum += mol.getAtomCharge(atom);
            formalChargeAbsSum += Math.abs(mol.getAtomCharge(atom));
        }

        return new ContextFeatures.ShellSummaryFeatures(
                mol.getAllAtoms(),
                mol.getAllBonds(),
                heteroAtomCount,
                aromaticAtomCount,
                ringAtomCount,
                formalChargeSum,
                formalChargeAbsSum,
                decoded.attachments().size()
        );
    }

    private static List<ContextFeatures.AttachmentAtomFeatures> computeAttachmentAtomFeatures(DecodedContextShell decoded) {
        StereoMolecule mol = decoded.mol();
        List<ContextFeatures.AttachmentAtomFeatures> features = new ArrayList<>(decoded.attachments().size());
        for (int attachmentIndex = 0; attachmentIndex < decoded.attachments().size(); attachmentIndex++) {
            int atom = decoded.attachments().get(attachmentIndex);
            ContextLabelCodec.ContextAtomData annotation = decoded.annotationByAtom().get(atom);
            int heteroNeighborCount = 0;
            for (int i = 0; i < mol.getConnAtoms(atom); i++) {
                int neighbor = mol.getConnAtom(atom, i);
                if (isHeteroAtom(mol.getAtomicNo(neighbor))) {
                    heteroNeighborCount++;
                }
            }

            features.add(new ContextFeatures.AttachmentAtomFeatures(
                    attachmentIndex,
                    mol.getAtomicNo(atom),
                    isAromatic(atom, annotation, mol),
                    isRingAtom(atom, annotation, mol),
                    mol.getAtomCharge(atom),
                    parentAtomPi(atom, annotation, mol),
                    smallestRingSize(atom, annotation, mol),
                    mol.getConnAtoms(atom),
                    heteroNeighborCount
            ));
        }
        return features;
    }

    private static AttachmentPairFeatureSet computeAttachmentPairFeatureSet(DecodedContextShell decoded) {
        StereoMolecule mol = decoded.mol();
        List<Integer> attachments = decoded.attachments();
        int[][] distanceMatrix = new int[attachments.size()][attachments.size()];
        for (int[] row : distanceMatrix) {
            Arrays.fill(row, -1);
        }
        for (int i = 0; i < attachments.size(); i++) {
            distanceMatrix[i][i] = 0;
        }

        List<Set<Integer>> ringMembership = collectRingMembership(decoded, false);
        List<Set<Integer>> fusedRingSystemMembership = collectRingSystemMembership(decoded, false);
        List<Set<Integer>> aromaticSystemMembership = collectRingSystemMembership(decoded, true);

        List<ContextFeatures.AttachmentPairFeatures> pairFeatures = new ArrayList<>();
        int minDistance = Integer.MAX_VALUE;
        int maxDistance = Integer.MIN_VALUE;
        boolean foundDistance = false;

        for (int i = 0; i < attachments.size(); i++) {
            for (int j = i + 1; j < attachments.size(); j++) {
                ContextFeatures.AttachmentPairFeatures pairFeature =
                        pairFeatureFromAnnotation(decoded, i, j, mol, ringMembership, aromaticSystemMembership, fusedRingSystemMembership);
                int distance = pairFeature.topologicalDistance();
                distanceMatrix[i][j] = distance;
                distanceMatrix[j][i] = distance;

                if (distance >= 0) {
                    foundDistance = true;
                    minDistance = Math.min(minDistance, distance);
                    maxDistance = Math.max(maxDistance, distance);
                }

                pairFeatures.add(pairFeature);
            }
        }

        if (!foundDistance) {
            minDistance = -1;
            maxDistance = -1;
        }

        return new AttachmentPairFeatureSet(pairFeatures, distanceMatrix, minDistance, maxDistance);
    }

    private static List<Set<Integer>> collectRingMembership(DecodedContextShell decoded, boolean aromaticOnly) {
        RingCollection rings = decoded.mol().getRingSet();
        List<Set<Integer>> membership = new ArrayList<>(decoded.attachments().size());
        for (int ignored : decoded.attachments()) {
            membership.add(new HashSet<>());
        }

        for (int ring = 0; ring < rings.getSize(); ring++) {
            if (aromaticOnly && !rings.isAromatic(ring)) {
                continue;
            }
            for (int attachmentIndex = 0; attachmentIndex < decoded.attachments().size(); attachmentIndex++) {
                if (rings.isAtomMember(ring, decoded.attachments().get(attachmentIndex))) {
                    membership.get(attachmentIndex).add(ring);
                }
            }
        }

        return membership;
    }

    /**
     * Recover pairwise attachment relations from parent-aware annotations when present, falling back to
     * shell-local graph calculations otherwise.
     */
    private static ContextFeatures.AttachmentPairFeatures pairFeatureFromAnnotation(
            DecodedContextShell decoded,
            int attachmentIndex1,
            int attachmentIndex2,
            StereoMolecule mol,
            List<Set<Integer>> ringMembership,
            List<Set<Integer>> aromaticSystemMembership,
            List<Set<Integer>> fusedRingSystemMembership
    ) {
        ContextLabelCodec.ContextAtomData annotation = decoded.annotationByAtom().get(decoded.attachments().get(attachmentIndex1));
        if (annotation != null) {
            for (ContextLabelCodec.PairRelation pairRelation : annotation.pairRelations()) {
                if (pairRelation.otherAttachmentIndex() == attachmentIndex2) {
                    return new ContextFeatures.AttachmentPairFeatures(
                            attachmentIndex1,
                            attachmentIndex2,
                            pairRelation.topologicalDistance(),
                            pairRelation.sameRing(),
                            pairRelation.sameAromaticSystem(),
                            pairRelation.sameFusedRingSystem()
                    );
                }
            }
        }

        return new ContextFeatures.AttachmentPairFeatures(
                attachmentIndex1,
                attachmentIndex2,
                mol.getPathLength(decoded.attachments().get(attachmentIndex1), decoded.attachments().get(attachmentIndex2)),
                intersects(ringMembership.get(attachmentIndex1), ringMembership.get(attachmentIndex2)),
                intersects(aromaticSystemMembership.get(attachmentIndex1), aromaticSystemMembership.get(attachmentIndex2)),
                intersects(fusedRingSystemMembership.get(attachmentIndex1), fusedRingSystemMembership.get(attachmentIndex2))
        );
    }

    private static List<Set<Integer>> collectRingSystemMembership(DecodedContextShell decoded, boolean aromaticOnly) {
        RingCollection rings = decoded.mol().getRingSet();
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

        List<Set<Integer>> membership = new ArrayList<>(decoded.attachments().size());
        for (int ignored : decoded.attachments()) {
            membership.add(new HashSet<>());
        }

        for (int ring = 0; ring < ringCount; ring++) {
            if (component[ring] == -1) {
                continue;
            }
            for (int attachmentIndex = 0; attachmentIndex < decoded.attachments().size(); attachmentIndex++) {
                if (rings.isAtomMember(ring, decoded.attachments().get(attachmentIndex))) {
                    membership.get(attachmentIndex).add(component[ring]);
                }
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

    private static boolean intersects(Set<Integer> a, Set<Integer> b) {
        for (int value : a) {
            if (b.contains(value)) {
                return true;
            }
        }
        return false;
    }

    private static boolean isHeteroAtom(int atomicNo) {
        return atomicNo > 1 && atomicNo != 6;
    }

    /**
     * Prefer parent-core aromaticity when encoded; otherwise fall back to shell-local helper arrays.
     */
    private static boolean isAromatic(int atom, ContextLabelCodec.ContextAtomData annotation, StereoMolecule mol) {
        return annotation != null ? annotation.parentAromaticAtom() : mol.isAromaticAtom(atom);
    }

    /** Prefer parent-core ring membership when encoded; otherwise fall back to shell-local helper arrays. */
    private static boolean isRingAtom(int atom, ContextLabelCodec.ContextAtomData annotation, StereoMolecule mol) {
        return annotation != null ? annotation.parentRingAtom() : mol.isRingAtom(atom);
    }

    /** Prefer parent-core pi value when encoded; otherwise fall back to shell-local helper arrays. */
    private static int parentAtomPi(int atom, ContextLabelCodec.ContextAtomData annotation, StereoMolecule mol) {
        return annotation != null ? annotation.parentAtomPi() : mol.getAtomPi(atom);
    }

    /** Prefer parent-core smallest ring size when encoded; otherwise fall back to shell-local helper arrays. */
    private static int smallestRingSize(int atom, ContextLabelCodec.ContextAtomData annotation, StereoMolecule mol) {
        return annotation != null ? annotation.parentSmallestRingSize() : mol.getAtomRingSize(atom);
    }

    /**
     * Decode the raw shell IDCode into a molecule plus ordered attachment references and any stored
     * per-atom annotations.
     */
    private static DecodedContextShell decodeContextShell(String contextShellIdcode) {
        Objects.requireNonNull(contextShellIdcode, "contextShellIdcode");
        if (contextShellIdcode.isEmpty()) {
            throw new IllegalArgumentException("contextShellIdcode must not be empty");
        }

        StereoMolecule mol = new StereoMolecule();
        new IDCodeParser().parse(mol, contextShellIdcode);
        mol.ensureHelperArrays(StereoMolecule.cHelperRings);

        List<ContextLabelCodec.ContextAtomData> annotationByAtom = new ArrayList<>(mol.getAllAtoms());
        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            annotationByAtom.add(null);
        }
        List<AttachmentRef> attachmentRefs = new ArrayList<>();
        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            ContextLabelCodec.ContextAtomData annotation = ContextLabelCodec.decode(mol.getAtomCustomLabel(atom));
            annotationByAtom.set(atom, annotation);
            if (annotation == null || annotation.attachmentIndex() == null) {
                continue;
            }
            attachmentRefs.add(new AttachmentRef(annotation.attachmentIndex(), atom));
        }

        attachmentRefs.sort((a, b) -> Integer.compare(a.attachmentIndex(), b.attachmentIndex()));
        List<Integer> attachments = new ArrayList<>(attachmentRefs.size());
        for (AttachmentRef attachmentRef : attachmentRefs) {
            attachments.add(attachmentRef.atom());
        }

        return new DecodedContextShell(mol, attachments, annotationByAtom);
    }

    private record DecodedContextShell(
            StereoMolecule mol,
            List<Integer> attachments,
            List<ContextLabelCodec.ContextAtomData> annotationByAtom
    ) {}

    private record AttachmentRef(int attachmentIndex, int atom) {}

    private record AttachmentPairFeatureSet(
            List<ContextFeatures.AttachmentPairFeatures> pairFeatures,
            int[][] distanceMatrix,
            int minDistance,
            int maxDistance
    ) {}
}

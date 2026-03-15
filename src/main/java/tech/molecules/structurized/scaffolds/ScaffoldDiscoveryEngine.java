package tech.molecules.structurized.scaffolds;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import tech.molecules.structurized.OpenChemLibUtil;
import tech.molecules.structurized.transforms.TransformationBench;
import tech.molecules.structurized.transforms.TransformationSplitter;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

/**
 * V1 scaffold discovery engine.
 *
 * <p>The engine uses fragment fingerprints only to build small local neighborhoods. Actual scaffold
 * candidates are generated from strict pairwise common cores inside those neighborhoods. Identical
 * scaffolds are merged canonically and then rescored using local scaffold-to-compound support.</p>
 */
public final class ScaffoldDiscoveryEngine {
    private ScaffoldDiscoveryEngine() {}

    public static ScaffoldDiscoveryResult discover(
            List<StereoMolecule> molecules,
            TransformationBench.StrictMCSProvider mcsProvider,
            ScaffoldDiscoveryConfig cfg
    ) {
        Objects.requireNonNull(molecules, "molecules");
        Objects.requireNonNull(mcsProvider, "mcsProvider");
        if (cfg == null) {
            cfg = new ScaffoldDiscoveryConfig();
        }

        DescriptorHandlerLongFFP512 fingerprintHandler = DescriptorHandlerLongFFP512.getDefaultInstance();
        List<PreparedCompound> preparedCompounds = prepareCompounds(molecules, fingerprintHandler, cfg);
        List<Integer> seedIndices = selectSeedIndices(preparedCompounds.size(), cfg);

        Map<PairKey, ScaffoldTemplate> pairwiseScaffoldCache = new HashMap<>();
        Map<String, ScaffoldDecomposition> scaffoldSupportCache = new HashMap<>();
        Map<String, MutableCandidate> candidatesById = new HashMap<>();
        int pairwiseCandidateCount = 0;

        for (int seedIndex : seedIndices) {
            PreparedCompound seed = preparedCompounds.get(seedIndex);
            LinkedHashSet<Integer> localCompounds = new LinkedHashSet<>();
            localCompounds.add(seedIndex);
            localCompounds.addAll(seed.neighborIndices);

            for (int neighborIndex : seed.neighborIndices) {
                PairKey pairKey = PairKey.of(seedIndex, neighborIndex);
                ScaffoldTemplate template = pairwiseScaffoldCache.get(pairKey);
                if (!pairwiseScaffoldCache.containsKey(pairKey)) {
                    template = createPairwiseScaffold(preparedCompounds.get(pairKey.i), preparedCompounds.get(pairKey.j), mcsProvider, cfg);
                    pairwiseScaffoldCache.put(pairKey, template);
                }
                if (template == null) {
                    continue;
                }

                pairwiseCandidateCount++;
                MutableCandidate candidate = candidatesById.get(template.idcode);
                if (candidate == null) {
                    candidate = new MutableCandidate(template);
                    candidatesById.put(template.idcode, candidate);
                }
                candidate.discoverySeedIndices.add(seedIndex);
                candidate.discoveryCompoundIndices.add(seedIndex);
                candidate.discoveryCompoundIndices.add(neighborIndex);
                candidate.localSearchCompoundIndices.addAll(localCompounds);
            }
        }

        List<ScaffoldCandidate> finalCandidates = new ArrayList<>();
        ScaffoldAnalyzer.Config analyzerConfig = new ScaffoldAnalyzer.Config();
        analyzerConfig.radiusR = cfg.radiusR;
        analyzerConfig.featureMask = cfg.featureMask;

        for (MutableCandidate candidate : candidatesById.values()) {
            evaluateCandidate(candidate, preparedCompounds, analyzerConfig, scaffoldSupportCache, cfg);
            if (candidate.supportCompoundIndices.size() < cfg.minSupport) {
                continue;
            }
            finalCandidates.add(candidate.toImmutable(cfg));
        }

        finalCandidates.sort(Comparator
                .comparingDouble((ScaffoldCandidate candidate) -> candidate.combinedScore).reversed()
                .thenComparingInt(candidate -> candidate.supportCount)
                .thenComparing(candidate -> candidate.template.idcode));

        List<CompoundRecord> compoundRecords = preparedCompounds.stream()
                .map(PreparedCompound::toRecord)
                .toList();

        return new ScaffoldDiscoveryResult(
                compoundRecords,
                finalCandidates,
                seedIndices.size(),
                pairwiseCandidateCount,
                candidatesById.size()
        );
    }

    private static List<PreparedCompound> prepareCompounds(
            List<StereoMolecule> molecules,
            DescriptorHandlerLongFFP512 fingerprintHandler,
            ScaffoldDiscoveryConfig cfg
    ) {
        List<PreparedCompound> prepared = new ArrayList<>(molecules.size());
        for (int index = 0; index < molecules.size(); index++) {
            StereoMolecule molecule = new StereoMolecule(molecules.get(index));
            molecule.ensureHelperArrays(Molecule.cHelperSymmetrySimple);
            prepared.add(new PreparedCompound(
                    index,
                    molecule,
                    new Canonizer(molecule).getIDCode(),
                    fingerprintHandler.createDescriptor(molecule),
                    List.of()
            ));
        }

        List<List<Integer>> neighborIndices = computeNeighborhoods(prepared, fingerprintHandler, cfg);
        List<PreparedCompound> updated = new ArrayList<>(prepared.size());
        for (int index = 0; index < prepared.size(); index++) {
            PreparedCompound compound = prepared.get(index);
            updated.add(new PreparedCompound(
                    compound.index,
                    compound.molecule,
                    compound.idcode,
                    compound.fingerprint,
                    neighborIndices.get(index)
            ));
        }
        return updated;
    }

    private static List<List<Integer>> computeNeighborhoods(
            List<PreparedCompound> compounds,
            DescriptorHandlerLongFFP512 fingerprintHandler,
            ScaffoldDiscoveryConfig cfg
    ) {
        List<List<Integer>> neighborhoods = new ArrayList<>(compounds.size());
        for (int i = 0; i < compounds.size(); i++) {
            List<NeighborScore> scores = new ArrayList<>();
            for (int j = 0; j < compounds.size(); j++) {
                if (i == j) {
                    continue;
                }
                float similarity = fingerprintHandler.getSimilarity(compounds.get(i).fingerprint, compounds.get(j).fingerprint);
                if (similarity < cfg.minNeighborSimilarity) {
                    continue;
                }
                scores.add(new NeighborScore(j, similarity));
            }
            scores.sort(Comparator
                    .comparingDouble((NeighborScore score) -> score.similarity).reversed()
                    .thenComparingInt(score -> score.index));

            List<Integer> neighbors = new ArrayList<>(Math.min(cfg.neighborCount, scores.size()));
            for (int k = 0; k < scores.size() && k < cfg.neighborCount; k++) {
                neighbors.add(scores.get(k).index);
            }
            neighborhoods.add(List.copyOf(neighbors));
        }
        return neighborhoods;
    }

    private static List<Integer> selectSeedIndices(int moleculeCount, ScaffoldDiscoveryConfig cfg) {
        List<Integer> seeds = new ArrayList<>(moleculeCount);
        for (int index = 0; index < moleculeCount; index++) {
            seeds.add(index);
        }
        if (cfg.shuffleSeeds) {
            Collections.shuffle(seeds, new Random(cfg.randomSeed));
        }
        if (cfg.maxSeeds < seeds.size()) {
            return List.copyOf(seeds.subList(0, cfg.maxSeeds));
        }
        return List.copyOf(seeds);
    }

    private static ScaffoldTemplate createPairwiseScaffold(
            PreparedCompound compoundA,
            PreparedCompound compoundB,
            TransformationBench.StrictMCSProvider mcsProvider,
            ScaffoldDiscoveryConfig cfg
    ) {
        StereoMolecule a = new StereoMolecule(compoundA.molecule);
        StereoMolecule b = new StereoMolecule(compoundB.molecule);
        a.ensureHelperArrays(Molecule.cHelperSymmetrySimple);
        b.ensureHelperArrays(Molecule.cHelperSymmetrySimple);

        TransformationSplitter.MCSMap mcs = mcsProvider.computeStrictMCS(a, b);
        if (mcs == null) {
            return null;
        }

        BitSet coreAtoms = new BitSet(a.getAtoms());
        for (int atom = 0; atom < a.getAtoms(); atom++) {
            if (mcs.mapAtoB[atom] >= 0) {
                coreAtoms.set(atom);
            }
        }
        if (coreAtoms.isEmpty()) {
            return null;
        }

        StereoMolecule core = new StereoMolecule();
        int[] oldToCore = new int[a.getAtoms()];
        Arrays.fill(oldToCore, -1);
        a.copyMoleculeByAtoms(core, OpenChemLibUtil.bitsetToBool(coreAtoms, a.getAtoms()), true, oldToCore);
        core.setFragment(false);
        core.ensureHelperArrays(Molecule.cHelperSymmetrySimple);

        if (OpenChemLibUtil.heavyAtomCount(core) < cfg.minScaffoldHeavyAtoms) {
            return null;
        }
        if (cfg.requireConnectedScaffold && !isConnected(core)) {
            return null;
        }

        return ScaffoldTemplate.create(core);
    }

    private static void evaluateCandidate(
            MutableCandidate candidate,
            List<PreparedCompound> compounds,
            ScaffoldAnalyzer.Config analyzerConfig,
            Map<String, ScaffoldDecomposition> scaffoldSupportCache,
            ScaffoldDiscoveryConfig cfg
    ) {
        candidate.scaffoldHeavyAtomCount = OpenChemLibUtil.heavyAtomCount(candidate.template.scaffold);

        for (int compoundIndex : candidate.localSearchCompoundIndices) {
            String cacheKey = candidate.template.idcode + "|" + compoundIndex;
            ScaffoldDecomposition decomposition = scaffoldSupportCache.get(cacheKey);
            if (decomposition == null) {
                decomposition = ScaffoldAnalyzer.analyze(candidate.template, compounds.get(compoundIndex).molecule, analyzerConfig);
                scaffoldSupportCache.put(cacheKey, decomposition);
            }
            if (decomposition.failure != null) {
                continue;
            }

            candidate.supportCompoundIndices.add(compoundIndex);
            candidate.explainedFractionSum += ((double) candidate.scaffoldHeavyAtomCount)
                    / Math.max(1, OpenChemLibUtil.heavyAtomCount(compounds.get(compoundIndex).molecule));

            for (SubstitutionEvent event : decomposition.substitutionEvents) {
                candidate.observedExitVectorAtoms.addAll(event.scaffoldAtoms);
                candidate.observedExitVectorSymmetryClasses.addAll(event.symmetryClasses);
            }
        }
    }

    private static boolean isConnected(StereoMolecule molecule) {
        if (molecule.getAtoms() <= 1) {
            return true;
        }

        boolean[] visited = new boolean[molecule.getAtoms()];
        ArrayDeque<Integer> queue = new ArrayDeque<>();
        queue.add(0);
        visited[0] = true;
        int seen = 0;

        while (!queue.isEmpty()) {
            int atom = queue.removeFirst();
            seen++;
            for (int i = 0; i < molecule.getConnAtoms(atom); i++) {
                int neighbor = molecule.getConnAtom(atom, i);
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    queue.add(neighbor);
                }
            }
        }

        return seen == molecule.getAtoms();
    }

    private record PairKey(int i, int j) {
        static PairKey of(int a, int b) {
            return a <= b ? new PairKey(a, b) : new PairKey(b, a);
        }
    }

    private record NeighborScore(int index, float similarity) {}

    private static final class PreparedCompound {
        private final int index;
        private final StereoMolecule molecule;
        private final String idcode;
        private final long[] fingerprint;
        private final List<Integer> neighborIndices;

        private PreparedCompound(int index, StereoMolecule molecule, String idcode, long[] fingerprint, List<Integer> neighborIndices) {
            this.index = index;
            this.molecule = molecule;
            this.idcode = idcode;
            this.fingerprint = Arrays.copyOf(fingerprint, fingerprint.length);
            this.neighborIndices = List.copyOf(neighborIndices);
        }

        private CompoundRecord toRecord() {
            return new CompoundRecord(index, molecule, idcode, fingerprint, neighborIndices);
        }
    }

    private static final class MutableCandidate {
        private final ScaffoldTemplate template;
        private final Set<Integer> discoverySeedIndices;
        private final Set<Integer> discoveryCompoundIndices;
        private final Set<Integer> localSearchCompoundIndices;
        private final Set<Integer> supportCompoundIndices;
        private final Set<Integer> observedExitVectorAtoms;
        private final Set<Integer> observedExitVectorSymmetryClasses;
        private double explainedFractionSum;
        private int scaffoldHeavyAtomCount;

        private MutableCandidate(ScaffoldTemplate template) {
            this.template = template;
            this.discoverySeedIndices = new TreeSet<>();
            this.discoveryCompoundIndices = new TreeSet<>();
            this.localSearchCompoundIndices = new TreeSet<>();
            this.supportCompoundIndices = new TreeSet<>();
            this.observedExitVectorAtoms = new TreeSet<>();
            this.observedExitVectorSymmetryClasses = new TreeSet<>();
        }

        private ScaffoldCandidate toImmutable(ScaffoldDiscoveryConfig cfg) {
            int supportCount = supportCompoundIndices.size();
            double averageExplainedFraction = supportCount == 0 ? 0.0 : explainedFractionSum / supportCount;
            int observedExitVectorCount = observedExitVectorAtoms.size();
            double combinedScore = cfg.combinedScore(
                    supportCount,
                    averageExplainedFraction,
                    scaffoldHeavyAtomCount,
                    observedExitVectorCount
            );

            return new ScaffoldCandidate(
                    template,
                    List.copyOf(discoverySeedIndices),
                    List.copyOf(discoveryCompoundIndices),
                    List.copyOf(localSearchCompoundIndices),
                    List.copyOf(supportCompoundIndices),
                    supportCount,
                    averageExplainedFraction,
                    scaffoldHeavyAtomCount,
                    observedExitVectorCount,
                    List.copyOf(observedExitVectorAtoms),
                    List.copyOf(observedExitVectorSymmetryClasses),
                    combinedScore
            );
        }
    }
}

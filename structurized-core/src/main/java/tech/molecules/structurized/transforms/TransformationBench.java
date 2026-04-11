package tech.molecules.structurized.transforms;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

import java.time.Duration;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * TransformationBench / ComputeBench
 * ----------------------------------
 *
 * Runs a full pairwise analysis over a set of molecules. For each pair (A,B), it
 *  1) obtains a strict MCS map (atoms+bonds identical inside the preserved core),
 *  2) splits A→B into independent {@link TransformationGroup}s using {@link TransformationSplitter},
 *  3) emits per-pair results and aggregates summary statistics.
 *
 * The bench is agnostic to how the strict MCS is computed: plug your own provider via {@link StrictMCSProvider}.
 */
public class TransformationBench {

    // ======================
    // Configuration & Types
    // ======================

    public static class Config {
        public int radiusR = 1;                 // context radius in core
        public int featureMask = TransformationSplitter.FeatureMask.DEFAULT;
        public boolean keepMultiCenter = true;  // include non-(1,1) edits
        public boolean symmetricPairs = false;  // if true, analyze both (i,j) and (j,i)
        public boolean verbose = true;          // print progress
        public int maxPairs = Integer.MAX_VALUE;// soft limit for large sets
    }

    /** Provide strict MCS mapping for A and B (atoms+bonds identical in core). Return null if none. */
    public interface StrictMCSProvider {
        TransformationSplitter.MCSMap computeStrictMCS(StereoMolecule A, StereoMolecule B);
    }

    /** Aggregate statistics over all pairs. */
    public static class Summary {
        public int nPairsTried;
        public int nPairsAnalyzed;
        public int nPairsWithTG;
        public int nPairsSingleTG;
        public int nPairsMultiTG;
        public Map<TransformationType,Integer> tgTypeCounts = new EnumMap<>(TransformationType.class);
        public Map<String,Integer> signatureFrequency = new LinkedHashMap<>(); // sigId -> count
        public Duration walltime = Duration.ZERO;
    }

    // ======================
    // API
    // ======================

    /** Analyze all (unordered) pairs from molecules and return top-level {@link PairTransformation} objects. */
    public static Result run(List<StereoMolecule> molecules, StrictMCSProvider mcsProvider, Config cfg) {
        Objects.requireNonNull(molecules, "molecules");
        Objects.requireNonNull(mcsProvider, "mcsProvider");
        if (cfg == null) cfg = new Config();

        Instant t0 = Instant.now();
        Summary summary = new Summary();
        List<PairTransformation> results = new ArrayList<>();
        Map<String,Integer> sigFreq = new ConcurrentHashMap<>();
        for (TransformationType t : TransformationType.values()) summary.tgTypeCounts.put(t, 0);

        int n = molecules.size();
        int limit = Math.min(cfg.maxPairs, n * (n - 1) / (cfg.symmetricPairs ? 1 : 2));
        int processed = 0;

        outer:
        for (int i = 0; i < n; i++) {
            for (int j = (cfg.symmetricPairs ? 0 : i + 1); j < n; j++) {
                if (i == j) {
                    continue;
                }
                if (processed >= limit) {
                    break outer;
                }

                processed++;
                summary.nPairsTried++;
                StereoMolecule A = cloneWithHelpers(molecules.get(i));
                StereoMolecule B = cloneWithHelpers(molecules.get(j));

                try {
                    TransformationSplitter.MCSMap mcs = mcsProvider.computeStrictMCS(A, B);
                    if (mcs == null) {
                        results.add(new PairTransformation(i, j, "no strict MCS"));
                        continue;
                    }
                    List<TransformationGroup> groups =
                            TransformationSplitter.splitIntoTransformations(A, B, mcs, cfg.radiusR, cfg.featureMask);

                    // optionally filter out non-(1,1)
                    if (!cfg.keepMultiCenter) {
                        groups = new ArrayList<>(groups);
                        groups.removeIf(g -> g.type != TransformationType.REPLACEMENT);
                    }

                    results.add(new PairTransformation(i, j, groups));
                    summary.nPairsAnalyzed++;
                    if (!groups.isEmpty()) {
                        summary.nPairsWithTG++;
                        if (groups.size() == 1) summary.nPairsSingleTG++; else summary.nPairsMultiTG++;
                        for (var g : groups) {
                            summary.tgTypeCounts.compute(g.type, (k, v) -> v + 1);
                            sigFreq.merge(g.signature.sigId, 1, Integer::sum);
                        }
                    }

                    if (cfg.verbose && processed % 100 == 0) {
                        System.out.printf("Processed %d/%d pairs...\n", processed, limit);
                    }
                } catch (Exception ex) {
                    results.add(new PairTransformation(i, j, "exception: " + ex.getClass().getSimpleName() + ": " + ex.getMessage()));
                }
            }
        }

        summary.signatureFrequency = sortByValueDesc(sigFreq);
        summary.walltime = Duration.between(t0, Instant.now());
        return new Result(results, summary);
    }

    /** Convenience: build molecules from IDCodes. */
    public static List<StereoMolecule> fromIDCodes(List<String> idcodes) {
        IDCodeParser p = new IDCodeParser();
        List<StereoMolecule> out = new ArrayList<>(idcodes.size());
        for (String id : idcodes) {
            StereoMolecule m = new StereoMolecule();
            p.parse(m, id);
            m.ensureHelperArrays(StereoMolecule.cHelperRings);
            out.add(m);
        }
        return out;
    }

    // ======================
    // Result container
    // ======================

    public static class Result {
        public final List<PairTransformation> pairs;
        public final Summary summary;
        public Result(List<PairTransformation> pairs, Summary summary) {
            this.pairs = pairs; this.summary = summary;
        }

        public void printSummary() {
            System.out.println("=== TransformationBench Summary ===");
            System.out.printf("Pairs tried:      %d\n", summary.nPairsTried);
            System.out.printf("Pairs analyzed:   %d\n", summary.nPairsAnalyzed);
            System.out.printf("Pairs with TG:    %d (single=%d, multi=%d)\n",
                    summary.nPairsWithTG, summary.nPairsSingleTG, summary.nPairsMultiTG);
            System.out.println("TG type counts:   " + summary.tgTypeCounts);
            System.out.printf("Unique signatures: %d\n", summary.signatureFrequency.size());
            System.out.printf("Wall time: %s\n", summary.walltime);
        }
    }

    // ======================
    // Helpers
    // ======================

    private static StereoMolecule cloneWithHelpers(StereoMolecule in) {
        StereoMolecule m = new StereoMolecule(in);
        m.ensureHelperArrays(StereoMolecule.cHelperRings);
        return m;
    }

    private static LinkedHashMap<String,Integer> sortByValueDesc(Map<String,Integer> map) {
        LinkedHashMap<String,Integer> out = new LinkedHashMap<>();
        map.entrySet().stream().sorted((a,b) -> Integer.compare(b.getValue(), a.getValue()))
                .forEach(e -> out.put(e.getKey(), e.getValue()));
        return out;
    }

    // ======================
    // Example MCS providers
    // ======================

    /**
     * Example provider stub: you should plug your strict MCS here.
     * This stub returns null (no MCS). Replace with your implementation based on OCL utilities.
     */
    public static class StubStrictMCSProvider implements StrictMCSProvider {
        @Override public TransformationSplitter.MCSMap computeStrictMCS(StereoMolecule A, StereoMolecule B) {
            return null; // TODO: implement with your strict MCS algorithm
        }
    }
}

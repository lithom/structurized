package tech.molecules.structurized.scaffolds;

import java.util.List;

/**
 * Result of one scaffold discovery run.
 */
public final class ScaffoldDiscoveryResult {
    public final List<CompoundRecord> compounds;
    public final List<ScaffoldCandidate> candidates;
    public final int seedCount;
    public final int pairwiseCandidateCount;
    public final int uniqueCandidateCount;

    public ScaffoldDiscoveryResult(
            List<CompoundRecord> compounds,
            List<ScaffoldCandidate> candidates,
            int seedCount,
            int pairwiseCandidateCount,
            int uniqueCandidateCount
    ) {
        this.compounds = List.copyOf(compounds);
        this.candidates = List.copyOf(candidates);
        this.seedCount = seedCount;
        this.pairwiseCandidateCount = pairwiseCandidateCount;
        this.uniqueCandidateCount = uniqueCandidateCount;
    }
}

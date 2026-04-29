package tech.molecules.structurized.mmp;

import java.util.List;

/**
 * Complete in-memory MMP mining output.
 */
public record MmpMiningResult(
        List<MmpFragmentationRecord> fragmentationRecords,
        List<MmpPair> pairs,
        List<MmpTransformStats> transformStats
) {
    public MmpMiningResult {
        fragmentationRecords = List.copyOf(fragmentationRecords == null ? List.of() : fragmentationRecords);
        pairs = List.copyOf(pairs == null ? List.of() : pairs);
        transformStats = List.copyOf(transformStats == null ? List.of() : transformStats);
    }
}

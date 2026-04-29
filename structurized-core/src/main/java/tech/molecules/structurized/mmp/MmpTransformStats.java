package tech.molecules.structurized.mmp;

import java.util.List;
import java.util.Objects;

/**
 * Numeric delta summary for one directed MMP transformation.
 */
public record MmpTransformStats(
        String transformId,
        String fromValueIdcode,
        String toValueIdcode,
        int cutCount,
        int supportCount,
        double meanDelta,
        double medianDelta,
        double standardDeviation,
        double minDelta,
        double maxDelta,
        double positiveFraction,
        List<MmpPair> examplePairs
) {
    public MmpTransformStats {
        Objects.requireNonNull(transformId, "transformId");
        Objects.requireNonNull(fromValueIdcode, "fromValueIdcode");
        Objects.requireNonNull(toValueIdcode, "toValueIdcode");
        examplePairs = List.copyOf(examplePairs == null ? List.of() : examplePairs);
    }
}

package tech.molecules.structurized.mmp;

import java.util.Objects;

/**
 * Directed matched molecular pair sharing the same constant key.
 */
public record MmpPair(
        String compoundIdA,
        String compoundIdB,
        Double valueA,
        Double valueB,
        Double delta,
        String keyIdcode,
        String fromValueIdcode,
        String toValueIdcode,
        String transformId,
        int cutCount
) {
    public MmpPair {
        if (compoundIdA == null || compoundIdA.isBlank()) {
            throw new IllegalArgumentException("compoundIdA must not be blank");
        }
        if (compoundIdB == null || compoundIdB.isBlank()) {
            throw new IllegalArgumentException("compoundIdB must not be blank");
        }
        Objects.requireNonNull(keyIdcode, "keyIdcode");
        Objects.requireNonNull(fromValueIdcode, "fromValueIdcode");
        Objects.requireNonNull(toValueIdcode, "toValueIdcode");
        if (cutCount < 1) {
            throw new IllegalArgumentException("cutCount must be positive");
        }
        if (transformId == null || transformId.isBlank()) {
            transformId = createTransformId(cutCount, fromValueIdcode, toValueIdcode);
        }
        if (delta == null && valueA != null && valueB != null) {
            delta = valueB - valueA;
        }
    }

    public static String createTransformId(int cutCount, String fromValueIdcode, String toValueIdcode) {
        return cutCount + "|" + fromValueIdcode + ">>" + toValueIdcode;
    }
}

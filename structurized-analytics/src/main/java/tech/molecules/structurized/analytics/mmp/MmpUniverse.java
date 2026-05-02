package tech.molecules.structurized.analytics.mmp;

import java.time.Instant;
import java.util.List;
import java.util.Objects;

/**
 * Named compound universe used for one structural MMP computation.
 */
public record MmpUniverse(
        String universeId,
        String name,
        List<String> subjectSetIds,
        String mmpConfigHash,
        Instant createdAt,
        String metadata
) {
    public MmpUniverse {
        universeId = requireText(universeId, "universeId");
        name = requireText(name, "name");
        subjectSetIds = List.copyOf(subjectSetIds == null ? List.of() : subjectSetIds);
        mmpConfigHash = requireText(mmpConfigHash, "mmpConfigHash");
        createdAt = Objects.requireNonNull(createdAt, "createdAt");
        metadata = normalize(metadata);
    }

    private static String requireText(String value, String field) {
        String normalized = normalize(value);
        if (normalized == null) {
            throw new IllegalArgumentException(field + " must not be blank");
        }
        return normalized;
    }

    private static String normalize(String value) {
        if (value == null) {
            return null;
        }
        String trimmed = value.trim();
        return trimmed.isEmpty() ? null : trimmed;
    }
}

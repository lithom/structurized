package tech.molecules.structurized.analytics.mmp;

import java.time.Instant;
import java.util.Objects;

/**
 * Metadata for one endpoint-specific MMP statistics computation.
 */
public record MmpEndpointStatsRun(
        String runId,
        String endpointId,
        String endpointSubjectSetId,
        String universeId,
        String mmpConfigHash,
        String statsConfigHash,
        Instant createdAt,
        int subjectCount,
        int valueCount,
        int pairCount,
        int statsCount,
        String metadata
) {
    public MmpEndpointStatsRun {
        runId = requireText(runId, "runId");
        endpointId = requireText(endpointId, "endpointId");
        endpointSubjectSetId = requireText(endpointSubjectSetId, "endpointSubjectSetId");
        universeId = requireText(universeId, "universeId");
        mmpConfigHash = requireText(mmpConfigHash, "mmpConfigHash");
        statsConfigHash = requireText(statsConfigHash, "statsConfigHash");
        createdAt = Objects.requireNonNull(createdAt, "createdAt");
        if (subjectCount < 0 || valueCount < 0 || pairCount < 0 || statsCount < 0) {
            throw new IllegalArgumentException("run counts must not be negative");
        }
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

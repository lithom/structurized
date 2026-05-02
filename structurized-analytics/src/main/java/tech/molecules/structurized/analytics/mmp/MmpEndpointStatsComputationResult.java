package tech.molecules.structurized.analytics.mmp;

import java.util.List;

/**
 * Summary returned by the PRISM-backed MMP statistics computation service.
 */
public record MmpEndpointStatsComputationResult(
        List<MmpUniverse> universes,
        List<MmpEndpointStatsRun> statsRuns,
        int requestedEndpointCount,
        int structuralSubjectCount,
        int missingStructureCount,
        int fragmentationRecordCount,
        int pairCount,
        List<String> warnings
) {
    public MmpEndpointStatsComputationResult {
        universes = List.copyOf(universes == null ? List.of() : universes);
        statsRuns = List.copyOf(statsRuns == null ? List.of() : statsRuns);
        warnings = List.copyOf(warnings == null ? List.of() : warnings);
    }
}

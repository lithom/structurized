package tech.molecules.structurized.analytics.mmp;

import com.actelion.research.chem.StereoMolecule;
import tech.molecules.structurized.mmp.MmpInputCompound;
import tech.molecules.structurized.mmp.MmpMiner;
import tech.molecules.structurized.mmp.MmpMiningConfig;
import tech.molecules.structurized.mmp.MmpMiningResult;
import tech.molecules.structurized.mmp.MmpPair;
import tech.molecules.structurized.mmp.MmpStatsAggregator;
import tech.molecules.structurized.mmp.MmpTransformStats;
import tech.molecules.structurized.prism.model.EndpointDataType;
import tech.molecules.structurized.prism.model.EndpointDefinition;
import tech.molecules.structurized.prism.provider.EndpointProvider;
import tech.molecules.structurized.prism.provider.SubjectSet;
import tech.molecules.structurized.prism.provider.SubjectSetProvider;
import tech.molecules.structurized.prism.query.EndpointFetchRequest;
import tech.molecules.structurized.prism.query.EndpointValueRecord;
import tech.molecules.structurized.prism.result.NumericResult;
import tech.molecules.structurized.prism.result.NumericState;
import tech.molecules.structurized.prism.result.OptionalNumericResult;
import tech.molecules.structurized.prism.result.OptionalNumericState;

import java.time.Instant;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;

/**
 * Computes MMP transform statistics for numeric PRISM endpoints and persists the results.
 */
public final class MmpEndpointStatsComputationService {
    private final EndpointProvider endpointProvider;
    private final SubjectSetProvider subjectSetProvider;
    private final StructureProvider structureProvider;
    private final MmpUniverseRepository universeRepository;
    private final MmpPairRepository pairRepository;
    private final MmpEndpointStatsRepository statsRepository;

    public MmpEndpointStatsComputationService(
            EndpointProvider endpointProvider,
            SubjectSetProvider subjectSetProvider,
            StructureProvider structureProvider,
            MmpUniverseRepository universeRepository,
            MmpPairRepository pairRepository,
            MmpEndpointStatsRepository statsRepository
    ) {
        this.endpointProvider = Objects.requireNonNull(endpointProvider, "endpointProvider");
        this.subjectSetProvider = Objects.requireNonNull(subjectSetProvider, "subjectSetProvider");
        this.structureProvider = Objects.requireNonNull(structureProvider, "structureProvider");
        this.universeRepository = Objects.requireNonNull(universeRepository, "universeRepository");
        this.pairRepository = Objects.requireNonNull(pairRepository, "pairRepository");
        this.statsRepository = Objects.requireNonNull(statsRepository, "statsRepository");
    }

    public MmpEndpointStatsComputationResult computeAndPersist(
            MmpEndpointStatsConfig statsConfig,
            MmpMiningConfig mmpConfig
    ) {
        Objects.requireNonNull(statsConfig, "statsConfig");
        Objects.requireNonNull(mmpConfig, "mmpConfig");

        List<EndpointDefinition> endpoints = selectNumericEndpoints(statsConfig);
        Map<String, EndpointSubjectSet> endpointSets = resolveEndpointSubjectSets(endpoints, statsConfig);
        String mmpConfigHash = MmpAnalyticsHashes.mmpConfigHash(mmpConfig);
        String statsConfigHash = MmpAnalyticsHashes.statsConfigHash(statsConfig);
        ArrayList<String> warnings = new ArrayList<>();

        if (statsConfig.universeMode() == MmpUniverseMode.PER_ENDPOINT) {
            return computePerEndpoint(endpoints, endpointSets, statsConfig, mmpConfig, mmpConfigHash, statsConfigHash, warnings);
        }
        return computeUnion(endpoints, endpointSets, statsConfig, mmpConfig, mmpConfigHash, statsConfigHash, warnings);
    }

    private MmpEndpointStatsComputationResult computeUnion(
            List<EndpointDefinition> endpoints,
            Map<String, EndpointSubjectSet> endpointSets,
            MmpEndpointStatsConfig statsConfig,
            MmpMiningConfig mmpConfig,
            String mmpConfigHash,
            String statsConfigHash,
            List<String> warnings
    ) {
        LinkedHashSet<String> unionSubjects = new LinkedHashSet<>();
        ArrayList<String> subjectSetIds = new ArrayList<>();
        for (EndpointDefinition endpoint : endpoints) {
            EndpointSubjectSet endpointSet = endpointSets.get(endpoint.getId());
            subjectSetIds.add(endpointSet.subjectSetId());
            unionSubjects.addAll(endpointSet.subjectIds());
        }

        ComputedUniverse computed = computeUniverse(
                statsConfig.universeId() != null
                        ? statsConfig.universeId()
                        : "mmp-union-" + shortHash(mmpConfigHash + "|" + String.join(",", subjectSetIds)),
                statsConfig.universeName() != null ? statsConfig.universeName() : "Union MMP universe",
                subjectSetIds,
                List.copyOf(unionSubjects),
                mmpConfig,
                mmpConfigHash,
                warnings
        );
        persistUniverse(computed);

        ArrayList<MmpEndpointStatsRun> runs = new ArrayList<>();
        for (EndpointDefinition endpoint : endpoints) {
            EndpointSubjectSet endpointSet = endpointSets.get(endpoint.getId());
            runs.add(computeAndPersistStatsRun(
                    endpoint,
                    endpointSet,
                    computed.universe(),
                    computed.result().pairs(),
                    statsConfig,
                    mmpConfig,
                    mmpConfigHash,
                    statsConfigHash
            ));
        }
        return new MmpEndpointStatsComputationResult(
                List.of(computed.universe()),
                runs,
                endpoints.size(),
                computed.structuralSubjectCount(),
                computed.missingStructureCount(),
                computed.result().fragmentationRecords().size(),
                computed.result().pairs().size(),
                warnings
        );
    }

    private MmpEndpointStatsComputationResult computePerEndpoint(
            List<EndpointDefinition> endpoints,
            Map<String, EndpointSubjectSet> endpointSets,
            MmpEndpointStatsConfig statsConfig,
            MmpMiningConfig mmpConfig,
            String mmpConfigHash,
            String statsConfigHash,
            List<String> warnings
    ) {
        ArrayList<MmpUniverse> universes = new ArrayList<>();
        ArrayList<MmpEndpointStatsRun> runs = new ArrayList<>();
        int structuralSubjects = 0;
        int missingStructures = 0;
        int fragmentationRecords = 0;
        int pairs = 0;

        for (EndpointDefinition endpoint : endpoints) {
            EndpointSubjectSet endpointSet = endpointSets.get(endpoint.getId());
            ComputedUniverse computed = computeUniverse(
                    "mmp-" + endpoint.getId() + "-" + shortHash(mmpConfigHash + "|" + endpointSet.subjectSetId()),
                    "MMP universe for " + endpoint.getName(),
                    List.of(endpointSet.subjectSetId()),
                    endpointSet.subjectIds(),
                    mmpConfig,
                    mmpConfigHash,
                    warnings
            );
            persistUniverse(computed);
            universes.add(computed.universe());
            structuralSubjects += computed.structuralSubjectCount();
            missingStructures += computed.missingStructureCount();
            fragmentationRecords += computed.result().fragmentationRecords().size();
            pairs += computed.result().pairs().size();
            runs.add(computeAndPersistStatsRun(
                    endpoint,
                    endpointSet,
                    computed.universe(),
                    computed.result().pairs(),
                    statsConfig,
                    mmpConfig,
                    mmpConfigHash,
                    statsConfigHash
            ));
        }

        return new MmpEndpointStatsComputationResult(
                universes,
                runs,
                endpoints.size(),
                structuralSubjects,
                missingStructures,
                fragmentationRecords,
                pairs,
                warnings
        );
    }

    private ComputedUniverse computeUniverse(
            String universeId,
            String universeName,
            List<String> subjectSetIds,
            List<String> subjectIds,
            MmpMiningConfig mmpConfig,
            String mmpConfigHash,
            List<String> warnings
    ) {
        Map<String, StereoMolecule> structures = structureProvider.fetchStructures(subjectIds);
        int missing = subjectIds.size() - structures.size();
        if (missing > 0) {
            warnings.add("Skipped " + missing + " subjects without structures for universe " + universeId);
        }

        ArrayList<MmpInputCompound> compounds = new ArrayList<>();
        structures.entrySet().stream()
                .sorted(Map.Entry.comparingByKey())
                .forEach(entry -> compounds.add(new MmpInputCompound(entry.getKey(), entry.getValue(), null)));
        MmpMiningResult result = MmpMiner.mine(compounds, mmpConfig);
        MmpUniverse universe = new MmpUniverse(
                universeId,
                universeName,
                subjectSetIds,
                mmpConfigHash,
                Instant.now(),
                "subjects=" + subjectIds.size() + ";structures=" + structures.size() + ";missingStructures=" + missing
        );
        return new ComputedUniverse(universe, result, structures.keySet().stream().sorted().toList(), missing);
    }

    private void persistUniverse(ComputedUniverse computed) {
        universeRepository.saveUniverse(computed.universe(), computed.structuralSubjectIds());
        pairRepository.replaceFragmentationRecords(computed.universe().universeId(), computed.result().fragmentationRecords());
        pairRepository.replacePairs(computed.universe().universeId(), computed.result().pairs());
    }

    private MmpEndpointStatsRun computeAndPersistStatsRun(
            EndpointDefinition endpoint,
            EndpointSubjectSet endpointSet,
            MmpUniverse universe,
            List<MmpPair> structuralPairs,
            MmpEndpointStatsConfig statsConfig,
            MmpMiningConfig mmpConfig,
            String mmpConfigHash,
            String statsConfigHash
    ) {
        Map<String, Double> values = fetchEndpointValues(endpoint.getId(), endpointSet.subjectIds(), statsConfig.batchSize());
        ArrayList<MmpPair> valuedPairs = new ArrayList<>();
        for (MmpPair pair : structuralPairs) {
            Double valueA = values.get(pair.compoundIdA());
            Double valueB = values.get(pair.compoundIdB());
            if (valueA == null || valueB == null) {
                continue;
            }
            valuedPairs.add(new MmpPair(
                    pair.compoundIdA(),
                    pair.compoundIdB(),
                    valueA,
                    valueB,
                    null,
                    pair.keyIdcode(),
                    pair.fromValueIdcode(),
                    pair.toValueIdcode(),
                    pair.transformId(),
                    pair.cutCount()
            ));
        }
        List<MmpTransformStats> stats = MmpStatsAggregator.aggregate(valuedPairs, mmpConfig);
        MmpEndpointStatsRun run = new MmpEndpointStatsRun(
                "mmp-stats-" + endpoint.getId() + "-" + shortHash(universe.universeId() + "|" + statsConfigHash + "|" + Instant.now()),
                endpoint.getId(),
                endpointSet.subjectSetId(),
                universe.universeId(),
                mmpConfigHash,
                statsConfigHash,
                Instant.now(),
                endpointSet.subjectIds().size(),
                values.size(),
                valuedPairs.size(),
                stats.size(),
                "endpointName=" + endpoint.getName()
        );
        statsRepository.saveStatsRun(run, stats);
        return run;
    }

    private List<EndpointDefinition> selectNumericEndpoints(MmpEndpointStatsConfig config) {
        LinkedHashMap<String, EndpointDefinition> byId = new LinkedHashMap<>();
        for (EndpointDefinition endpoint : endpointProvider.listEndpointDefinitions()) {
            if (isNumericEndpoint(endpoint)) {
                byId.put(endpoint.getId(), endpoint);
            }
        }
        if (config.endpointIds().isEmpty()) {
            return byId.values().stream()
                    .sorted(Comparator.comparing(EndpointDefinition::getId))
                    .toList();
        }

        ArrayList<EndpointDefinition> selected = new ArrayList<>();
        for (String endpointId : config.endpointIds()) {
            EndpointDefinition endpoint = byId.get(endpointId);
            if (endpoint == null) {
                throw new IllegalArgumentException("unknown or non-numeric endpoint '" + endpointId + "'");
            }
            selected.add(endpoint);
        }
        return List.copyOf(selected);
    }

    private static boolean isNumericEndpoint(EndpointDefinition endpoint) {
        return endpoint.getDatatype() == EndpointDataType.NUMERIC
                || endpoint.getDatatype() == EndpointDataType.OPTIONAL_NUMERIC;
    }

    private Map<String, EndpointSubjectSet> resolveEndpointSubjectSets(
            List<EndpointDefinition> endpoints,
            MmpEndpointStatsConfig config
    ) {
        List<SubjectSet> allSets = subjectSetProvider.listSubjectSets();
        LinkedHashMap<String, EndpointSubjectSet> result = new LinkedHashMap<>();
        for (EndpointDefinition endpoint : endpoints) {
            String subjectSetId = config.endpointSubjectSetIds().get(endpoint.getId());
            if (subjectSetId == null) {
                subjectSetId = inferSubjectSetId(endpoint, allSets, config)
                        .orElseThrow(() -> new IllegalArgumentException(
                                "no measured subject set configured or inferred for endpoint '" + endpoint.getId() + "'"));
            }
            result.put(endpoint.getId(), new EndpointSubjectSet(subjectSetId, fetchAllSubjects(subjectSetId, config.batchSize())));
        }
        return Map.copyOf(result);
    }

    private Optional<String> inferSubjectSetId(
            EndpointDefinition endpoint,
            List<SubjectSet> allSets,
            MmpEndpointStatsConfig config
    ) {
        List<String> candidates = List.of(
                endpoint.getId(),
                "assay:" + endpoint.getId() + ":measured",
                "assay-measured:" + endpoint.getId(),
                "endpoint:" + endpoint.getId() + ":measured",
                endpoint.getPath(),
                "assay:" + endpoint.getPath() + ":measured"
        );
        for (String candidate : candidates) {
            for (SubjectSet subjectSet : allSets) {
                if (candidate.equals(subjectSet.getId())) {
                    return Optional.of(candidate);
                }
            }
        }

        return allSets.stream()
                .filter(set -> config.measuredSubjectSetType() == null
                        || config.measuredSubjectSetType().equals(set.getSetType()))
                .filter(set -> config.measuredSubjectSetScope() == null
                        || config.measuredSubjectSetScope().equals(set.getSubjectSetScope()))
                .filter(set -> set.getId().contains(endpoint.getId()) || set.getName().contains(endpoint.getName()))
                .map(SubjectSet::getId)
                .findFirst();
    }

    private List<String> fetchAllSubjects(String subjectSetId, int batchSize) {
        long count = subjectSetProvider.countSubjects(subjectSetId);
        ArrayList<String> subjects = new ArrayList<>();
        for (int offset = 0; offset < count; offset += batchSize) {
            subjects.addAll(subjectSetProvider.listSubjects(subjectSetId, offset, batchSize));
        }
        return subjects.stream().distinct().toList();
    }

    private Map<String, Double> fetchEndpointValues(String endpointId, List<String> subjectIds, int batchSize) {
        LinkedHashMap<String, Double> values = new LinkedHashMap<>();
        for (int offset = 0; offset < subjectIds.size(); offset += batchSize) {
            List<String> batch = subjectIds.subList(offset, Math.min(offset + batchSize, subjectIds.size()));
            EndpointFetchRequest request = EndpointFetchRequest.builder()
                    .subjectIds(batch)
                    .endpointIds(List.of(endpointId))
                    .build();
            for (EndpointValueRecord record : endpointProvider.fetchEndpointValues(request)) {
                extractNumericValue(record).ifPresent(value -> values.put(record.getSubjectId(), value));
            }
        }
        return Map.copyOf(values);
    }

    private static Optional<Double> extractNumericValue(EndpointValueRecord record) {
        if (record.getResult() instanceof NumericResult numeric
                && numeric.getState() == NumericState.VALUE
                && numeric.getMean() != null) {
            return Optional.of(numeric.getMean());
        }
        if (record.getResult() instanceof OptionalNumericResult optional
                && optional.getState() == OptionalNumericState.VALUE
                && optional.getMean() != null) {
            return Optional.of(optional.getMean());
        }
        return Optional.empty();
    }

    private static String shortHash(String value) {
        return Integer.toUnsignedString(value.hashCode(), 36);
    }

    private record EndpointSubjectSet(String subjectSetId, List<String> subjectIds) {
        private EndpointSubjectSet {
            subjectIds = List.copyOf(subjectIds);
        }
    }

    private record ComputedUniverse(
            MmpUniverse universe,
            MmpMiningResult result,
            List<String> structuralSubjectIds,
            int missingStructureCount
    ) {
        int structuralSubjectCount() {
            return structuralSubjectIds.size();
        }
    }
}

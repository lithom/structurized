package tech.molecules.structurized.analytics.mmp;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Configuration for PRISM-backed MMP statistics computation.
 */
public final class MmpEndpointStatsConfig {
    private final List<String> endpointIds;
    private final Map<String, String> endpointSubjectSetIds;
    private final MmpUniverseMode universeMode;
    private final String universeId;
    private final String universeName;
    private final String measuredSubjectSetType;
    private final String measuredSubjectSetScope;
    private final int batchSize;

    private MmpEndpointStatsConfig(Builder builder) {
        this.endpointIds = List.copyOf(new ArrayList<>(builder.endpointIds));
        this.endpointSubjectSetIds = Map.copyOf(new LinkedHashMap<>(builder.endpointSubjectSetIds));
        this.universeMode = Objects.requireNonNull(builder.universeMode, "universeMode");
        this.universeId = normalize(builder.universeId);
        this.universeName = normalize(builder.universeName);
        this.measuredSubjectSetType = normalize(builder.measuredSubjectSetType);
        this.measuredSubjectSetScope = normalize(builder.measuredSubjectSetScope);
        this.batchSize = builder.batchSize;
        if (batchSize < 1) {
            throw new IllegalArgumentException("batchSize must be positive");
        }
    }

    public static Builder builder() {
        return new Builder();
    }

    public List<String> endpointIds() {
        return endpointIds;
    }

    public Map<String, String> endpointSubjectSetIds() {
        return endpointSubjectSetIds;
    }

    public MmpUniverseMode universeMode() {
        return universeMode;
    }

    public String universeId() {
        return universeId;
    }

    public String universeName() {
        return universeName;
    }

    public String measuredSubjectSetType() {
        return measuredSubjectSetType;
    }

    public String measuredSubjectSetScope() {
        return measuredSubjectSetScope;
    }

    public int batchSize() {
        return batchSize;
    }

    public static final class Builder {
        private List<String> endpointIds = List.of();
        private Map<String, String> endpointSubjectSetIds = Map.of();
        private MmpUniverseMode universeMode = MmpUniverseMode.UNION_OF_ENDPOINT_SUBJECTS;
        private String universeId;
        private String universeName;
        private String measuredSubjectSetType = "ASSAY_MEASURED";
        private String measuredSubjectSetScope = "ASSAYS";
        private int batchSize = 500;

        private Builder() {}

        public Builder endpointIds(List<String> endpointIds) {
            this.endpointIds = endpointIds == null ? List.of() : List.copyOf(endpointIds);
            return this;
        }

        public Builder addEndpointId(String endpointId) {
            Objects.requireNonNull(endpointId, "endpointId");
            ArrayList<String> next = new ArrayList<>(endpointIds);
            next.add(endpointId);
            endpointIds = List.copyOf(next);
            return this;
        }

        public Builder endpointSubjectSetIds(Map<String, String> endpointSubjectSetIds) {
            this.endpointSubjectSetIds = endpointSubjectSetIds == null ? Map.of() : Map.copyOf(endpointSubjectSetIds);
            return this;
        }

        public Builder putEndpointSubjectSetId(String endpointId, String subjectSetId) {
            LinkedHashMap<String, String> next = new LinkedHashMap<>(endpointSubjectSetIds);
            next.put(requireText(endpointId, "endpointId"), requireText(subjectSetId, "subjectSetId"));
            endpointSubjectSetIds = Map.copyOf(next);
            return this;
        }

        public Builder universeMode(MmpUniverseMode universeMode) {
            this.universeMode = universeMode;
            return this;
        }

        public Builder universeId(String universeId) {
            this.universeId = universeId;
            return this;
        }

        public Builder universeName(String universeName) {
            this.universeName = universeName;
            return this;
        }

        public Builder measuredSubjectSetType(String measuredSubjectSetType) {
            this.measuredSubjectSetType = measuredSubjectSetType;
            return this;
        }

        public Builder measuredSubjectSetScope(String measuredSubjectSetScope) {
            this.measuredSubjectSetScope = measuredSubjectSetScope;
            return this;
        }

        public Builder batchSize(int batchSize) {
            this.batchSize = batchSize;
            return this;
        }

        public MmpEndpointStatsConfig build() {
            return new MmpEndpointStatsConfig(this);
        }
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

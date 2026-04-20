package tech.molecules.structurized.scaffolds;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.IsomericSmilesCreator;
import com.actelion.research.chem.StereoMolecule;
import tech.molecules.structurized.OpenChemLibUtil;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Full-dataset scaffold-relative decomposition plus simple 1D/2D projection views.
 */
public final class ScaffoldDatasetDecomposition {
    public final ScaffoldTemplate template;
    public final List<CompoundDecompositionRecord> records;
    public final List<Integer> observedExitVectorAtoms;
    public final List<Integer> observedExitVectorSymmetryClasses;
    public final int matchedCompoundCount;
    public final int unmatchedCompoundCount;
    public final int multiAttachmentCompoundCount;

    private ScaffoldDatasetDecomposition(
            ScaffoldTemplate template,
            List<CompoundDecompositionRecord> records,
            List<Integer> observedExitVectorAtoms,
            List<Integer> observedExitVectorSymmetryClasses,
            int matchedCompoundCount,
            int unmatchedCompoundCount,
            int multiAttachmentCompoundCount
    ) {
        this.template = template;
        this.records = List.copyOf(records);
        this.observedExitVectorAtoms = List.copyOf(observedExitVectorAtoms);
        this.observedExitVectorSymmetryClasses = List.copyOf(observedExitVectorSymmetryClasses);
        this.matchedCompoundCount = matchedCompoundCount;
        this.unmatchedCompoundCount = unmatchedCompoundCount;
        this.multiAttachmentCompoundCount = multiAttachmentCompoundCount;
    }

    public static ScaffoldDatasetDecomposition analyze(
            List<CompoundRecord> compounds,
            ScaffoldTemplate template,
            ScaffoldAnalyzer.Config cfg
    ) {
        Objects.requireNonNull(compounds, "compounds");
        Objects.requireNonNull(template, "template");
        if (cfg == null) {
            cfg = new ScaffoldAnalyzer.Config();
        }

        int scaffoldHeavyAtomCount = OpenChemLibUtil.heavyAtomCount(template.scaffold);
        Map<String, String> fragmentSmilesCache = new TreeMap<>();
        List<CompoundDecompositionRecord> records = new ArrayList<>(compounds.size());
        Set<Integer> observedAtoms = new TreeSet<>();
        Set<Integer> multiAttachmentCompounds = new TreeSet<>();
        int matchedCount = 0;

        for (CompoundRecord compound : compounds) {
            ScaffoldDecomposition decomposition = ScaffoldAnalyzer.analyze(template, compound.molecule, cfg);
            boolean matched = decomposition.failure == null;
            double explainedFraction = matched
                    ? ((double) scaffoldHeavyAtomCount) / Math.max(1, OpenChemLibUtil.heavyAtomCount(compound.molecule))
                    : 0.0;

            List<Integer> occupiedAtoms = new ArrayList<>();
            List<Integer> occupiedSymmetryClasses = new ArrayList<>();
            Map<Integer, ExitVectorAssignment> assignments = new LinkedHashMap<>();
            Set<Integer> multiAttachmentAtoms = new TreeSet<>();
            Set<Integer> ambiguousAtoms = new TreeSet<>();
            boolean hasMultiAttachment = false;

            if (matched) {
                matchedCount++;
                for (SubstitutionEvent event : decomposition.substitutionEvents) {
                    occupiedAtoms.addAll(event.scaffoldAtoms);
                    occupiedSymmetryClasses.addAll(event.symmetryClasses);
                    observedAtoms.addAll(event.scaffoldAtoms);

                    if (event.eventType == SubstitutionEventType.SINGLE_ATTACHMENT && event.scaffoldAtoms.size() == 1) {
                        int scaffoldAtom = event.scaffoldAtoms.getFirst();
                        ExitVectorAssignment assignment = new ExitVectorAssignment(
                                scaffoldAtom,
                                event.symmetryClasses.getFirst(),
                                event.addedFragmentIdcode,
                                fragmentSmiles(fragmentSmilesCache, event.addedFragmentIdcode)
                        );

                        ExitVectorAssignment existing = assignments.get(scaffoldAtom);
                        if (existing != null && !Objects.equals(existing.fragmentIdcode, assignment.fragmentIdcode)) {
                            assignments.remove(scaffoldAtom);
                            ambiguousAtoms.add(scaffoldAtom);
                            continue;
                        }
                        if (!ambiguousAtoms.contains(scaffoldAtom)) {
                            assignments.put(scaffoldAtom, assignment);
                        }
                    } else {
                        hasMultiAttachment = true;
                        multiAttachmentAtoms.addAll(event.scaffoldAtoms);
                    }
                }
            }

            if (hasMultiAttachment) {
                multiAttachmentCompounds.add(compound.index);
            }

            records.add(new CompoundDecompositionRecord(
                    compound,
                    decomposition,
                    matched,
                    explainedFraction,
                    occupiedAtoms.stream().distinct().sorted().toList(),
                    occupiedSymmetryClasses.stream().distinct().sorted().toList(),
                    assignments,
                    multiAttachmentAtoms,
                    ambiguousAtoms,
                    hasMultiAttachment
            ));
        }

        List<Integer> observedExitVectorAtoms = List.copyOf(observedAtoms);
        List<Integer> observedExitVectorSymmetryClasses = observedExitVectorAtoms.stream()
                .map(atom -> template.atomSymmetryClasses[atom])
                .toList();

        return new ScaffoldDatasetDecomposition(
                template,
                records,
                observedExitVectorAtoms,
                observedExitVectorSymmetryClasses,
                matchedCount,
                compounds.size() - matchedCount,
                multiAttachmentCompounds.size()
        );
    }

    public String exitVectorLabel(int scaffoldAtom) {
        return "Atom " + (scaffoldAtom + 1) + " (sym " + template.atomSymmetryClasses[scaffoldAtom] + ")";
    }

    public OneDimProjection oneDimProjection(int scaffoldAtom) {
        return oneDimProjection(scaffoldAtom, true);
    }

    public OneDimProjection oneDimProjection(int scaffoldAtom, boolean includeUnmatched) {
        Map<ProjectionBucket, List<Integer>> indicesByBucket = new LinkedHashMap<>();
        for (CompoundDecompositionRecord record : records) {
            if (!includeUnmatched && !record.matched) {
                continue;
            }
            ProjectionBucket bucket = bucketFor(record, scaffoldAtom);
            indicesByBucket.computeIfAbsent(bucket, ignored -> new ArrayList<>()).add(record.compound.index);
        }

        List<OneDimProjectionRow> rows = indicesByBucket.entrySet().stream()
                .map(entry -> new OneDimProjectionRow(entry.getKey(), List.copyOf(entry.getValue())))
                .sorted(Comparator
                        .comparingInt((OneDimProjectionRow row) -> row.compoundIndices.size()).reversed()
                        .thenComparing(row -> row.bucket.sortKey))
                .toList();

        return new OneDimProjection(scaffoldAtom, exitVectorLabel(scaffoldAtom), rows);
    }

    public TwoDimProjection twoDimProjection(int rowScaffoldAtom, int columnScaffoldAtom) {
        return twoDimProjection(rowScaffoldAtom, columnScaffoldAtom, true);
    }

    public TwoDimProjection twoDimProjection(int rowScaffoldAtom, int columnScaffoldAtom, boolean includeUnmatched) {
        Map<ProjectionBucket, Integer> rowOrder = new LinkedHashMap<>();
        Map<ProjectionBucket, Integer> columnOrder = new LinkedHashMap<>();
        Map<CellKey, List<Integer>> compoundsByCell = new LinkedHashMap<>();

        for (CompoundDecompositionRecord record : records) {
            if (!includeUnmatched && !record.matched) {
                continue;
            }
            ProjectionBucket rowBucket = bucketFor(record, rowScaffoldAtom);
            ProjectionBucket columnBucket = bucketFor(record, columnScaffoldAtom);
            rowOrder.computeIfAbsent(rowBucket, ignored -> rowOrder.size());
            columnOrder.computeIfAbsent(columnBucket, ignored -> columnOrder.size());
        }

        List<ProjectionBucket> rowBuckets = sortBuckets(rowOrder.keySet());
        List<ProjectionBucket> columnBuckets = sortBuckets(columnOrder.keySet());
        Map<ProjectionBucket, Integer> rowIndex = indexBuckets(rowBuckets);
        Map<ProjectionBucket, Integer> columnIndex = indexBuckets(columnBuckets);
        int[][] counts = new int[rowBuckets.size()][columnBuckets.size()];

        for (CompoundDecompositionRecord record : records) {
            if (!includeUnmatched && !record.matched) {
                continue;
            }
            ProjectionBucket rowBucket = bucketFor(record, rowScaffoldAtom);
            ProjectionBucket columnBucket = bucketFor(record, columnScaffoldAtom);
            int r = rowIndex.get(rowBucket);
            int c = columnIndex.get(columnBucket);
            counts[r][c]++;
            compoundsByCell.computeIfAbsent(new CellKey(r, c), ignored -> new ArrayList<>()).add(record.compound.index);
        }

        return new TwoDimProjection(
                rowScaffoldAtom,
                columnScaffoldAtom,
                exitVectorLabel(rowScaffoldAtom),
                exitVectorLabel(columnScaffoldAtom),
                rowBuckets,
                columnBuckets,
                counts,
                compoundsByCell.entrySet().stream()
                        .collect(LinkedHashMap::new,
                                (map, entry) -> map.put(entry.getKey(), List.copyOf(entry.getValue())),
                                Map::putAll)
        );
    }

    private static Map<ProjectionBucket, Integer> indexBuckets(List<ProjectionBucket> buckets) {
        Map<ProjectionBucket, Integer> index = new LinkedHashMap<>();
        for (int i = 0; i < buckets.size(); i++) {
            index.put(buckets.get(i), i);
        }
        return index;
    }

    private static List<ProjectionBucket> sortBuckets(Set<ProjectionBucket> buckets) {
        return buckets.stream()
                .sorted(Comparator.comparing(bucket -> bucket.sortKey))
                .toList();
    }

    private ProjectionBucket bucketFor(CompoundDecompositionRecord record, int scaffoldAtom) {
        if (!record.matched) {
            return ProjectionBucket.unmatched();
        }
        if (record.multiAttachmentAtoms.contains(scaffoldAtom) || record.ambiguousAtoms.contains(scaffoldAtom)) {
            return ProjectionBucket.multiAttachment();
        }
        ExitVectorAssignment assignment = record.assignmentsByScaffoldAtom.get(scaffoldAtom);
        if (assignment != null) {
            return ProjectionBucket.substituent(assignment.fragmentIdcode, assignment.fragmentSmiles);
        }
        return ProjectionBucket.unsubstituted();
    }

    private static String fragmentSmiles(Map<String, String> cache, String fragmentIdcode) {
        if (fragmentIdcode == null || fragmentIdcode.isBlank()) {
            return "";
        }
        return cache.computeIfAbsent(fragmentIdcode, idcode -> {
            StereoMolecule fragment = new StereoMolecule();
            new IDCodeParser().parse(fragment, idcode);
            return new IsomericSmilesCreator(fragment).getSmiles();
        });
    }

    public enum ProjectionBucketType {
        SUBSTITUENT,
        UNSUBSTITUTED,
        MULTI_ATTACHMENT,
        UNMATCHED
    }

    public static final class ProjectionBucket {
        public final ProjectionBucketType type;
        public final String fragmentIdcode;
        public final String displayLabel;
        private final String sortKey;

        private ProjectionBucket(ProjectionBucketType type, String fragmentIdcode, String displayLabel, String sortKey) {
            this.type = type;
            this.fragmentIdcode = fragmentIdcode;
            this.displayLabel = displayLabel;
            this.sortKey = sortKey;
        }

        public static ProjectionBucket substituent(String fragmentIdcode, String fragmentSmiles) {
            return new ProjectionBucket(
                    ProjectionBucketType.SUBSTITUENT,
                    fragmentIdcode,
                    fragmentSmiles.isBlank() ? "[fragment]" : fragmentSmiles,
                    "1|" + fragmentSmiles + "|" + fragmentIdcode
            );
        }

        public static ProjectionBucket unsubstituted() {
            return new ProjectionBucket(ProjectionBucketType.UNSUBSTITUTED, null, "[unsubstituted]", "0|");
        }

        public static ProjectionBucket multiAttachment() {
            return new ProjectionBucket(ProjectionBucketType.MULTI_ATTACHMENT, null, "[multi-attachment]", "2|");
        }

        public static ProjectionBucket unmatched() {
            return new ProjectionBucket(ProjectionBucketType.UNMATCHED, null, "[unmatched]", "3|");
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof ProjectionBucket other)) {
                return false;
            }
            return type == other.type
                    && Objects.equals(fragmentIdcode, other.fragmentIdcode)
                    && Objects.equals(displayLabel, other.displayLabel);
        }

        @Override
        public int hashCode() {
            return Objects.hash(type, fragmentIdcode, displayLabel);
        }
    }

    public static final class OneDimProjection {
        public final int scaffoldAtom;
        public final String exitVectorLabel;
        public final List<OneDimProjectionRow> rows;

        private OneDimProjection(int scaffoldAtom, String exitVectorLabel, List<OneDimProjectionRow> rows) {
            this.scaffoldAtom = scaffoldAtom;
            this.exitVectorLabel = exitVectorLabel;
            this.rows = List.copyOf(rows);
        }
    }

    public static final class OneDimProjectionRow {
        public final ProjectionBucket bucket;
        public final List<Integer> compoundIndices;

        private OneDimProjectionRow(ProjectionBucket bucket, List<Integer> compoundIndices) {
            this.bucket = bucket;
            this.compoundIndices = List.copyOf(compoundIndices);
        }
    }

    public static final class TwoDimProjection {
        public final int rowScaffoldAtom;
        public final int columnScaffoldAtom;
        public final String rowExitVectorLabel;
        public final String columnExitVectorLabel;
        public final List<ProjectionBucket> rowBuckets;
        public final List<ProjectionBucket> columnBuckets;
        public final int[][] counts;
        private final Map<CellKey, List<Integer>> compoundIndicesByCell;

        private TwoDimProjection(
                int rowScaffoldAtom,
                int columnScaffoldAtom,
                String rowExitVectorLabel,
                String columnExitVectorLabel,
                List<ProjectionBucket> rowBuckets,
                List<ProjectionBucket> columnBuckets,
                int[][] counts,
                Map<CellKey, List<Integer>> compoundIndicesByCell
        ) {
            this.rowScaffoldAtom = rowScaffoldAtom;
            this.columnScaffoldAtom = columnScaffoldAtom;
            this.rowExitVectorLabel = rowExitVectorLabel;
            this.columnExitVectorLabel = columnExitVectorLabel;
            this.rowBuckets = List.copyOf(rowBuckets);
            this.columnBuckets = List.copyOf(columnBuckets);
            this.counts = counts;
            this.compoundIndicesByCell = Map.copyOf(compoundIndicesByCell);
        }

        public List<Integer> compoundIndices(int rowIndex, int columnIndex) {
            return compoundIndicesByCell.getOrDefault(new CellKey(rowIndex, columnIndex), List.of());
        }
    }

    private record CellKey(int rowIndex, int columnIndex) {}
}

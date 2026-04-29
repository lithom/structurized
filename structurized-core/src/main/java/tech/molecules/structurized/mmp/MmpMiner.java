package tech.molecules.structurized.mmp;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * In-memory matched molecular pair miner.
 */
public final class MmpMiner {
    private MmpMiner() {}

    public static MmpMiningResult mine(List<MmpInputCompound> compounds, MmpMiningConfig config) {
        Objects.requireNonNull(compounds, "compounds");
        Objects.requireNonNull(config, "config");

        Map<String, MmpInputCompound> inputById = new HashMap<>();
        List<MmpFragmentationRecord> records = new ArrayList<>();
        for (MmpInputCompound compound : compounds) {
            inputById.put(compound.compoundId(), compound);
            records.addAll(MmpFragmenter.fragment(compound, config));
        }
        records.sort(Comparator.comparing(MmpFragmentationRecord::canonicalRecordId));

        List<MmpPair> pairs = buildPairs(records, inputById, config);
        List<MmpTransformStats> stats = MmpStatsAggregator.aggregate(pairs, config);
        return new MmpMiningResult(records, pairs, stats);
    }

    private static List<MmpPair> buildPairs(
            List<MmpFragmentationRecord> records,
            Map<String, MmpInputCompound> inputById,
            MmpMiningConfig config
    ) {
        Map<KeyGroup, List<MmpFragmentationRecord>> byKey = new HashMap<>();
        for (MmpFragmentationRecord record : records) {
            byKey.computeIfAbsent(new KeyGroup(record.cutCount(), record.keyIdcode()), ignored -> new ArrayList<>()).add(record);
        }

        List<MmpPair> pairs = new ArrayList<>();
        List<KeyGroup> sortedGroups = byKey.keySet().stream().sorted().toList();
        for (KeyGroup group : sortedGroups) {
            List<MmpFragmentationRecord> groupRecords = byKey.get(group).stream()
                    .sorted(Comparator.comparing(MmpFragmentationRecord::compoundId)
                            .thenComparing(MmpFragmentationRecord::valueIdcode))
                    .toList();
            int emittedForKey = 0;
            for (int i = 0; i < groupRecords.size(); i++) {
                for (int j = i + 1; j < groupRecords.size(); j++) {
                    MmpFragmentationRecord a = groupRecords.get(i);
                    MmpFragmentationRecord b = groupRecords.get(j);
                    if (a.compoundId().equals(b.compoundId()) || a.valueIdcode().equals(b.valueIdcode())) {
                        continue;
                    }
                    pairs.add(toPair(a, b, inputById));
                    emittedForKey++;
                    if (config.emitReverseTransforms()) {
                        pairs.add(toPair(b, a, inputById));
                        emittedForKey++;
                    }
                    if (emittedForKey >= config.maxPairsPerKey()) {
                        break;
                    }
                }
                if (emittedForKey >= config.maxPairsPerKey()) {
                    break;
                }
            }
        }
        pairs.sort(Comparator.comparing(MmpPair::transformId)
                .thenComparing(MmpPair::compoundIdA)
                .thenComparing(MmpPair::compoundIdB)
                .thenComparing(MmpPair::keyIdcode));
        return pairs;
    }

    private static MmpPair toPair(
            MmpFragmentationRecord from,
            MmpFragmentationRecord to,
            Map<String, MmpInputCompound> inputById
    ) {
        Double valueA = inputById.get(from.compoundId()).value();
        Double valueB = inputById.get(to.compoundId()).value();
        return new MmpPair(
                from.compoundId(),
                to.compoundId(),
                valueA,
                valueB,
                null,
                from.keyIdcode(),
                from.valueIdcode(),
                to.valueIdcode(),
                null,
                from.cutCount()
        );
    }

    private record KeyGroup(int cutCount, String keyIdcode) implements Comparable<KeyGroup> {
        @Override
        public int compareTo(KeyGroup other) {
            int cutCmp = Integer.compare(cutCount, other.cutCount);
            if (cutCmp != 0) {
                return cutCmp;
            }
            return keyIdcode.compareTo(other.keyIdcode);
        }
    }
}

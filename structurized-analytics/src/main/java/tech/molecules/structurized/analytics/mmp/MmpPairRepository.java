package tech.molecules.structurized.analytics.mmp;

import tech.molecules.structurized.mmp.MmpFragmentationRecord;
import tech.molecules.structurized.mmp.MmpPair;

import java.util.List;

/**
 * Persistence API for structure-derived MMP fragmentations and pairs.
 */
public interface MmpPairRepository {
    void replaceFragmentationRecords(String universeId, List<MmpFragmentationRecord> records);

    void replacePairs(String universeId, List<MmpPair> pairs);

    List<MmpPair> listPairs(String universeId);
}

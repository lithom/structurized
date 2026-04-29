package tech.molecules.structurized.mmp;

import java.util.List;
import java.util.Objects;

/**
 * One canonical MMP fragmentation for one compound.
 */
public record MmpFragmentationRecord(
        String compoundId,
        int cutCount,
        String keyIdcode,
        String valueIdcode,
        int keyHeavyAtomCount,
        int valueHeavyAtomCount,
        List<Integer> cutBondIndices,
        String canonicalRecordId
) {
    public MmpFragmentationRecord {
        if (compoundId == null || compoundId.isBlank()) {
            throw new IllegalArgumentException("compoundId must not be blank");
        }
        if (cutCount < 1) {
            throw new IllegalArgumentException("cutCount must be positive");
        }
        Objects.requireNonNull(keyIdcode, "keyIdcode");
        Objects.requireNonNull(valueIdcode, "valueIdcode");
        cutBondIndices = List.copyOf(cutBondIndices == null ? List.of() : cutBondIndices);
        if (canonicalRecordId == null || canonicalRecordId.isBlank()) {
            canonicalRecordId = compoundId + "|" + cutCount + "|" + keyIdcode + "|" + valueIdcode;
        }
    }
}

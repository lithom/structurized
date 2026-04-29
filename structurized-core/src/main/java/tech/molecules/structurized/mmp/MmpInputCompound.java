package tech.molecules.structurized.mmp;

import com.actelion.research.chem.StereoMolecule;

import java.util.Objects;

/**
 * Input record for in-memory MMP mining.
 *
 * <p>The molecule is copied on construction because OpenChemLib molecules are mutable.</p>
 */
public record MmpInputCompound(String compoundId, StereoMolecule molecule, Double value) {
    public MmpInputCompound {
        if (compoundId == null || compoundId.isBlank()) {
            throw new IllegalArgumentException("compoundId must not be blank");
        }
        Objects.requireNonNull(molecule, "molecule");
        molecule = new StereoMolecule(molecule);
    }

    @Override
    public StereoMolecule molecule() {
        return new StereoMolecule(molecule);
    }
}

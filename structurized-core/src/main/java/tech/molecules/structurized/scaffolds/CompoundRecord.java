package tech.molecules.structurized.scaffolds;

import com.actelion.research.chem.StereoMolecule;

import java.util.Arrays;
import java.util.List;

/**
 * One input compound together with its precomputed neighborhood information.
 */
public final class CompoundRecord {
    public final int index;
    public final StereoMolecule molecule;
    public final String idcode;
    public final long[] fingerprint;
    public final List<Integer> neighborIndices;

    public CompoundRecord(
            int index,
            StereoMolecule molecule,
            String idcode,
            long[] fingerprint,
            List<Integer> neighborIndices
    ) {
        this.index = index;
        this.molecule = new StereoMolecule(molecule);
        this.idcode = idcode;
        this.fingerprint = Arrays.copyOf(fingerprint, fingerprint.length);
        this.neighborIndices = List.copyOf(neighborIndices);
    }
}

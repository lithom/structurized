package tech.molecules.structurized;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;

import java.util.BitSet;

/**
 * Minimal utility showcasing OpenChemLib usage.
 */
public final class OpenChemLibUtil {

    private OpenChemLibUtil() {}

    /**
     * Parses the given SMILES and returns the heavy-atom count.
     *
     * @param smiles SMILES string (e.g., "CCO", "c1ccccc1")
     * @return number of atoms in the parsed molecule (implicit H not counted)
     * @throws IllegalArgumentException if parsing fails
     */
    public static int atomCountFromSmiles(String smiles) {
        if (smiles == null || smiles.isEmpty()) {
            throw new IllegalArgumentException("SMILES must not be null or empty");
        }

        StereoMolecule mol = new StereoMolecule();
        try {
            SmilesParser parser = new SmilesParser();
            parser.parse(mol, smiles);
            // Ensure helper arrays are populated for downstream operations if needed
            mol.ensureHelperArrays(StereoMolecule.cHelperCIP);
            return mol.getAllAtoms();
        } catch (Exception e) {
            throw new IllegalArgumentException("Invalid SMILES: " + smiles, e);
        }
    }


    public static boolean[] bitsetToBool(BitSet bs, int length) {
        boolean b[] = new boolean[length];
        for(int i=0;i<length;i++) {b[i]=bs.get(i);}
        return b;
    }

    /**
     * Counts non-hydrogen atoms, irrespective of whether hydrogens are explicit or implicit.
     */
    public static int heavyAtomCount(StereoMolecule mol) {
        int count = 0;
        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            if (mol.getAtomicNo(atom) > 1) {
                count++;
            }
        }
        return count;
    }

}

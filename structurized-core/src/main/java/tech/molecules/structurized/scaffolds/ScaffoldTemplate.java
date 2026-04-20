package tech.molecules.structurized.scaffolds;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Canonical scaffold graph together with candidate exit vectors and atom symmetry classes.
 *
 * <p>This first scaffold-mode implementation intentionally keeps the template minimal:
 * it stores only the chemical graph of the unsubstituted scaffold. Exit vectors are not
 * predefined manually. Instead, every scaffold atom is a candidate exit vector and real
 * occupied positions are inferred from matched compounds.</p>
 */
public final class ScaffoldTemplate {
    public final StereoMolecule scaffold;
    public final String idcode;
    public final int[] atomSymmetryClasses;
    public final List<ExitVector> candidateExitVectors;

    private ScaffoldTemplate(
            StereoMolecule scaffold,
            String idcode,
            int[] atomSymmetryClasses,
            List<ExitVector> candidateExitVectors
    ) {
        this.scaffold = scaffold;
        this.idcode = idcode;
        this.atomSymmetryClasses = atomSymmetryClasses;
        this.candidateExitVectors = List.copyOf(candidateExitVectors);
    }

    public static ScaffoldTemplate create(StereoMolecule source) {
        StereoMolecule scaffold = new StereoMolecule(source);
        scaffold.setFragment(false);
        scaffold.ensureHelperArrays(Molecule.cHelperSymmetrySimple);

        int[] atomSymmetryClasses = new int[scaffold.getAtoms()];
        List<ExitVector> exitVectors = new ArrayList<>(scaffold.getAtoms());
        for (int atom = 0; atom < scaffold.getAtoms(); atom++) {
            atomSymmetryClasses[atom] = scaffold.getSymmetryRank(atom);
            exitVectors.add(new ExitVector(atom, atomSymmetryClasses[atom], scaffold.getAtomicNo(atom)));
        }

        return new ScaffoldTemplate(
                scaffold,
                new Canonizer(scaffold).getIDCode(),
                Arrays.copyOf(atomSymmetryClasses, atomSymmetryClasses.length),
                exitVectors
        );
    }

    /**
     * Creates a display-only scaffold copy with labeled exit-vector pseudo atoms attached to the
     * provided scaffold atoms.
     *
     * <p>The canonical scaffold stored in this template remains unchanged. The added pseudo atoms
     * exist only to make scaffold depictions easier to interpret in validation GUIs.</p>
     */
    public StereoMolecule createDisplayMoleculeWithExitVectors(List<Integer> exitVectorAtoms) {
        StereoMolecule display = new StereoMolecule(scaffold);
        display.setFragment(false);
        int scaffoldAtomCount = scaffold.getAtoms();

        Set<Integer> orderedAtoms = new LinkedHashSet<>(exitVectorAtoms);
        int labelIndex = 1;
        for (int atom : orderedAtoms) {
            if (atom < 0 || atom >= scaffoldAtomCount) {
                continue;
            }

            int exitVectorAtom = display.addAtom(0);
            display.setAtomCustomLabel(exitVectorAtom, "R" + labelIndex);
            display.addBond(atom, exitVectorAtom, Molecule.cBondTypeSingle);
            labelIndex++;
        }

        new CoordinateInventor(0).invent(display);
        display.ensureHelperArrays(Molecule.cHelperSymmetrySimple);
        return display;
    }
}

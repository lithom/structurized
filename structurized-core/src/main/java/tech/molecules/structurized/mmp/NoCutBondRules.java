package tech.molecules.structurized.mmp;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.List;

/**
 * Built-in chemistry veto rules for MMP bond cutting.
 */
public final class NoCutBondRules {
    private NoCutBondRules() {}

    /**
     * Conservative defaults intended to suppress common non-medicinal-chemistry cuts.
     */
    public static List<NoCutBondRule> defaultRules() {
        return List.of(new AmideBondRule());
    }

    /**
     * Creates a generic substructure-backed rule that vetoes selected bonds from every query match.
     */
    public static NoCutBondRule substructureRule(String id, StereoMolecule query, int... queryBondIndicesToVeto) {
        return new SubstructureNoCutBondRule(id, query, queryBondIndicesToVeto);
    }

    /**
     * Blocks amide C-N bonds, i.e. a single C-N bond where the carbon is a carbonyl carbon.
     */
    public static final class AmideBondRule implements NoCutBondRule {
        @Override
        public void markNoCutBonds(StereoMolecule mol, boolean[] noCutBond) {
            mol.ensureHelperArrays(Molecule.cHelperRings);
            for (int bond = 0; bond < mol.getBonds(); bond++) {
                if (mol.getBondOrder(bond) != 1) {
                    continue;
                }
                int a1 = mol.getBondAtom(0, bond);
                int a2 = mol.getBondAtom(1, bond);
                if (isCarbonylCarbon(mol, a1) && mol.getAtomicNo(a2) == 7) {
                    noCutBond[bond] = true;
                } else if (isCarbonylCarbon(mol, a2) && mol.getAtomicNo(a1) == 7) {
                    noCutBond[bond] = true;
                }
            }
        }

        private boolean isCarbonylCarbon(StereoMolecule mol, int atom) {
            if (mol.getAtomicNo(atom) != 6) {
                return false;
            }
            for (int i = 0; i < mol.getConnAtoms(atom); i++) {
                int bond = mol.getConnBond(atom, i);
                int neighbor = mol.getConnAtom(atom, i);
                if (mol.getAtomicNo(neighbor) == 8 && mol.getBondOrder(bond) == 2) {
                    return true;
                }
            }
            return false;
        }
    }
}

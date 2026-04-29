package tech.molecules.structurized.mmp;

import com.actelion.research.chem.StereoMolecule;

/**
 * Marks bonds that must not be considered as MMP cut bonds.
 */
public interface NoCutBondRule {
    /**
     * Set {@code noCutBond[bond] = true} for all bonds blocked by this rule.
     */
    void markNoCutBonds(StereoMolecule mol, boolean[] noCutBond);
}

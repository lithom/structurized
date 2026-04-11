package tech.molecules.structurized;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class OpenChemLibUtilTest {

    @Test
    void countsAtomsForSimpleMolecules() {
        assertEquals(3, OpenChemLibUtil.atomCountFromSmiles("CCO")); // ethanol skeleton
        assertEquals(6, OpenChemLibUtil.atomCountFromSmiles("c1ccccc1")); // benzene ring
    }

    @Test
    void rejectsInvalidSmiles() {
        assertThrows(IllegalArgumentException.class, () -> OpenChemLibUtil.atomCountFromSmiles("not_a_smiles"));
        assertThrows(IllegalArgumentException.class, () -> OpenChemLibUtil.atomCountFromSmiles(""));
    }
}


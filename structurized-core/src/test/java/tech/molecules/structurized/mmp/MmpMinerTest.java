package tech.molecules.structurized.mmp;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class MmpMinerTest {

    @Test
    void minesPairsForSharedArylCore() throws Exception {
        MmpMiningConfig config = MmpMiningConfig.builder()
                .maxCuts(1)
                .minKeyHeavyAtoms(6)
                .maxVariableHeavyAtoms(4)
                .maxVariableToMolHeavyAtomFraction(1.0)
                .minTransformSupport(1)
                .build();

        MmpMiningResult result = MmpMiner.mine(List.of(
                new MmpInputCompound("toluene", parse("Cc1ccccc1"), 1.0),
                new MmpInputCompound("ethylbenzene", parse("CCc1ccccc1"), 3.5)
        ), config);

        assertFalse(result.fragmentationRecords().isEmpty());
        assertFalse(result.pairs().isEmpty());
        assertTrue(result.pairs().stream().anyMatch(pair ->
                pair.compoundIdA().equals("toluene")
                        && pair.compoundIdB().equals("ethylbenzene")
                        && pair.delta() == 2.5));
        assertFalse(result.transformStats().isEmpty());
    }

    @Test
    void aggregateComputesBasicDeltaStatistics() {
        MmpPair p1 = pair("a", "b", 1.0, 3.0, "K1", "V1", "V2");
        MmpPair p2 = pair("c", "d", 2.0, 6.0, "K2", "V1", "V2");
        MmpMiningConfig config = MmpMiningConfig.builder().minTransformSupport(2).build();

        List<MmpTransformStats> stats = MmpStatsAggregator.aggregate(List.of(p1, p2), config);

        assertEquals(1, stats.size());
        MmpTransformStats stat = stats.getFirst();
        assertEquals(2, stat.supportCount());
        assertEquals(3.0, stat.meanDelta(), 1.0e-12);
        assertEquals(3.0, stat.medianDelta(), 1.0e-12);
        assertEquals(1.0, stat.standardDeviation(), 1.0e-12);
        assertEquals(2.0, stat.minDelta(), 1.0e-12);
        assertEquals(4.0, stat.maxDelta(), 1.0e-12);
        assertEquals(1.0, stat.positiveFraction(), 1.0e-12);
    }

    private static MmpPair pair(String a, String b, double valueA, double valueB, String key, String from, String to) {
        return new MmpPair(a, b, valueA, valueB, null, key, from, to, null, 1);
    }

    private static StereoMolecule parse(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(Molecule.cHelperRings);
        return molecule;
    }
}

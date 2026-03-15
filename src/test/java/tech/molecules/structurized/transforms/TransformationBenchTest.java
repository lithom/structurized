package tech.molecules.structurized.transforms;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class TransformationBenchTest {

    @Test
    void demoProviderDoesNotCreateBogusZeroChangeGroups() throws Exception {
        StereoMolecule toluene = parse("c1ccccc1C");
        StereoMolecule fluorobenzene = parse("c1ccccc1F");

        TransformationSplitter.MCSMap mcs =
                new TransformationBenchDemo.OCLMCSFastProvider().computeStrictMCS(toluene, fluorobenzene);

        List<TransformationGroup> groups =
                TransformationSplitter.splitIntoTransformations(
                        toluene,
                        fluorobenzene,
                        mcs,
                        1,
                        TransformationSplitter.FeatureMask.DEFAULT
                );

        assertEquals(1, groups.size());
        assertEquals(TransformationType.REPLACEMENT, groups.getFirst().type);
        assertEquals(List.of(5), groups.getFirst().attachmentsA);
    }

    @Test
    void symmetricPairsEnumeratesBothDirectionsButNoSelfPairs() throws Exception {
        TransformationBench.Config cfg = new TransformationBench.Config();
        cfg.symmetricPairs = true;
        cfg.verbose = false;

        TransformationBench.Result result = TransformationBench.run(
                List.of(parse("CC"), parse("CCC")),
                new TransformationBench.StubStrictMCSProvider(),
                cfg
        );

        assertEquals(2, result.summary.nPairsTried);
        assertEquals(2, result.pairs.size());
        assertTrue(result.pairs.stream().anyMatch(p -> p.i == 0 && p.j == 1));
        assertTrue(result.pairs.stream().anyMatch(p -> p.i == 1 && p.j == 0));
        assertFalse(result.pairs.stream().anyMatch(p -> p.i == p.j));
    }

    @Test
    void maxPairsCapsAttemptedPairsEvenWhenMcsProviderReturnsNull() throws Exception {
        TransformationBench.Config cfg = new TransformationBench.Config();
        cfg.maxPairs = 1;
        cfg.verbose = false;

        TransformationBench.Result result = TransformationBench.run(
                List.of(parse("CC"), parse("CCC"), parse("CCCC")),
                new TransformationBench.StubStrictMCSProvider(),
                cfg
        );

        assertEquals(1, result.summary.nPairsTried);
        assertEquals(1, result.pairs.size());
    }

    private static StereoMolecule parse(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(StereoMolecule.cHelperRings);
        return molecule;
    }
}

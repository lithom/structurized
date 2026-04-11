package tech.molecules.structurized.scaffolds;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import org.junit.jupiter.api.Test;
import tech.molecules.structurized.transforms.TransformationBenchDemo;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ScaffoldDiscoveryEngineTest {

    @Test
    void discoversBenzeneScaffoldAcrossAromaticSeries() throws Exception {
        List<StereoMolecule> molecules = List.of(
                parse("c1ccccc1C"),
                parse("c1ccccc1F"),
                parse("COc1ccccc1"),
                parse("CCc1ccccc1"),
                parse("Cc1ccccc1C"),
                parse("CCO")
        );

        ScaffoldDiscoveryConfig cfg = new ScaffoldDiscoveryConfig();
        cfg.neighborCount = 4;
        cfg.minNeighborSimilarity = 0.10f;
        cfg.minScaffoldHeavyAtoms = 6;
        cfg.minSupport = 3;

        ScaffoldDiscoveryResult result = ScaffoldDiscoveryEngine.discover(
                molecules,
                new TransformationBenchDemo.OCLMCSFastProvider(),
                cfg
        );

        String benzeneIdcode = ScaffoldTemplate.create(parse("c1ccccc1")).idcode;
        ScaffoldCandidate benzene = result.candidates.stream()
                .filter(candidate -> candidate.template.idcode.equals(benzeneIdcode))
                .findFirst()
                .orElseThrow();

        assertEquals(6, result.compounds.size());
        assertFalse(result.candidates.isEmpty());
        assertEquals(5, benzene.supportCount);
        assertEquals(List.of(0, 1, 2, 3, 4), benzene.supportCompoundIndices);
        assertTrue(benzene.averageExplainedFraction > 0.70);
        assertEquals(6, benzene.scaffoldHeavyAtomCount);
        assertTrue(benzene.observedExitVectorCount >= 1);
    }

    @Test
    void neighborhoodsRespectRequestedNeighborCount() throws Exception {
        List<StereoMolecule> molecules = List.of(
                parse("c1ccccc1C"),
                parse("c1ccccc1F"),
                parse("COc1ccccc1"),
                parse("CCc1ccccc1")
        );

        ScaffoldDiscoveryConfig cfg = new ScaffoldDiscoveryConfig();
        cfg.neighborCount = 2;
        cfg.minNeighborSimilarity = 0.0f;
        cfg.minScaffoldHeavyAtoms = 6;

        ScaffoldDiscoveryResult result = ScaffoldDiscoveryEngine.discover(
                molecules,
                new TransformationBenchDemo.OCLMCSFastProvider(),
                cfg
        );

        assertTrue(result.compounds.stream().allMatch(record -> record.neighborIndices.size() <= 2));
    }

    private static StereoMolecule parse(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(StereoMolecule.cHelperSymmetrySimple);
        return molecule;
    }
}

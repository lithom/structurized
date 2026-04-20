package tech.molecules.structurized.scaffolds;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ScaffoldDatasetDecompositionTest {

    @Test
    void fullDatasetDecompositionBuildsOneDimAndTwoDimViews() throws Exception {
        ScaffoldTemplate template = ScaffoldTemplate.create(parse("c1ccccc1"));
        List<CompoundRecord> compounds = List.of(
                compound(0, "c1ccccc1C"),
                compound(1, "c1ccccc1F"),
                compound(2, "Cc1ccccc1C"),
                compound(3, "CCO")
        );

        ScaffoldDatasetDecomposition dataset = ScaffoldDatasetDecomposition.analyze(
                compounds,
                template,
                new ScaffoldAnalyzer.Config()
        );

        assertEquals(3, dataset.matchedCompoundCount);
        assertEquals(1, dataset.unmatchedCompoundCount);
        assertTrue(dataset.observedExitVectorAtoms.size() >= 2);
        assertEquals(4, dataset.records.size());

        int firstExitVector = dataset.observedExitVectorAtoms.getFirst();
        ScaffoldDatasetDecomposition.OneDimProjection oneDim = dataset.oneDimProjection(firstExitVector);
        assertEquals(4, oneDim.rows.stream().mapToInt(row -> row.compoundIndices.size()).sum());
        assertTrue(oneDim.rows.stream().anyMatch(row -> row.bucket.type == ScaffoldDatasetDecomposition.ProjectionBucketType.SUBSTITUENT));
        assertTrue(oneDim.rows.stream().anyMatch(row -> row.bucket.type == ScaffoldDatasetDecomposition.ProjectionBucketType.UNMATCHED));

        int secondExitVector = dataset.observedExitVectorAtoms.get(1);
        ScaffoldDatasetDecomposition.TwoDimProjection twoDim = dataset.twoDimProjection(firstExitVector, secondExitVector);
        int total = 0;
        for (int row = 0; row < twoDim.counts.length; row++) {
            for (int col = 0; col < twoDim.counts[row].length; col++) {
                total += twoDim.counts[row][col];
            }
        }
        assertEquals(4, total);
        assertFalse(twoDim.rowBuckets.isEmpty());
        assertFalse(twoDim.columnBuckets.isEmpty());
    }

    private static CompoundRecord compound(int index, String smiles) throws Exception {
        StereoMolecule molecule = parse(smiles);
        return new CompoundRecord(index, molecule, new Canonizer(molecule).getIDCode(), new long[8], List.of());
    }

    private static StereoMolecule parse(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(StereoMolecule.cHelperSymmetrySimple);
        return molecule;
    }
}

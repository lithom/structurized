package tech.molecules.structurized.scaffolds;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import org.junit.jupiter.api.Test;

import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ScaffoldAnalyzerTest {

    @Test
    void benzeneTemplateCarriesOneSymmetryClass() throws Exception {
        ScaffoldTemplate template = ScaffoldTemplate.create(parse("c1ccccc1"));

        assertEquals(6, template.scaffold.getAtoms());
        assertEquals(6, template.candidateExitVectors.size());
        Set<Integer> symmetryClasses = template.candidateExitVectors.stream()
                .map(exitVector -> exitVector.symmetryClass)
                .collect(Collectors.toSet());
        assertEquals(1, symmetryClasses.size());
    }

    @Test
    void benzeneToTolueneProducesSingleAttachmentSubstitution() throws Exception {
        ScaffoldTemplate template = ScaffoldTemplate.create(parse("c1ccccc1"));
        ScaffoldDecomposition decomposition =
                ScaffoldAnalyzer.analyze(template, parse("c1ccccc1C"), new ScaffoldAnalyzer.Config());

        assertNull(decomposition.failure);
        assertNotNull(decomposition.match);
        assertEquals(1, decomposition.transformationGroups.size());
        assertEquals(1, decomposition.substitutionEvents.size());
        assertEquals(SubstitutionEventType.SINGLE_ATTACHMENT, decomposition.substitutionEvents.getFirst().eventType);
        assertEquals(1, decomposition.substitutionEvents.getFirst().exitVectors.size());
        assertEquals(1, decomposition.match.uniqueMatchCount);
    }

    @Test
    void benzeneToXyleneProducesTwoIndependentSubstitutionEvents() throws Exception {
        ScaffoldTemplate template = ScaffoldTemplate.create(parse("c1ccccc1"));
        ScaffoldDecomposition decomposition =
                ScaffoldAnalyzer.analyze(template, parse("Cc1ccccc1C"), new ScaffoldAnalyzer.Config());

        assertNull(decomposition.failure);
        assertEquals(2, decomposition.substitutionEvents.size());
        assertTrue(decomposition.substitutionEvents.stream()
                .allMatch(event -> event.eventType == SubstitutionEventType.SINGLE_ATTACHMENT));
    }

    @Test
    void benzeneToNaphthaleneProducesMultiAttachmentEvent() throws Exception {
        ScaffoldTemplate template = ScaffoldTemplate.create(parse("c1ccccc1"));
        ScaffoldDecomposition decomposition =
                ScaffoldAnalyzer.analyze(template, parse("c1ccc2ccccc2c1"), new ScaffoldAnalyzer.Config());

        assertNull(decomposition.failure);
        assertFalse(decomposition.substitutionEvents.isEmpty());
        assertTrue(decomposition.substitutionEvents.stream()
                .anyMatch(event -> event.eventType == SubstitutionEventType.MULTI_ATTACHMENT));
    }

    private static StereoMolecule parse(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(StereoMolecule.cHelperSymmetrySimple);
        return molecule;
    }
}

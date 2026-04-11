package tech.molecules.structurized.context;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import org.junit.jupiter.api.Test;
import tech.molecules.structurized.transforms.TransformationBenchDemo;
import tech.molecules.structurized.transforms.TransformationGroup;
import tech.molecules.structurized.transforms.TransformationSplitter;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ContextFeatureExtractorTest {

    @Test
    void extractsSingleAttachmentContextFeatures() throws Exception {
        TransformationGroup group = firstGroupFor("c1ccccc1C", "c1ccccc1F");

        ContextFeatures features = ContextFeatureExtractor.computeAll(group.signature.contextShellIdcode);

        assertEquals(3, features.shellSummary().atomCount());
        assertEquals(2, features.shellSummary().bondCount());
        assertEquals(0, features.shellSummary().heteroAtomCount());
        assertEquals(3, features.shellSummary().aromaticAtomCount());
        assertEquals(3, features.shellSummary().ringAtomCount());
        assertEquals(1, features.shellSummary().attachmentCount());

        assertEquals(1, features.attachmentAtomFeatures().size());
        ContextFeatures.AttachmentAtomFeatures attachment = features.attachmentAtomFeatures().getFirst();
        assertEquals(0, attachment.attachmentIndex());
        assertEquals(6, attachment.atomicNo());
        assertTrue(attachment.aromatic());
        assertTrue(attachment.ringAtom());
        assertEquals(2, attachment.degreeWithinShell());
        assertEquals(0, attachment.heteroNeighborCountWithinShell());
        assertEquals(6, attachment.smallestRingSize());
        assertTrue(features.attachmentPairFeatures().isEmpty());
        assertEquals(-1, features.minAttachmentDistance());
        assertEquals(-1, features.maxAttachmentDistance());
    }

    @Test
    void extractsPairwiseAttachmentFeaturesFromMultiAttachmentContext() throws Exception {
        TransformationGroup replacement = firstGroupFor("c1ccccc1C", "c1ccncc1");

        ContextFeatures features = ContextFeatureExtractor.computeAll(replacement.signature.contextShellIdcode);

        assertEquals(2, features.shellSummary().attachmentCount());
        assertEquals(4, features.shellSummary().atomCount());
        assertEquals(2, features.shellSummary().bondCount());
        assertEquals(2, features.attachmentAtomFeatures().size());
        assertEquals(1, features.attachmentPairFeatures().size());

        ContextFeatures.AttachmentPairFeatures pair = features.attachmentPairFeatures().getFirst();
        assertEquals(0, pair.attachmentIndex1());
        assertEquals(1, pair.attachmentIndex2());
        assertEquals(4, pair.topologicalDistance());
        assertFalse(pair.sameRing());
        assertFalse(pair.sameAromaticSystem());
        assertFalse(pair.sameFusedRingSystem());
        assertEquals(4, features.minAttachmentDistance());
        assertEquals(4, features.maxAttachmentDistance());
        assertEquals(4, features.attachmentDistanceMatrix()[0][1]);
        assertEquals(4, features.attachmentDistanceMatrix()[1][0]);
    }

    @Test
    void contextShellEncodingPreservesAttachmentLabels() throws Exception {
        TransformationGroup group = firstGroupFor("c1ccccc1C", "c1ccncc1");
        ContextFeatures features = ContextFeatureExtractor.computeAll(group.signature.contextShellIdcode);

        assertFalse(group.signature.contextShellIdcode.isEmpty());
        assertEquals(List.of(0, 1), features.attachmentAtomFeatures().stream().map(ContextFeatures.AttachmentAtomFeatures::attachmentIndex).toList());
    }

    @Test
    void expandedRawContextCapturesBroaderNeighborhood() throws Exception {
        TransformationGroup group = firstGroupFor("c1ccccc1C", "c1ccccc1F");

        ContextFeatures canonical = ContextFeatureExtractor.computeAll(group.signature.contextShellIdcode);
        ContextFeatures expanded = ContextFeatureExtractor.computeAll(group.signature.expandedRawContextIdcode);

        assertEquals(1, group.signature.radiusR);
        assertEquals(3, group.signature.expandedRawContextRadius);
        assertFalse(group.signature.expandedRawContextIdcode.isEmpty());
        assertEquals(3, canonical.shellSummary().atomCount());
        assertEquals(6, expanded.shellSummary().atomCount());
        assertEquals(6, expanded.shellSummary().aromaticAtomCount());
        assertEquals(6, expanded.shellSummary().ringAtomCount());
    }

    private static TransformationGroup firstGroupFor(String smilesA, String smilesB) throws Exception {
        StereoMolecule a = parse(smilesA);
        StereoMolecule b = parse(smilesB);
        TransformationSplitter.MCSMap mcs =
                new TransformationBenchDemo.OCLMCSFastProvider().computeStrictMCS(a, b);
        List<TransformationGroup> groups =
                TransformationSplitter.splitIntoTransformations(
                        a,
                        b,
                        mcs,
                        1,
                        TransformationSplitter.FeatureMask.DEFAULT
                );
        return groups.getFirst();
    }

    private static StereoMolecule parse(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(StereoMolecule.cHelperRings);
        return molecule;
    }
}

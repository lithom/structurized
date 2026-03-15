package tech.molecules.structurized.context;

import java.util.List;

public record ContextFeatures(
        String contextShellIdcode,
        ShellSummaryFeatures shellSummary,
        List<AttachmentAtomFeatures> attachmentAtomFeatures,
        List<AttachmentPairFeatures> attachmentPairFeatures,
        int[][] attachmentDistanceMatrix,
        int minAttachmentDistance,
        int maxAttachmentDistance
) {
    public ContextFeatures {
        attachmentAtomFeatures = List.copyOf(attachmentAtomFeatures);
        attachmentPairFeatures = List.copyOf(attachmentPairFeatures);
        attachmentDistanceMatrix = copyMatrix(attachmentDistanceMatrix);
    }

    private static int[][] copyMatrix(int[][] matrix) {
        int[][] copy = new int[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            copy[i] = matrix[i].clone();
        }
        return copy;
    }

    public record ShellSummaryFeatures(
            int atomCount,
            int bondCount,
            int heteroAtomCount,
            int aromaticAtomCount,
            int ringAtomCount,
            int formalChargeSum,
            int formalChargeAbsSum,
            int attachmentCount
    ) {}

    public record AttachmentAtomFeatures(
            int attachmentIndex,
            int atomicNo,
            boolean aromatic,
            boolean ringAtom,
            int formalCharge,
            int atomPi,
            int smallestRingSize,
            int degreeWithinShell,
            int heteroNeighborCountWithinShell
    ) {}

    public record AttachmentPairFeatures(
            int attachmentIndex1,
            int attachmentIndex2,
            int topologicalDistance,
            boolean sameRing,
            boolean sameAromaticSystem,
            boolean sameFusedRingSystem
    ) {}
}

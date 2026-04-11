package tech.molecules.structurized.scaffolds;

import java.util.List;

/**
 * One discovered scaffold candidate with local support statistics.
 */
public final class ScaffoldCandidate {
    public final ScaffoldTemplate template;
    public final List<Integer> discoverySeedIndices;
    public final List<Integer> discoveryCompoundIndices;
    public final List<Integer> localSearchCompoundIndices;
    public final List<Integer> supportCompoundIndices;
    public final int supportCount;
    public final double averageExplainedFraction;
    public final int scaffoldHeavyAtomCount;
    public final int observedExitVectorCount;
    public final List<Integer> observedExitVectorAtoms;
    public final List<Integer> observedExitVectorSymmetryClasses;
    public final double combinedScore;

    public ScaffoldCandidate(
            ScaffoldTemplate template,
            List<Integer> discoverySeedIndices,
            List<Integer> discoveryCompoundIndices,
            List<Integer> localSearchCompoundIndices,
            List<Integer> supportCompoundIndices,
            int supportCount,
            double averageExplainedFraction,
            int scaffoldHeavyAtomCount,
            int observedExitVectorCount,
            List<Integer> observedExitVectorAtoms,
            List<Integer> observedExitVectorSymmetryClasses,
            double combinedScore
    ) {
        this.template = template;
        this.discoverySeedIndices = List.copyOf(discoverySeedIndices);
        this.discoveryCompoundIndices = List.copyOf(discoveryCompoundIndices);
        this.localSearchCompoundIndices = List.copyOf(localSearchCompoundIndices);
        this.supportCompoundIndices = List.copyOf(supportCompoundIndices);
        this.supportCount = supportCount;
        this.averageExplainedFraction = averageExplainedFraction;
        this.scaffoldHeavyAtomCount = scaffoldHeavyAtomCount;
        this.observedExitVectorCount = observedExitVectorCount;
        this.observedExitVectorAtoms = List.copyOf(observedExitVectorAtoms);
        this.observedExitVectorSymmetryClasses = List.copyOf(observedExitVectorSymmetryClasses);
        this.combinedScore = combinedScore;
    }
}

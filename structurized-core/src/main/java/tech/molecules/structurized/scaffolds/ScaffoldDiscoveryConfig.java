package tech.molecules.structurized.scaffolds;

import tech.molecules.structurized.transforms.TransformationSplitter;

/**
 * Configuration for the V1 scaffold discovery engine.
 */
public final class ScaffoldDiscoveryConfig {
    public int neighborCount = 4;
    public float minNeighborSimilarity = 0.15f;
    public int maxSeeds = Integer.MAX_VALUE;
    public boolean shuffleSeeds = false;
    public long randomSeed = 0L;
    public int minScaffoldHeavyAtoms = 5;
    public boolean requireConnectedScaffold = true;
    public int minSupport = 2;
    public int radiusR = 1;
    public int featureMask = TransformationSplitter.FeatureMask.DEFAULT;

    double combinedScore(int supportCount, double averageExplainedFraction, int scaffoldHeavyAtomCount, int observedExitVectorCount) {
        return supportCount * 1000.0
                + averageExplainedFraction * 100.0
                + scaffoldHeavyAtomCount * 10.0
                + observedExitVectorCount;
    }
}

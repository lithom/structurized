package tech.molecules.structurized.transforms;

import java.util.List;

/**
 * Full pairwise relation between two indexed compounds, including all independent local edit groups.
 */
public final class PairTransformation {
    public final int i;
    public final int j;
    public final int nGroups;
    public final List<TransformationGroup> groups;
    public final String failure;

    public PairTransformation(int i, int j, List<TransformationGroup> groups) {
        this.i = i;
        this.j = j;
        this.groups = List.copyOf(groups);
        this.nGroups = groups.size();
        this.failure = null;
    }

    public PairTransformation(int i, int j, String failure) {
        this.i = i;
        this.j = j;
        this.groups = List.of();
        this.nGroups = 0;
        this.failure = failure;
    }
}

package tech.molecules.structurized.scaffolds;

import tech.molecules.structurized.transforms.TransformationGroup;

import java.util.List;

/**
 * Scaffold-relative decomposition of one compound.
 */
public final class ScaffoldDecomposition {
    public final ScaffoldTemplate template;
    public final ScaffoldMatch match;
    public final List<TransformationGroup> transformationGroups;
    public final List<SubstitutionEvent> substitutionEvents;
    public final String failure;

    public ScaffoldDecomposition(
            ScaffoldTemplate template,
            ScaffoldMatch match,
            List<TransformationGroup> transformationGroups,
            List<SubstitutionEvent> substitutionEvents
    ) {
        this.template = template;
        this.match = match;
        this.transformationGroups = List.copyOf(transformationGroups);
        this.substitutionEvents = List.copyOf(substitutionEvents);
        this.failure = null;
    }

    public ScaffoldDecomposition(ScaffoldTemplate template, String failure) {
        this.template = template;
        this.match = null;
        this.transformationGroups = List.of();
        this.substitutionEvents = List.of();
        this.failure = failure;
    }
}

package tech.molecules.structurized.scaffolds;

import tech.molecules.structurized.transforms.TransformationGroup;
import tech.molecules.structurized.transforms.TransformationType;

import java.util.List;

/**
 * One scaffold-relative substitution event inferred from a scaffold/compound decomposition.
 */
public final class SubstitutionEvent {
    public final TransformationGroup transformationGroup;
    public final List<ExitVector> exitVectors;
    public final List<Integer> scaffoldAtoms;
    public final List<Integer> symmetryClasses;
    public final String addedFragmentIdcode;
    public final TransformationType transformationType;
    public final SubstitutionEventType eventType;

    public SubstitutionEvent(
            TransformationGroup transformationGroup,
            List<ExitVector> exitVectors,
            List<Integer> scaffoldAtoms,
            List<Integer> symmetryClasses,
            String addedFragmentIdcode,
            TransformationType transformationType,
            SubstitutionEventType eventType
    ) {
        this.transformationGroup = transformationGroup;
        this.exitVectors = List.copyOf(exitVectors);
        this.scaffoldAtoms = List.copyOf(scaffoldAtoms);
        this.symmetryClasses = List.copyOf(symmetryClasses);
        this.addedFragmentIdcode = addedFragmentIdcode;
        this.transformationType = transformationType;
        this.eventType = eventType;
    }
}

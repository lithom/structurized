package tech.molecules.structurized.scaffolds;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Full scaffold-relative decomposition result for one dataset compound.
 */
public final class CompoundDecompositionRecord {
    public final CompoundRecord compound;
    public final ScaffoldDecomposition decomposition;
    public final boolean matched;
    public final double explainedFraction;
    public final List<Integer> occupiedExitVectorAtoms;
    public final List<Integer> occupiedExitVectorSymmetryClasses;
    public final Map<Integer, ExitVectorAssignment> assignmentsByScaffoldAtom;
    public final Set<Integer> multiAttachmentAtoms;
    public final Set<Integer> ambiguousAtoms;
    public final boolean hasMultiAttachment;

    public CompoundDecompositionRecord(
            CompoundRecord compound,
            ScaffoldDecomposition decomposition,
            boolean matched,
            double explainedFraction,
            List<Integer> occupiedExitVectorAtoms,
            List<Integer> occupiedExitVectorSymmetryClasses,
            Map<Integer, ExitVectorAssignment> assignmentsByScaffoldAtom,
            Set<Integer> multiAttachmentAtoms,
            Set<Integer> ambiguousAtoms,
            boolean hasMultiAttachment
    ) {
        this.compound = compound;
        this.decomposition = decomposition;
        this.matched = matched;
        this.explainedFraction = explainedFraction;
        this.occupiedExitVectorAtoms = List.copyOf(occupiedExitVectorAtoms);
        this.occupiedExitVectorSymmetryClasses = List.copyOf(occupiedExitVectorSymmetryClasses);
        this.assignmentsByScaffoldAtom = Map.copyOf(new TreeMap<>(assignmentsByScaffoldAtom));
        this.multiAttachmentAtoms = Set.copyOf(new TreeSet<>(multiAttachmentAtoms));
        this.ambiguousAtoms = Set.copyOf(new TreeSet<>(ambiguousAtoms));
        this.hasMultiAttachment = hasMultiAttachment;
    }
}

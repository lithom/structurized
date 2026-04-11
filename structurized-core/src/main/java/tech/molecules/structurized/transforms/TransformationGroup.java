package tech.molecules.structurized.transforms;

import java.util.Collections;
import java.util.List;

/**
 * One independent local edit group within a pairwise A->B relation.
 */
public final class TransformationGroup {
    public final TransformationSignature signature;
    public final List<Integer> attachmentsA;
    public final TransformationType type;

    public TransformationGroup(TransformationSignature signature, List<Integer> attachmentsA, TransformationType type) {
        this.signature = signature;
        this.attachmentsA = Collections.unmodifiableList(attachmentsA);
        this.type = type;
    }
}

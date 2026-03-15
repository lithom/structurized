package tech.molecules.structurized.scaffolds;

/**
 * Candidate scaffold exit vector represented by one scaffold atom and its symmetry class.
 */
public final class ExitVector {
    public final int atomIndex;
    public final int symmetryClass;
    public final int atomicNo;

    public ExitVector(int atomIndex, int symmetryClass, int atomicNo) {
        this.atomIndex = atomIndex;
        this.symmetryClass = symmetryClass;
        this.atomicNo = atomicNo;
    }
}

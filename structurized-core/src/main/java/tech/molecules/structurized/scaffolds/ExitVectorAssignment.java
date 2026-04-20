package tech.molecules.structurized.scaffolds;

/**
 * Single-attachment substituent assignment for one scaffold exit vector in one compound.
 */
public final class ExitVectorAssignment {
    public final int scaffoldAtom;
    public final int symmetryClass;
    public final String fragmentIdcode;
    public final String fragmentSmiles;

    public ExitVectorAssignment(int scaffoldAtom, int symmetryClass, String fragmentIdcode, String fragmentSmiles) {
        this.scaffoldAtom = scaffoldAtom;
        this.symmetryClass = symmetryClass;
        this.fragmentIdcode = fragmentIdcode;
        this.fragmentSmiles = fragmentSmiles;
    }
}

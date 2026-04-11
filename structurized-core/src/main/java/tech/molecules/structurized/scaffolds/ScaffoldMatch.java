package tech.molecules.structurized.scaffolds;

import java.util.Arrays;

/**
 * Selected mapping of a scaffold template onto one compound.
 */
public final class ScaffoldMatch {
    public final ScaffoldTemplate template;
    public final int[] scaffoldToCompoundAtom;
    public final int[] compoundToScaffoldAtom;
    public final int uniqueMatchCount;
    public final int selectedMatchIndex;

    public ScaffoldMatch(
            ScaffoldTemplate template,
            int[] scaffoldToCompoundAtom,
            int[] compoundToScaffoldAtom,
            int uniqueMatchCount,
            int selectedMatchIndex
    ) {
        this.template = template;
        this.scaffoldToCompoundAtom = Arrays.copyOf(scaffoldToCompoundAtom, scaffoldToCompoundAtom.length);
        this.compoundToScaffoldAtom = Arrays.copyOf(compoundToScaffoldAtom, compoundToScaffoldAtom.length);
        this.uniqueMatchCount = uniqueMatchCount;
        this.selectedMatchIndex = selectedMatchIndex;
    }
}

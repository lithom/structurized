package tech.molecules.structurized.mmp;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;

import java.util.Arrays;
import java.util.Objects;

/**
 * Query-backed no-cut rule that vetoes selected query bonds in every substructure match.
 */
public final class SubstructureNoCutBondRule implements NoCutBondRule {
    private final String id;
    private final StereoMolecule query;
    private final int[] queryBondIndicesToVeto;

    public SubstructureNoCutBondRule(String id, StereoMolecule query, int... queryBondIndicesToVeto) {
        if (id == null || id.isBlank()) {
            throw new IllegalArgumentException("id must not be blank");
        }
        this.id = id;
        Objects.requireNonNull(query, "query");
        this.query = new StereoMolecule(query);
        this.query.setFragment(true);
        this.query.ensureHelperArrays(Molecule.cHelperRings);
        this.queryBondIndicesToVeto = Arrays.copyOf(queryBondIndicesToVeto, queryBondIndicesToVeto.length);
        for (int queryBond : this.queryBondIndicesToVeto) {
            if (queryBond < 0 || queryBond >= this.query.getBonds()) {
                throw new IllegalArgumentException("query bond index out of range: " + queryBond);
            }
        }
    }

    public String id() {
        return id;
    }

    @Override
    public void markNoCutBonds(StereoMolecule mol, boolean[] noCutBond) {
        SSSearcher searcher = new SSSearcher();
        searcher.setMol(query, mol);
        if (searcher.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode) == 0) {
            return;
        }
        for (int[] match : searcher.getMatchList()) {
            for (int queryBond : queryBondIndicesToVeto) {
                int qA = query.getBondAtom(0, queryBond);
                int qB = query.getBondAtom(1, queryBond);
                int tA = match[qA];
                int tB = match[qB];
                if (tA < 0 || tB < 0) {
                    continue;
                }
                int targetBond = mol.getBond(tA, tB);
                if (targetBond >= 0) {
                    noCutBond[targetBond] = true;
                }
            }
        }
    }
}

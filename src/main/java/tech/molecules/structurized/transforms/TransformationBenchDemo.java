package tech.molecules.structurized.transforms;

import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.SmilesParser;

import java.util.*;

/**
 * Toy demo that wires TransformationBench with an MCS provider based on OpenChemLib's MCSFast
 * and runs a full pairwise analysis on a small set of SMILES.
 *
 * NOTE: The exact MCSFast API varies slightly across OCL versions. The provider below
 * shows a STRICT usage pattern and validates bond identity inside the preserved core.
 * If a method name differs in your build, adjust the few marked lines accordingly
 * (search for "// <--- CHECK YOUR OCL").
 */
public class TransformationBenchDemo {

    public static void main(String[] args) {
        // --- 1) Small toy set (edit freely) ---
        List<String> smiles = List.of(
                "c1ccccc1C",        // toluene
                "c1ccccc1F",        // fluorobenzene
                "c1ccccc1Cl",       // chlorobenzene
                "COc1ccccc1",       // anisole
                "CCc1ccccc1",       // ethylbenzene
                "c1ccncc1"          // pyridine
        );

        List<StereoMolecule> mols = parseSmiles(smiles);

        // --- 2) Configure bench ---
        TransformationBench.Config cfg = new TransformationBench.Config();
        cfg.radiusR = 1;
        cfg.featureMask = TransformationSplitter.FeatureMask.DEFAULT;
        cfg.keepMultiCenter = true;
        cfg.symmetricPairs = false; // unordered pairs only
        cfg.verbose = true;
        cfg.maxPairs = Integer.MAX_VALUE;

        // --- 3) Run ---
        TransformationBench.Result res = TransformationBench.run(
                mols,
                new OCLMCSFastProvider(),
                cfg
        );

        // --- 4) Print summary & a few details ---
        res.printSummary();
        int shown = 0;
        for (PairTransformation pair : res.pairs) {
            if (pair.failure != null) continue;
            if (pair.groups.isEmpty()) continue;
            System.out.printf("Pair (%d,%d): %d TG(s)\n", pair.i, pair.j, pair.groups.size());
            for (TransformationGroup g : pair.groups) {
                System.out.printf("  - %s  sigId=%s  removedId=%s  addedId=%s  attachments=%s\n",
                        g.type,
                        g.signature.sigId.substring(0, 12),
                        shorten(g.signature.removedIdcode),
                        shorten(g.signature.addedIdcode),
                        g.attachmentsA);
            }
            if (++shown >= 8) break; // don't spam
        }
    }

    private static String shorten(String s) {
        if (s == null) return "-";
        return s.length() <= 18 ? s : s.substring(0, 18) + "…";
    }

    private static List<StereoMolecule> parseSmiles(List<String> smiles) {
        var sp = new SmilesParser();
        List<StereoMolecule> out = new ArrayList<>();
        for (String smi : smiles) {
            StereoMolecule m = new StereoMolecule();
            try { sp.parse(m, smi); } catch (Exception e) { throw new RuntimeException("SMILES parse failed: "+smi, e); }
            m.ensureHelperArrays(StereoMolecule.cHelperRings);
            out.add(m);
        }
        return out;
    }

    /** Strict MCS provider based on OpenChemLib's MCSFast. */
    public static class OCLMCSFastProvider implements TransformationBench.StrictMCSProvider {
        @Override
        public TransformationSplitter.MCSMap computeStrictMCS(StereoMolecule A, StereoMolecule B) {
            try {
                // ---- Construct MCSFast (API may vary by version) ----
                // <--- CHECK YOUR OCL: if constructor differs, adjust here
                com.actelion.research.chem.mcs.MCSFast mcs = new com.actelion.research.chem.mcs.MCSFast();
                mcs.set(A,B);
                StereoMolecule mcsMol = mcs.getMCS();
                if (mcsMol == null || mcsMol.getAtoms() == 0) {
                    return null;
                }


                // Enforce strict element & bond-order matching if your version supports flags
                // <--- CHECK YOUR OCL: uncomment & set strict modes if available
                // mcs.setAtomTypeMode(MCSFast.ATOM_ELEMENT);
                // mcs.setBondTypeMode(MCSFast.BOND_ORDER);

                // Run and obtain an atom mapping from A->B (unmapped = -1)
                // <--- CHECK YOUR OCL: adjust method name to your version
                int[] mapAtoB = firstMappingAtoB(mcsMol,A,B);
                if (mapAtoB == null) return null;

                // Validate STRICTNESS: all mapped bonds must be identical in order & aromaticity
                if (!atomsStrictlyMatch(A, B, mapAtoB)) return null;
                if (!bondsStrictlyMatch(A, B, mapAtoB)) return null;

                return new TransformationSplitter.MCSMap(mapAtoB, invertMap(mapAtoB,B.getAtoms()));
            } catch (Throwable t) {
                // If anything goes wrong, we signal "no strict MCS" for this pair
                return null;
            }
        }

        // --- helpers ---
        private static int[] firstMappingAtoB(StereoMolecule mcs, StereoMolecule A, StereoMolecule B) {
            SSSearcher ss = new SSSearcher();
            mcs.setFragment(true);
            ss.setMol(mcs,A);
            if (ss.findFragmentInMolecule(SSSearcher.cCountModeFirstMatch, SSSearcher.cDefaultMatchMode) == 0) {
                return null;
            }
            int[] mcsToA = ss.getMatchList().getFirst();

            ss.setMol(mcs,B);
            if (ss.findFragmentInMolecule(SSSearcher.cCountModeFirstMatch, SSSearcher.cDefaultMatchMode) == 0) {
                return null;
            }
            int[] mcsToB = ss.getMatchList().getFirst();

            int[] AtoMCS = invertMap(mcsToA,A.getAtoms());

            int[] AtoB = new int[AtoMCS.length];
            for(int a=0;a<AtoB.length;a++) {
                if(AtoMCS[a]>=0) {
                    AtoB[a] = mcsToB[ AtoMCS[a] ];
                }
                else {
                    AtoB[a] = -1;
                }
            }
            return AtoB;
        }

        private static boolean atomsStrictlyMatch(StereoMolecule A, StereoMolecule B, int[] a2b) {
            for (int a = 0; a < A.getAllAtoms(); a++) {
                int b = a2b[a];
                if (b < 0) {
                    continue;
                }
                if (A.getAtomicNo(a) != B.getAtomicNo(b)) return false;
                if (A.getAtomCharge(a) != B.getAtomCharge(b)) return false;
                if (A.getAtomMass(a) != B.getAtomMass(b)) return false;
                if (A.isAromaticAtom(a) != B.isAromaticAtom(b)) return false;
            }
            return true;
        }

        private static int[] invertMap(int[] mapAtoB, int lengthB) {
            //int maxB = -1;
            //for (int v : mapAtoB) maxB = Math.max(maxB, v);
            //int[] inv = new int[maxB + 1];
            int[] inv = new int[lengthB];
            Arrays.fill(inv, -1);
            for (int a = 0; a < mapAtoB.length; a++) {
                int b = mapAtoB[a]; if (b >= 0) inv[b] = a;
            }
            return inv;
        }

        private static boolean bondsStrictlyMatch(StereoMolecule A, StereoMolecule B, int[] a2b) {
            for (int a = 0; a < A.getAllAtoms(); a++) {
                int b = a2b[a];
                if (b < 0) continue;
                for (int i = 0; i < A.getConnAtoms(a); i++) {
                    int aN = A.getConnAtom(a, i);
                    int bN = a2b[aN];
                    if (bN < 0) continue; // neighbor not mapped → outside core
                    int bondA = A.getConnBond(a, i);
                    int bondB = B.getBond(b, bN);
                    if (bondB == -1) return false; // must exist in core
                    // strict order & type equality
                    if (A.getBondOrder(bondA) != B.getBondOrder(bondB)) return false;
                    if (A.getBondType(bondA) != B.getBondType(bondB)) return false;
                }
            }
            return true;
        }
    }
}


package tech.molecules.structurized.context;

import java.util.ArrayList;
import java.util.List;

/**
 * Encodes parent-aware context-shell annotations into OpenChemLib custom atom labels.
 *
 * <p>The context shell itself is stored as an IDCode of a copied shell fragment. Since that copied shell
 * may truncate larger parent-core ring systems, the custom label carries the parent-derived semantics that
 * must survive canonization and later decoding.</p>
 *
 * <p>Every encoded atom label may contain:
 * - attachment index, if the atom is one of the ordered attachments
 * - parent-core ring membership
 * - parent-core aromaticity
 * - parent-core atom pi value
 * - parent-core smallest ring size
 * - pairwise relation data from this attachment to later attachment indices</p>
 */
public final class ContextLabelCodec {
    private static final String PREFIX = "CTX";

    private ContextLabelCodec() {}

    /**
     * Convert structured per-atom annotation data into a compact custom-label payload.
     */
    public static String encode(ContextAtomData atomData) {
        StringBuilder sb = new StringBuilder(PREFIX);
        if (atomData.attachmentIndex() != null) {
            sb.append(";ATT=").append(atomData.attachmentIndex());
        }
        sb.append(";PR=").append(atomData.parentRingAtom() ? 1 : 0);
        sb.append(";PA=").append(atomData.parentAromaticAtom() ? 1 : 0);
        sb.append(";PP=").append(atomData.parentAtomPi());
        sb.append(";PS=").append(atomData.parentSmallestRingSize());
        if (!atomData.pairRelations().isEmpty()) {
            sb.append(";PAIR=");
            for (int i = 0; i < atomData.pairRelations().size(); i++) {
                PairRelation relation = atomData.pairRelations().get(i);
                if (i > 0) {
                    sb.append("|");
                }
                sb.append(relation.otherAttachmentIndex()).append(":")
                        .append(relation.sameRing() ? 1 : 0).append(":")
                        .append(relation.sameAromaticSystem() ? 1 : 0).append(":")
                        .append(relation.sameFusedRingSystem() ? 1 : 0).append(":")
                        .append(relation.topologicalDistance());
            }
        }
        return sb.toString();
    }

    /**
     * Decode a custom-label payload back into structured context annotation data.
     *
     * @return decoded annotation, or {@code null} if the label is absent or not a context label
     */
    public static ContextAtomData decode(String label) {
        if (label == null || !label.startsWith(PREFIX)) {
            return null;
        }

        Integer attachmentIndex = null;
        boolean parentRingAtom = false;
        boolean parentAromaticAtom = false;
        int parentAtomPi = 0;
        int parentSmallestRingSize = 0;
        List<PairRelation> pairRelations = new ArrayList<>();

        String[] parts = label.split(";");
        for (int i = 1; i < parts.length; i++) {
            String part = parts[i];
            if (part.startsWith("ATT=")) {
                attachmentIndex = Integer.parseInt(part.substring(4));
            } else if (part.startsWith("PR=")) {
                parentRingAtom = "1".equals(part.substring(3));
            } else if (part.startsWith("PA=")) {
                parentAromaticAtom = "1".equals(part.substring(3));
            } else if (part.startsWith("PP=")) {
                parentAtomPi = Integer.parseInt(part.substring(3));
            } else if (part.startsWith("PS=")) {
                parentSmallestRingSize = Integer.parseInt(part.substring(3));
            } else if (part.startsWith("PAIR=")) {
                String pairData = part.substring(5);
                if (!pairData.isEmpty()) {
                    String[] pairEntries = pairData.split("\\|");
                    for (String entry : pairEntries) {
                        String[] fields = entry.split(":");
                        pairRelations.add(new PairRelation(
                                Integer.parseInt(fields[0]),
                                "1".equals(fields[1]),
                                "1".equals(fields[2]),
                                "1".equals(fields[3]),
                                Integer.parseInt(fields[4])
                        ));
                    }
                }
            }
        }

        return new ContextAtomData(
                attachmentIndex,
                parentRingAtom,
                parentAromaticAtom,
                parentAtomPi,
                parentSmallestRingSize,
                pairRelations
        );
    }

    /**
     * Parent-aware annotation for one atom of a copied context shell.
     */
    public record ContextAtomData(
            Integer attachmentIndex,
            boolean parentRingAtom,
            boolean parentAromaticAtom,
            int parentAtomPi,
            int parentSmallestRingSize,
            List<PairRelation> pairRelations
    ) {
        public ContextAtomData {
            pairRelations = List.copyOf(pairRelations);
        }
    }

    /**
     * Pairwise relation from one attachment atom to another ordered attachment atom in the same shell.
     */
    public record PairRelation(
            int otherAttachmentIndex,
            boolean sameRing,
            boolean sameAromaticSystem,
            boolean sameFusedRingSystem,
            int topologicalDistance
    ) {}
}

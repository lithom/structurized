package tech.molecules.structurized.transforms;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.nio.charset.StandardCharsets;

/**
 * Canonical identity of one local transformation group at a chosen context granularity.
 */
public final class TransformationSignature {
    public final String removedIdcode;
    public final String addedIdcode;
    public final String attachmentPattern;
    /**
     * Canonical raw context shell around the ordered attachment atoms.
     *
     * This is an OpenChemLib IDCode of the copied core-only shell fragment, with custom atom labels
     * enabled. The labels encode attachment identity plus parent-core annotations such as ring and
     * aromatic membership, atom pi value, smallest ring size, and attachment-pair relation data.
     */
    public final String contextShellIdcode;
    /**
     * Larger raw context shell intended for downstream analysis and feature derivation.
     *
     * This uses the same attachment-aware and parent-aware encoding as {@link #contextShellIdcode},
     * but is built with a larger fixed neighborhood radius. It is not included in {@link #sigId}.
     */
    public final String expandedRawContextIdcode;
    public final int radiusR;
    public final int expandedRawContextRadius;
    public final int featureMask;
    public final String rxnClass;
    public final String sigId;

    public TransformationSignature(
            String removedIdcode,
            String addedIdcode,
            String attachmentPattern,
            String contextShellIdcode,
            String expandedRawContextIdcode,
            int radiusR,
            int expandedRawContextRadius,
            int featureMask,
            String rxnClass
    ) {
        this.removedIdcode = removedIdcode;
        this.addedIdcode = addedIdcode;
        this.attachmentPattern = attachmentPattern;
        this.contextShellIdcode = contextShellIdcode;
        this.expandedRawContextIdcode = expandedRawContextIdcode;
        this.radiusR = radiusR;
        this.expandedRawContextRadius = expandedRawContextRadius;
        this.featureMask = featureMask;
        this.rxnClass = rxnClass;
        this.sigId = hashFields(
                removedIdcode,
                addedIdcode,
                attachmentPattern,
                contextShellIdcode,
                Integer.toString(radiusR),
                Integer.toString(featureMask),
                rxnClass
        );
    }

    private static String hashFields(String... fields) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            for (String field : fields) {
                md.update((field == null ? "" : field).getBytes(StandardCharsets.UTF_8));
                md.update((byte) '|');
            }
            byte[] hash = md.digest();
            StringBuilder sb = new StringBuilder();
            for (byte b : hash) {
                sb.append(String.format("%02x", b));
            }
            return sb.toString();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
    }
}

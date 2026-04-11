package tech.molecules.structurized.gui;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Minimal reader for simple SMILES text files.
 *
 * <p>Accepted format is intentionally simple:
 * every non-empty non-comment line contributes the first token as the SMILES string.
 * This supports common .smi files with optional name columns.</p>
 */
public final class SmilesInputReader {
    private SmilesInputReader() {}

    public static List<String> readSmilesFile(Path path) throws IOException {
        return parseSmilesLines(Files.readAllLines(path));
    }

    static List<String> parseSmilesLines(List<String> lines) {
        List<String> smiles = new ArrayList<>();
        boolean firstNonEmpty = true;
        for (String rawLine : lines) {
            if (rawLine == null) {
                continue;
            }
            String line = rawLine.trim();
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }

            String firstToken = firstToken(line);
            if (firstToken.isEmpty()) {
                continue;
            }

            if (firstNonEmpty && firstToken.equalsIgnoreCase("smiles")) {
                firstNonEmpty = false;
                continue;
            }

            firstNonEmpty = false;
            smiles.add(firstToken);
        }
        return smiles;
    }

    private static String firstToken(String line) {
        int limit = line.length();
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (Character.isWhitespace(c) || c == ',' || c == ';') {
                limit = i;
                break;
            }
        }
        return line.substring(0, limit).trim();
    }
}

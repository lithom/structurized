package tech.molecules.structurized.gui;

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class SmilesInputReaderTest {

    @Test
    void readsSimpleSmilesLinesWithOptionalNamesAndHeader() {
        List<String> smiles = SmilesInputReader.parseSmilesLines(List.of(
                "smiles name",
                "# comment",
                "",
                "c1ccccc1 benzene",
                "CCO,ethanol",
                "CCN"
        ));

        assertEquals(List.of("c1ccccc1", "CCO", "CCN"), smiles);
    }
}

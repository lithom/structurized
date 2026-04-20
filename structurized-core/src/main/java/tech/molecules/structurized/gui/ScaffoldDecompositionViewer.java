package tech.molecules.structurized.gui;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.IsomericSmilesCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JStructureView;
import tech.molecules.structurized.scaffolds.CompoundDecompositionRecord;
import tech.molecules.structurized.scaffolds.ExitVectorAssignment;
import tech.molecules.structurized.scaffolds.ScaffoldDatasetDecomposition;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JToggleButton;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableRowSorter;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Tabbed internal viewer for full-dataset scaffold-relative decomposition.
 */
public final class ScaffoldDecompositionViewer {
    private final ScaffoldDatasetDecomposition dataset;
    private final JFrame frame;
    private final JStructureView scaffoldView;
    private final StereoMolecule scaffoldDisplayMolecule;
    private final Map<Integer, Integer> exitVectorMarkerAtoms;
    private final JTextArea summaryArea;
    private final CompoundTableModel compoundTableModel;
    private final JTable compoundTable;
    private final JStructureView compoundView;
    private final JTextArea compoundDetailArea;
    private final List<ExitVectorChoice> exitVectorChoices;
    private final List<JToggleButton> oneDimButtons;
    private final ButtonGroup oneDimButtonGroup;
    private final JPanel oneDimButtonPanel;
    private final OneDimProjectionTableModel oneDimTableModel;
    private final JTable oneDimTable;
    private final JStructureView oneDimFragmentView;
    private final JTextArea oneDimDetailArea;
    private final List<JToggleButton> twoDimButtons;
    private final JPanel twoDimButtonPanel;
    private final JLabel twoDimSelectionLabel;
    private final TwoDimProjectionTableModel twoDimTableModel;
    private final JTable twoDimTable;
    private final JTextArea twoDimDetailArea;
    private ExitVectorChoice selectedOneDimChoice;
    private ExitVectorChoice hoveredOneDimChoice;
    private final List<ExitVectorChoice> selectedTwoDimChoices;
    private ExitVectorChoice hoveredTwoDimChoice;
    private JTabbedPane tabs;

    private ScaffoldDecompositionViewer(ScaffoldDatasetDecomposition dataset) {
        this.dataset = dataset;
        this.frame = new JFrame("Scaffold Decomposition Viewer");
        this.scaffoldView = new JStructureView();
        this.scaffoldDisplayMolecule = dataset.template.createDisplayMoleculeWithExitVectors(dataset.observedExitVectorAtoms);
        this.exitVectorMarkerAtoms = createExitVectorMarkerMap();
        this.summaryArea = new JTextArea();
        this.compoundTableModel = new CompoundTableModel(dataset.records.stream().filter(record -> record.matched).toList());
        this.compoundTable = new JTable(compoundTableModel);
        this.compoundView = new JStructureView();
        this.compoundDetailArea = new JTextArea();
        this.exitVectorChoices = exitVectorChoices(dataset);
        this.oneDimButtons = new ArrayList<>();
        this.oneDimButtonGroup = new ButtonGroup();
        this.oneDimButtonPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        this.oneDimTableModel = new OneDimProjectionTableModel();
        this.oneDimTable = new JTable(oneDimTableModel);
        this.oneDimFragmentView = new JStructureView();
        this.oneDimDetailArea = new JTextArea();
        this.twoDimButtons = new ArrayList<>();
        this.twoDimButtonPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        this.twoDimSelectionLabel = new JLabel("Select two exit vectors.");
        this.twoDimTableModel = new TwoDimProjectionTableModel();
        this.twoDimTable = new JTable(twoDimTableModel);
        this.twoDimDetailArea = new JTextArea();
        this.selectedTwoDimChoices = new ArrayList<>();
    }

    public static void show(JFrame owner, ScaffoldDatasetDecomposition dataset) {
        SwingUtilities.invokeLater(() -> new ScaffoldDecompositionViewer(dataset).show(owner));
    }

    private void show(JFrame owner) {
        frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        frame.setLayout(new BorderLayout(8, 8));
        frame.add(buildHeaderPanel(), BorderLayout.NORTH);
        frame.add(buildTabbedPane(), BorderLayout.CENTER);
        frame.setSize(new Dimension(1500, 900));
        frame.setLocationRelativeTo(owner);
        frame.setVisible(true);

        if (!compoundTableModel.records.isEmpty()) {
            compoundTable.setRowSelectionInterval(0, 0);
        }
        if (!exitVectorChoices.isEmpty()) {
            selectOneDimChoice(exitVectorChoices.getFirst());
            updateOneDimProjection();
        }
        if (!exitVectorChoices.isEmpty()) {
            selectTwoDimChoice(exitVectorChoices.getFirst());
            if (exitVectorChoices.size() > 1) {
                selectTwoDimChoice(exitVectorChoices.get(1));
            }
            updateTwoDimProjection();
        }
        applyHighlightForCurrentTab();
    }

    private JPanel buildHeaderPanel() {
        scaffoldView.setBorder(BorderFactory.createTitledBorder("Selected Scaffold"));
        scaffoldView.setPreferredSize(new Dimension(420, 220));
        scaffoldView.structureChanged(scaffoldDisplayMolecule);

        summaryArea.setEditable(false);
        summaryArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        summaryArea.setLineWrap(true);
        summaryArea.setWrapStyleWord(true);
        summaryArea.setText("""
                Scaffold SMILES: %s
                Scaffold IDCode: %s
                Matched Compounds: %d
                Unmatched Compounds: %d
                Compounds With Multi-Attachment Events: %d
                Observed Exit Vectors: %s
                """.formatted(
                new IsomericSmilesCreator(dataset.template.scaffold).getSmiles(),
                dataset.template.idcode,
                dataset.matchedCompoundCount,
                dataset.unmatchedCompoundCount,
                dataset.multiAttachmentCompoundCount,
                dataset.observedExitVectorAtoms.stream().map(dataset::exitVectorLabel).toList()
        ));
        summaryArea.setCaretPosition(0);

        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBorder(BorderFactory.createEmptyBorder(8, 8, 0, 8));
        panel.add(scaffoldView, BorderLayout.WEST);
        panel.add(new JScrollPane(summaryArea), BorderLayout.CENTER);
        return panel;
    }

    private JTabbedPane buildTabbedPane() {
        tabs = new JTabbedPane();
        tabs.addTab("Compounds", buildCompoundsTab());
        tabs.addTab("1D Projection", buildOneDimTab());
        tabs.addTab("2D Projection", buildTwoDimTab());
        tabs.addChangeListener(event -> applyHighlightForCurrentTab());
        return tabs;
    }

    private JSplitPane buildCompoundsTab() {
        compoundTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        compoundTable.setAutoCreateRowSorter(true);
        compoundTable.setRowHeight(24);
        compoundTable.setRowSorter(new TableRowSorter<>(compoundTableModel));
        setNumericRenderers(compoundTable, new int[]{0, 2, 3});
        compoundTable.getColumnModel().getColumn(1).setCellRenderer(new NumberRenderer(new DecimalFormat("0.000")));
        compoundTable.getSelectionModel().addListSelectionListener(event -> {
            if (!event.getValueIsAdjusting()) {
                updateCompoundDetail();
            }
        });

        compoundView.setBorder(BorderFactory.createTitledBorder("Compound"));
        compoundView.setPreferredSize(new Dimension(420, 280));
        configureTextArea(compoundDetailArea);

        JPanel detailPanel = new JPanel(new BorderLayout(8, 8));
        detailPanel.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 8));
        detailPanel.add(compoundView, BorderLayout.NORTH);
        detailPanel.add(new JScrollPane(compoundDetailArea), BorderLayout.CENTER);

        JSplitPane split = new JSplitPane(
                JSplitPane.HORIZONTAL_SPLIT,
                new JScrollPane(compoundTable),
                detailPanel
        );
        split.setResizeWeight(0.62);
        return split;
    }

    private JPanel buildOneDimTab() {
        JPanel top = new JPanel(new BorderLayout(8, 8));
        top.add(new JLabel("Exit Vector"), BorderLayout.WEST);
        buildOneDimButtons();
        top.add(oneDimButtonPanel, BorderLayout.CENTER);

        oneDimTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        oneDimTable.setAutoCreateRowSorter(true);
        oneDimTable.setRowHeight(24);
        oneDimTable.setRowSorter(new TableRowSorter<>(oneDimTableModel));
        setNumericRenderers(oneDimTable, new int[]{1});
        oneDimTable.getSelectionModel().addListSelectionListener(event -> {
            if (!event.getValueIsAdjusting()) {
                updateOneDimDetail();
            }
        });

        oneDimFragmentView.setBorder(BorderFactory.createTitledBorder("Substituent"));
        oneDimFragmentView.setPreferredSize(new Dimension(360, 240));
        configureTextArea(oneDimDetailArea);

        JPanel detailPanel = new JPanel(new BorderLayout(8, 8));
        detailPanel.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 8));
        detailPanel.add(oneDimFragmentView, BorderLayout.NORTH);
        detailPanel.add(new JScrollPane(oneDimDetailArea), BorderLayout.CENTER);

        JSplitPane split = new JSplitPane(
                JSplitPane.HORIZONTAL_SPLIT,
                new JScrollPane(oneDimTable),
                detailPanel
        );
        split.setResizeWeight(0.62);

        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));
        panel.add(top, BorderLayout.NORTH);
        panel.add(split, BorderLayout.CENTER);
        return panel;
    }

    private JPanel buildTwoDimTab() {
        JPanel selectors = new JPanel(new BorderLayout(8, 8));
        buildTwoDimButtons();
        selectors.add(twoDimButtonPanel, BorderLayout.CENTER);
        selectors.add(twoDimSelectionLabel, BorderLayout.SOUTH);

        twoDimTable.setCellSelectionEnabled(true);
        twoDimTable.setAutoCreateRowSorter(false);
        twoDimTable.setRowSelectionAllowed(true);
        twoDimTable.setColumnSelectionAllowed(true);
        twoDimTable.setRowHeight(24);
        twoDimTable.getSelectionModel().addListSelectionListener(event -> {
            if (!event.getValueIsAdjusting()) {
                updateTwoDimDetail();
            }
        });
        twoDimTable.getColumnModel().getSelectionModel().addListSelectionListener(event -> {
            if (!event.getValueIsAdjusting()) {
                updateTwoDimDetail();
            }
        });

        configureTextArea(twoDimDetailArea);

        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));
        panel.add(selectors, BorderLayout.NORTH);
        panel.add(new JScrollPane(twoDimTable), BorderLayout.CENTER);
        panel.add(new JScrollPane(twoDimDetailArea), BorderLayout.SOUTH);
        return panel;
    }

    private void updateCompoundDetail() {
        int viewRow = compoundTable.getSelectedRow();
        if (viewRow < 0) {
            compoundView.structureChanged();
            compoundDetailArea.setText("");
            return;
        }

        int modelRow = compoundTable.convertRowIndexToModel(viewRow);
        CompoundDecompositionRecord record = compoundTableModel.getRecordAt(modelRow);
        compoundView.structureChanged(record.compound.molecule);
        compoundDetailArea.setText("""
                Compound Index: %d
                Compound IDCode: %s
                Matched: %s
                Failure: %s
                Explained Fraction: %.3f
                Occupied Exit Vectors: %s
                Multi-Attachment: %s
                Assignments: %s
                Event Count: %d
                """.formatted(
                record.compound.index,
                record.compound.idcode,
                record.matched,
                record.decomposition.failure,
                record.explainedFraction,
                record.occupiedExitVectorAtoms.stream().map(dataset::exitVectorLabel).toList(),
                record.hasMultiAttachment,
                assignmentSummary(record),
                record.decomposition.substitutionEvents.size()
        ));
        compoundDetailArea.setCaretPosition(0);
    }

    private void updateOneDimProjection() {
        ExitVectorChoice choice = selectedOneDimChoice;
        if (choice == null) {
            oneDimTableModel.setProjection(null);
            oneDimFragmentView.structureChanged();
            oneDimDetailArea.setText("");
            applyHighlightForCurrentTab();
            return;
        }

        oneDimTableModel.setProjection(dataset.oneDimProjection(choice.atomIndex, false));
        setNumericRenderers(oneDimTable, new int[]{1});
        if (oneDimTableModel.getRowCount() > 0) {
            oneDimTable.setRowSelectionInterval(0, 0);
        } else {
            oneDimFragmentView.structureChanged();
            oneDimDetailArea.setText("");
        }
        applyHighlightForCurrentTab();
    }

    private void updateOneDimDetail() {
        int viewRow = oneDimTable.getSelectedRow();
        if (viewRow < 0 || oneDimTableModel.projection == null) {
            oneDimFragmentView.structureChanged();
            oneDimDetailArea.setText("");
            return;
        }

        int modelRow = oneDimTable.convertRowIndexToModel(viewRow);
        ScaffoldDatasetDecomposition.OneDimProjectionRow row = oneDimTableModel.projection.rows.get(modelRow);
        setFragmentView(oneDimFragmentView, row.bucket.fragmentIdcode);
        oneDimDetailArea.setText("""
                Exit Vector: %s
                Bucket: %s
                Compound Count: %d
                Compound Indices: %s
                """.formatted(
                oneDimTableModel.projection.exitVectorLabel,
                row.bucket.displayLabel,
                row.compoundIndices.size(),
                row.compoundIndices
        ));
        oneDimDetailArea.setCaretPosition(0);
    }

    private void updateTwoDimProjection() {
        if (selectedTwoDimChoices.size() < 2) {
            twoDimTableModel.setProjection(null);
            twoDimDetailArea.setText("Select two exit vectors to build the 2D projection.");
            twoDimSelectionLabel.setText(buildTwoDimSelectionText());
            applyHighlightForCurrentTab();
            return;
        }

        ExitVectorChoice rowChoice = selectedTwoDimChoices.get(0);
        ExitVectorChoice columnChoice = selectedTwoDimChoices.get(1);
        twoDimTableModel.setProjection(dataset.twoDimProjection(rowChoice.atomIndex, columnChoice.atomIndex));
        for (int column = 1; column < twoDimTable.getColumnModel().getColumnCount(); column++) {
            twoDimTable.getColumnModel().getColumn(column).setCellRenderer(new NumberRenderer(new DecimalFormat("0")));
        }
        twoDimSelectionLabel.setText(buildTwoDimSelectionText());
        if (twoDimTableModel.getRowCount() > 0 && twoDimTableModel.getColumnCount() > 1) {
            twoDimTable.setRowSelectionInterval(0, 0);
            twoDimTable.setColumnSelectionInterval(1, 1);
        } else {
            twoDimDetailArea.setText("");
        }
        applyHighlightForCurrentTab();
    }

    private void updateTwoDimDetail() {
        int row = twoDimTable.getSelectedRow();
        int column = twoDimTable.getSelectedColumn();
        ScaffoldDatasetDecomposition.TwoDimProjection projection = twoDimTableModel.projection;
        if (projection == null || row < 0 || column <= 0) {
            twoDimDetailArea.setText("");
            return;
        }

        int modelRow = twoDimTable.convertRowIndexToModel(row);
        int modelColumn = twoDimTable.convertColumnIndexToModel(column);
        int bucketColumn = modelColumn - 1;
        List<Integer> compoundIndices = projection.compoundIndices(modelRow, bucketColumn);
        twoDimDetailArea.setText("""
                Row Exit Vector: %s
                Column Exit Vector: %s
                Row Bucket: %s
                Column Bucket: %s
                Compound Count: %d
                Compound Indices: %s
                """.formatted(
                projection.rowExitVectorLabel,
                projection.columnExitVectorLabel,
                projection.rowBuckets.get(modelRow).displayLabel,
                projection.columnBuckets.get(bucketColumn).displayLabel,
                compoundIndices.size(),
                compoundIndices
        ));
        twoDimDetailArea.setCaretPosition(0);
    }

    private void buildOneDimButtons() {
        oneDimButtonPanel.removeAll();
        oneDimButtons.clear();
        for (ExitVectorChoice choice : exitVectorChoices) {
            JToggleButton button = new JToggleButton(choice.shortLabel());
            button.setToolTipText(choice.label);
            button.addActionListener(event -> {
                if (button.isSelected()) {
                    selectOneDimChoice(choice);
                    updateOneDimProjection();
                }
            });
            button.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseEntered(MouseEvent e) {
                    hoveredOneDimChoice = choice;
                    applyHighlightForCurrentTab();
                }

                @Override
                public void mouseExited(MouseEvent e) {
                    hoveredOneDimChoice = null;
                    applyHighlightForCurrentTab();
                }
            });
            oneDimButtonGroup.add(button);
            oneDimButtons.add(button);
            oneDimButtonPanel.add(button);
        }
    }

    private void buildTwoDimButtons() {
        twoDimButtonPanel.removeAll();
        twoDimButtons.clear();
        for (ExitVectorChoice choice : exitVectorChoices) {
            JToggleButton button = new JToggleButton(choice.shortLabel());
            button.setToolTipText(choice.label);
            button.addActionListener(event -> {
                toggleTwoDimChoice(choice, button.isSelected());
                updateTwoDimProjection();
            });
            button.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseEntered(MouseEvent e) {
                    hoveredTwoDimChoice = choice;
                    applyHighlightForCurrentTab();
                }

                @Override
                public void mouseExited(MouseEvent e) {
                    hoveredTwoDimChoice = null;
                    applyHighlightForCurrentTab();
                }
            });
            twoDimButtons.add(button);
            twoDimButtonPanel.add(button);
        }
        twoDimSelectionLabel.setText(buildTwoDimSelectionText());
    }

    private void selectOneDimChoice(ExitVectorChoice choice) {
        selectedOneDimChoice = choice;
        for (int i = 0; i < exitVectorChoices.size(); i++) {
            oneDimButtons.get(i).setSelected(exitVectorChoices.get(i).equals(choice));
        }
    }

    private void selectTwoDimChoice(ExitVectorChoice choice) {
        toggleTwoDimChoice(choice, true);
    }

    private void toggleTwoDimChoice(ExitVectorChoice choice, boolean selected) {
        if (selected) {
            selectedTwoDimChoices.remove(choice);
            selectedTwoDimChoices.add(choice);
            if (selectedTwoDimChoices.size() > 2) {
                ExitVectorChoice removed = selectedTwoDimChoices.remove(0);
                int removedIndex = exitVectorChoices.indexOf(removed);
                if (removedIndex >= 0) {
                    twoDimButtons.get(removedIndex).setSelected(false);
                }
            }
        } else {
            selectedTwoDimChoices.remove(choice);
        }

        for (int i = 0; i < exitVectorChoices.size(); i++) {
            twoDimButtons.get(i).setSelected(selectedTwoDimChoices.contains(exitVectorChoices.get(i)));
        }
        twoDimSelectionLabel.setText(buildTwoDimSelectionText());
    }

    private String buildTwoDimSelectionText() {
        String rowText = selectedTwoDimChoices.size() >= 1 ? selectedTwoDimChoices.get(0).shortLabel() : "-";
        String columnText = selectedTwoDimChoices.size() >= 2 ? selectedTwoDimChoices.get(1).shortLabel() : "-";
        return "Row = " + rowText + "    Column = " + columnText;
    }

    private void applyHighlightForCurrentTab() {
        if (tabs == null) {
            return;
        }
        int selectedTab = tabs.getSelectedIndex();
        if (selectedTab == 1) {
            applyScaffoldHighlight(
                    selectedOneDimChoice == null ? null : selectedOneDimChoice.atomIndex,
                    null,
                    hoveredOneDimChoice == null ? null : hoveredOneDimChoice.atomIndex
            );
        } else if (selectedTab == 2) {
            Integer rowAtom = selectedTwoDimChoices.size() >= 1 ? selectedTwoDimChoices.get(0).atomIndex : null;
            Integer columnAtom = selectedTwoDimChoices.size() >= 2 ? selectedTwoDimChoices.get(1).atomIndex : null;
            Integer hoverAtom = hoveredTwoDimChoice == null ? null : hoveredTwoDimChoice.atomIndex;
            applyScaffoldHighlight(rowAtom, columnAtom, hoverAtom);
        } else {
            scaffoldView.setAtomHighlightColors(null, null);
            scaffoldView.repaint();
        }
    }

    private void applyScaffoldHighlight(Integer rowAtom, Integer columnAtom, Integer hoverAtom) {
        int atomCount = scaffoldDisplayMolecule.getAllAtoms();
        int[] colors = new int[atomCount];
        float[] radii = new float[atomCount];

        highlightAtom(colors, radii, rowAtom, 0x88ffb000);
        highlightAtom(colors, radii, columnAtom, 0x8800b7ff);
        if (hoverAtom != null && !hoverAtom.equals(rowAtom) && !hoverAtom.equals(columnAtom)) {
            highlightAtom(colors, radii, hoverAtom, 0x8896ff00);
        }

        if (isAnyHighlight(colors)) {
            scaffoldView.setAtomHighlightColors(colors, radii);
        } else {
            scaffoldView.setAtomHighlightColors(null, null);
        }
        scaffoldView.repaint();
    }

    private void highlightAtom(int[] colors, float[] radii, Integer scaffoldAtom, int argb) {
        if (scaffoldAtom == null) {
            return;
        }
        if (scaffoldAtom >= 0 && scaffoldAtom < colors.length) {
            colors[scaffoldAtom] = argb;
            radii[scaffoldAtom] = 0.42f;
        }
        Integer markerAtom = exitVectorMarkerAtoms.get(scaffoldAtom);
        if (markerAtom != null && markerAtom >= 0 && markerAtom < colors.length) {
            colors[markerAtom] = argb;
            radii[markerAtom] = 0.52f;
        }
    }

    private static boolean isAnyHighlight(int[] colors) {
        for (int color : colors) {
            if (color != 0) {
                return true;
            }
        }
        return false;
    }

    private Map<Integer, Integer> createExitVectorMarkerMap() {
        Map<Integer, Integer> map = new HashMap<>();
        for (int atom = 0; atom < scaffoldDisplayMolecule.getAllAtoms(); atom++) {
            String label = scaffoldDisplayMolecule.getAtomCustomLabel(atom);
            if (label == null || !label.startsWith("R")) {
                continue;
            }
            try {
                int index = Integer.parseInt(label.substring(1)) - 1;
                if (index >= 0 && index < dataset.observedExitVectorAtoms.size()) {
                    map.put(dataset.observedExitVectorAtoms.get(index), atom);
                }
            } catch (NumberFormatException ignored) {
                // ignore malformed display-only labels
            }
        }
        return map;
    }

    private static List<ExitVectorChoice> exitVectorChoices(ScaffoldDatasetDecomposition dataset) {
        List<ExitVectorChoice> choices = new ArrayList<>();
        for (int atom : dataset.observedExitVectorAtoms) {
            choices.add(new ExitVectorChoice(atom, dataset.exitVectorLabel(atom), "R" + (choices.size() + 1)));
        }
        return choices;
    }

    private static void setNumericRenderers(JTable table, int[] columns) {
        DefaultTableCellRenderer renderer = new NumberRenderer(new DecimalFormat("0"));
        for (int column : columns) {
            table.getColumnModel().getColumn(column).setCellRenderer(renderer);
        }
    }

    private static void configureTextArea(JTextArea textArea) {
        textArea.setEditable(false);
        textArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        textArea.setLineWrap(true);
        textArea.setWrapStyleWord(true);
    }

    private static void setFragmentView(JStructureView view, String fragmentIdcode) {
        if (fragmentIdcode == null || fragmentIdcode.isBlank()) {
            view.structureChanged();
            return;
        }
        StereoMolecule mol = new StereoMolecule();
        new IDCodeParser().parse(mol, fragmentIdcode);
        view.structureChanged(mol);
    }

    private static String assignmentSummary(CompoundDecompositionRecord record) {
        if (!record.matched) {
            return "[unmatched]";
        }
        List<String> parts = new ArrayList<>();
        for (ExitVectorAssignmentView assignment : record.assignmentsByScaffoldAtom.entrySet().stream()
                .map(entry -> new ExitVectorAssignmentView(entry.getKey(), entry.getValue()))
                .sorted()
                .toList()) {
            parts.add("Atom " + (assignment.atomIndex + 1) + "=" + assignment.assignment.fragmentSmiles);
        }
        if (record.hasMultiAttachment) {
            parts.add("[multi-attachment]");
        }
        if (parts.isEmpty()) {
            return "[unsubstituted]";
        }
        return String.join(", ", parts);
    }

    private record ExitVectorAssignmentView(int atomIndex, ExitVectorAssignment assignment) implements Comparable<ExitVectorAssignmentView> {
        @Override
        public int compareTo(ExitVectorAssignmentView other) {
            return Integer.compare(atomIndex, other.atomIndex);
        }
    }

    private record ExitVectorChoice(int atomIndex, String label, String shortLabel) {
        public String shortLabel() {
            return shortLabel;
        }

        @Override
        public String toString() {
            return label;
        }
    }

    private static final class NumberRenderer extends DefaultTableCellRenderer {
        private final DecimalFormat format;

        private NumberRenderer(DecimalFormat format) {
            this.format = format;
            setHorizontalAlignment(JLabel.RIGHT);
        }

        @Override
        protected void setValue(Object value) {
            if (value instanceof Number number) {
                super.setValue(format.format(number));
                return;
            }
            super.setValue(value);
        }
    }

    private static final class CompoundTableModel extends AbstractTableModel {
        private static final String[] COLUMNS = {
                "Compound", "Explained", "Events", "Occupied EVs", "Matched", "Multi", "Assignments"
        };
        private final List<CompoundDecompositionRecord> records;

        private CompoundTableModel(List<CompoundDecompositionRecord> records) {
            this.records = List.copyOf(records);
        }

        public CompoundDecompositionRecord getRecordAt(int rowIndex) {
            return records.get(rowIndex);
        }

        @Override
        public int getRowCount() {
            return records.size();
        }

        @Override
        public int getColumnCount() {
            return COLUMNS.length;
        }

        @Override
        public String getColumnName(int column) {
            return COLUMNS[column];
        }

        @Override
        public Class<?> getColumnClass(int columnIndex) {
            return switch (columnIndex) {
                case 0, 2, 3 -> Integer.class;
                case 1 -> Double.class;
                case 4, 5, 6 -> String.class;
                default -> Object.class;
            };
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            CompoundDecompositionRecord record = records.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> record.compound.index;
                case 1 -> record.explainedFraction;
                case 2 -> record.decomposition.substitutionEvents.size();
                case 3 -> record.occupiedExitVectorAtoms.size();
                case 4 -> record.matched ? "yes" : "no";
                case 5 -> record.hasMultiAttachment ? "yes" : "no";
                case 6 -> assignmentSummary(record);
                default -> "";
            };
        }
    }

    private static final class OneDimProjectionTableModel extends AbstractTableModel {
        private static final String[] COLUMNS = {"Bucket", "Count", "Compounds"};
        private ScaffoldDatasetDecomposition.OneDimProjection projection;

        private void setProjection(ScaffoldDatasetDecomposition.OneDimProjection projection) {
            this.projection = projection;
            fireTableStructureChanged();
        }

        @Override
        public int getRowCount() {
            return projection == null ? 0 : projection.rows.size();
        }

        @Override
        public int getColumnCount() {
            return COLUMNS.length;
        }

        @Override
        public String getColumnName(int column) {
            return COLUMNS[column];
        }

        @Override
        public Class<?> getColumnClass(int columnIndex) {
            return switch (columnIndex) {
                case 1 -> Integer.class;
                default -> String.class;
            };
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            ScaffoldDatasetDecomposition.OneDimProjectionRow row = projection.rows.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> row.bucket.displayLabel;
                case 1 -> row.compoundIndices.size();
                case 2 -> row.compoundIndices.stream().limit(8).map(Object::toString).collect(Collectors.joining(", "));
                default -> "";
            };
        }
    }

    private static final class TwoDimProjectionTableModel extends AbstractTableModel {
        private ScaffoldDatasetDecomposition.TwoDimProjection projection;

        private void setProjection(ScaffoldDatasetDecomposition.TwoDimProjection projection) {
            this.projection = projection;
            fireTableStructureChanged();
        }

        @Override
        public int getRowCount() {
            return projection == null ? 0 : projection.rowBuckets.size();
        }

        @Override
        public int getColumnCount() {
            return projection == null ? 0 : projection.columnBuckets.size() + 1;
        }

        @Override
        public String getColumnName(int column) {
            if (projection == null) {
                return "";
            }
            if (column == 0) {
                return projection.rowExitVectorLabel;
            }
            return projection.columnBuckets.get(column - 1).displayLabel;
        }

        @Override
        public Class<?> getColumnClass(int columnIndex) {
            return columnIndex == 0 ? String.class : Integer.class;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            if (columnIndex == 0) {
                return projection.rowBuckets.get(rowIndex).displayLabel;
            }
            return projection.counts[rowIndex][columnIndex - 1];
        }
    }
}

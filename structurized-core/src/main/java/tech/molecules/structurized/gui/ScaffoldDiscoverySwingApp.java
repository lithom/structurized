package tech.molecules.structurized.gui;

import com.actelion.research.chem.IsomericSmilesCreator;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JStructureView;
import tech.molecules.structurized.scaffolds.ScaffoldCandidate;
import tech.molecules.structurized.scaffolds.ScaffoldDiscoveryConfig;
import tech.molecules.structurized.scaffolds.ScaffoldDiscoveryEngine;
import tech.molecules.structurized.scaffolds.ScaffoldDiscoveryResult;
import tech.molecules.structurized.transforms.TransformationBenchDemo;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.UIManager;
import javax.swing.WindowConstants;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableRowSorter;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Minimal Swing GUI for loading a SMILES file and inspecting discovered scaffold candidates.
 *
 * <p>This is intentionally a developer-facing validation tool rather than a polished user GUI.</p>
 */
public final class ScaffoldDiscoverySwingApp {
    private final JFrame frame;
    private final JTextField fileField;
    private final JButton browseButton;
    private final JButton runButton;
    private final JSpinner neighborCountSpinner;
    private final JSpinner minScaffoldSizeSpinner;
    private final JSpinner minSupportSpinner;
    private final JLabel statusLabel;
    private final JProgressBar progressBar;
    private final CandidateTableModel tableModel;
    private final JTable candidateTable;
    private final JStructureView structureView;
    private final JTextArea detailArea;

    private ScaffoldDiscoverySwingApp() {
        frame = new JFrame("structurized Scaffold Discovery");
        fileField = new JTextField(40);
        browseButton = new JButton("Browse");
        runButton = new JButton("Load And Run");
        neighborCountSpinner = new JSpinner(new SpinnerNumberModel(4, 1, 20, 1));
        minScaffoldSizeSpinner = new JSpinner(new SpinnerNumberModel(6, 1, 50, 1));
        minSupportSpinner = new JSpinner(new SpinnerNumberModel(2, 1, 1000, 1));
        statusLabel = new JLabel("Choose a SMILES file to start.");
        progressBar = new JProgressBar();
        tableModel = new CandidateTableModel();
        candidateTable = new JTable(tableModel);
        structureView = new JStructureView();
        detailArea = new JTextArea();
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception ignored) {
                // default look and feel is acceptable for this internal validation tool
            }
            new ScaffoldDiscoverySwingApp().show();
        });
    }

    private void show() {
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        frame.setLayout(new BorderLayout(8, 8));
        frame.add(buildTopPanel(), BorderLayout.NORTH);
        frame.add(buildCenterPanel(), BorderLayout.CENTER);
        frame.add(buildStatusPanel(), BorderLayout.SOUTH);
        frame.setSize(new Dimension(1300, 760));
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    private JPanel buildTopPanel() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBorder(BorderFactory.createEmptyBorder(8, 8, 0, 8));

        JPanel filePanel = new JPanel(new BorderLayout(8, 8));
        filePanel.add(new JLabel("Input SMILES File"), BorderLayout.WEST);
        filePanel.add(fileField, BorderLayout.CENTER);

        JPanel buttons = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        buttons.add(browseButton);
        buttons.add(runButton);
        filePanel.add(buttons, BorderLayout.EAST);

        JPanel options = new JPanel(new FlowLayout(FlowLayout.LEFT, 12, 0));
        options.add(new JLabel("Neighbors"));
        options.add(neighborCountSpinner);
        options.add(new JLabel("Min Scaffold Heavy Atoms"));
        options.add(minScaffoldSizeSpinner);
        options.add(new JLabel("Min Support"));
        options.add(minSupportSpinner);

        browseButton.addActionListener(event -> chooseFile());
        runButton.addActionListener(event -> runDiscovery());

        panel.add(filePanel, BorderLayout.NORTH);
        panel.add(options, BorderLayout.SOUTH);
        return panel;
    }

    private JSplitPane buildCenterPanel() {
        candidateTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        candidateTable.setAutoCreateRowSorter(true);
        candidateTable.setRowHeight(24);
        TableRowSorter<CandidateTableModel> sorter = new TableRowSorter<>(tableModel);
        candidateTable.setRowSorter(sorter);

        DefaultTableCellRenderer rightRenderer = new DefaultTableCellRenderer();
        rightRenderer.setHorizontalAlignment(JLabel.RIGHT);
        for (int column : new int[]{0, 1, 2, 3, 4, 5}) {
            candidateTable.getColumnModel().getColumn(column).setCellRenderer(rightRenderer);
        }

        candidateTable.getSelectionModel().addListSelectionListener(event -> {
            if (!event.getValueIsAdjusting()) {
                updateDetailView();
            }
        });

        JScrollPane tableScroll = new JScrollPane(candidateTable);

        structureView.setBorder(BorderFactory.createTitledBorder("Scaffold"));
        detailArea.setEditable(false);
        detailArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        detailArea.setLineWrap(true);
        detailArea.setWrapStyleWord(true);

        JPanel detailPanel = new JPanel(new BorderLayout(8, 8));
        detailPanel.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 8));
        detailPanel.add(structureView, BorderLayout.NORTH);
        detailPanel.add(new JScrollPane(detailArea), BorderLayout.CENTER);
        structureView.setPreferredSize(new Dimension(380, 280));

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, tableScroll, detailPanel);
        splitPane.setResizeWeight(0.65);
        return splitPane;
    }

    private JPanel buildStatusPanel() {
        JPanel panel = new JPanel(new BorderLayout(8, 8));
        panel.setBorder(BorderFactory.createEmptyBorder(0, 8, 8, 8));
        progressBar.setIndeterminate(false);
        progressBar.setStringPainted(false);
        panel.add(statusLabel, BorderLayout.CENTER);
        panel.add(progressBar, BorderLayout.EAST);
        return panel;
    }

    private void chooseFile() {
        JFileChooser chooser = new JFileChooser();
        if (!fileField.getText().isBlank()) {
            chooser.setSelectedFile(Path.of(fileField.getText()).toFile());
        }
        if (chooser.showOpenDialog(frame) == JFileChooser.APPROVE_OPTION) {
            fileField.setText(chooser.getSelectedFile().getAbsolutePath());
        }
    }

    private void runDiscovery() {
        if (fileField.getText().isBlank()) {
            JOptionPane.showMessageDialog(frame, "Please choose an input file first.", "No Input File", JOptionPane.WARNING_MESSAGE);
            return;
        }

        setBusy(true, "Loading input file...");
        Path path = Path.of(fileField.getText());
        ScaffoldDiscoveryConfig cfg = new ScaffoldDiscoveryConfig();
        cfg.neighborCount = (Integer) neighborCountSpinner.getValue();
        cfg.minScaffoldHeavyAtoms = (Integer) minScaffoldSizeSpinner.getValue();
        cfg.minSupport = (Integer) minSupportSpinner.getValue();

        SwingWorker<ScaffoldDiscoveryResult, Void> worker = new SwingWorker<>() {
            @Override
            protected ScaffoldDiscoveryResult doInBackground() throws Exception {
                List<String> smiles = SmilesInputReader.readSmilesFile(path);
                List<StereoMolecule> molecules = new ArrayList<>(smiles.size());
                SmilesParser parser = new SmilesParser();
                for (String smi : smiles) {
                    StereoMolecule molecule = new StereoMolecule();
                    parser.parse(molecule, smi);
                    molecule.ensureHelperArrays(StereoMolecule.cHelperSymmetrySimple);
                    molecules.add(molecule);
                }
                return ScaffoldDiscoveryEngine.discover(
                        molecules,
                        new TransformationBenchDemo.OCLMCSFastProvider(),
                        cfg
                );
            }

            @Override
            protected void done() {
                try {
                    ScaffoldDiscoveryResult result = get();
                    tableModel.setCandidates(result.candidates);
                    if (!result.candidates.isEmpty()) {
                        candidateTable.setRowSelectionInterval(0, 0);
                    } else {
                        structureView.structureChanged();
                        detailArea.setText("");
                    }
                    statusLabel.setText("Loaded " + result.compounds.size()
                            + " compounds, found " + result.candidates.size()
                            + " scaffold candidates from " + result.seedCount + " seeds.");
                } catch (Exception ex) {
                    tableModel.setCandidates(List.of());
                    detailArea.setText("");
                    JOptionPane.showMessageDialog(
                            frame,
                            "Failed to run scaffold discovery:\n" + ex.getMessage(),
                            "Discovery Failed",
                            JOptionPane.ERROR_MESSAGE
                    );
                    statusLabel.setText("Discovery failed.");
                } finally {
                    setBusy(false, statusLabel.getText());
                }
            }
        };
        worker.execute();
    }

    private void updateDetailView() {
        int viewRow = candidateTable.getSelectedRow();
        if (viewRow < 0) {
            structureView.structureChanged();
            detailArea.setText("");
            return;
        }

        int modelRow = candidateTable.convertRowIndexToModel(viewRow);
        ScaffoldCandidate candidate = tableModel.getCandidateAt(modelRow);
        structureView.setIDCode(candidate.template.idcode);

        String scaffoldSmiles = new IsomericSmilesCreator(candidate.template.scaffold).getSmiles();
        detailArea.setText("""
                Scaffold SMILES: %s
                Scaffold IDCode: %s
                Combined Score: %.2f
                Support Count: %d
                Average Explained Fraction: %.3f
                Scaffold Heavy Atom Count: %d
                Observed Exit Vector Count: %d
                Observed Exit Vector Atoms: %s
                Observed Exit Vector Symmetry Classes: %s
                Discovery Seeds: %s
                Discovery Compounds: %s
                Support Compounds: %s
                """.formatted(
                scaffoldSmiles,
                candidate.template.idcode,
                candidate.combinedScore,
                candidate.supportCount,
                candidate.averageExplainedFraction,
                candidate.scaffoldHeavyAtomCount,
                candidate.observedExitVectorCount,
                candidate.observedExitVectorAtoms,
                candidate.observedExitVectorSymmetryClasses,
                candidate.discoverySeedIndices,
                candidate.discoveryCompoundIndices,
                candidate.supportCompoundIndices
        ));
        detailArea.setCaretPosition(0);
    }

    private void setBusy(boolean busy, String message) {
        browseButton.setEnabled(!busy);
        runButton.setEnabled(!busy);
        neighborCountSpinner.setEnabled(!busy);
        minScaffoldSizeSpinner.setEnabled(!busy);
        minSupportSpinner.setEnabled(!busy);
        progressBar.setIndeterminate(busy);
        statusLabel.setText(message);
    }

    private static final class CandidateTableModel extends AbstractTableModel {
        private static final String[] COLUMNS = {
                "Rank", "Score", "Support", "Avg Explained", "Size", "ExitVecs", "Scaffold SMILES"
        };
        private final DecimalFormat scoreFormat = new DecimalFormat("0.00");
        private final DecimalFormat fractionFormat = new DecimalFormat("0.000");
        private List<ScaffoldCandidate> candidates = List.of();

        void setCandidates(List<ScaffoldCandidate> candidates) {
            this.candidates = List.copyOf(candidates);
            fireTableDataChanged();
        }

        ScaffoldCandidate getCandidateAt(int rowIndex) {
            return candidates.get(rowIndex);
        }

        @Override
        public int getRowCount() {
            return candidates.size();
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
        public Object getValueAt(int rowIndex, int columnIndex) {
            ScaffoldCandidate candidate = candidates.get(rowIndex);
            return switch (columnIndex) {
                case 0 -> rowIndex + 1;
                case 1 -> scoreFormat.format(candidate.combinedScore);
                case 2 -> candidate.supportCount;
                case 3 -> fractionFormat.format(candidate.averageExplainedFraction);
                case 4 -> candidate.scaffoldHeavyAtomCount;
                case 5 -> candidate.observedExitVectorCount;
                case 6 -> new IsomericSmilesCreator(candidate.template.scaffold).getSmiles();
                default -> "";
            };
        }
    }
}

package br.unifesp.fisicamedica;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagLayout;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import ij.ImagePlus;
import net.imagej.ImageJ;

@Plugin(type = Command.class, menuPath = "Plugins>Detectabilidade")
public class DetectabilidadePlugin implements Command {

    private static JFrame frame;
    private JTree tree;
    private DefaultTreeModel treeModel;
    
    private DefaultMutableTreeNode root;
    private DicomProcessor dicomProcessor;
    private JPanel leftSidePanel;
    private CardLayout cardLayout;

    public DetectabilidadePlugin() {
        this.root = new DefaultMutableTreeNode("");
        this.treeModel = new DefaultTreeModel(root);
        this.dicomProcessor = new DicomProcessor();
    }

    @Override
    public void run() {
    	try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (ClassNotFoundException | InstantiationException |
                 IllegalAccessException | UnsupportedLookAndFeelException e) {
            e.printStackTrace();
        }
        createAndShowGUI();
    }

    private void createAndShowGUI() {
        frame = new JFrame("");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setSize(new Dimension(1300,700));

        JPanel mainPanel = new JPanel(new BorderLayout());
        
        cardLayout = new CardLayout();
        leftSidePanel = new JPanel(cardLayout);
        leftSidePanel.setBackground(Color.WHITE);
        leftSidePanel.setBorder(BorderFactory.createLineBorder(new Color(173, 216, 230), 2));

        JPanel buttonPanel = new JPanel(new GridBagLayout());
        buttonPanel.setBackground(Color.WHITE);
        JButton selectDirButton = new JButton("Select Directory");
        selectDirButton.addActionListener(e -> selectAndAnalyzeDirectory());
        buttonPanel.add(selectDirButton);

        leftSidePanel.add(buttonPanel, "BUTTON");
        leftSidePanel.add(createTreePanel(), "TREE");

        cardLayout.show(leftSidePanel, "BUTTON");

        JPanel rightPanel = TopPanel.create();

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftSidePanel, rightPanel);
        splitPane.setDividerLocation(300);

        mainPanel.add(splitPane, BorderLayout.CENTER);
        frame.add(mainPanel);
        frame.setVisible(true);
    }

    private JPanel createTreePanel() {
        tree = new JTree(treeModel);
        tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        ImageOpener imageOpener = new ImageOpener();

        tree.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    TreePath path = tree.getPathForLocation(e.getX(), e.getY());
                    if (path != null) {
                        DefaultMutableTreeNode selectedNode = (DefaultMutableTreeNode) path.getLastPathComponent();
                        if (selectedNode.isLeaf()) {
                            Object userObject = selectedNode.getUserObject();
                            if (userObject instanceof SeriesData) {
                                SeriesData serieNode = (SeriesData) userObject;
                                ImagePlus stackImage = imageOpener.createImageStack(serieNode);
                                if (stackImage != null) {
                                    stackImage.show();
                                }
                            }
                        }
                    }
                }
            }
        });

        JScrollPane treeView = new JScrollPane(tree);
        JPanel treePanel = new JPanel(new BorderLayout());
        treePanel.add(treeView, BorderLayout.CENTER);
        return treePanel;
    }

    private void selectAndAnalyzeDirectory() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        if (fileChooser.showOpenDialog(frame) == JFileChooser.APPROVE_OPTION) {
            File dir = fileChooser.getSelectedFile();
            
            JDialog loadingDialog = createLoadingDialog("Loading images...");
            new Thread(() -> {
                try {
                    dicomProcessor.analyzeDirectory(dir, loadingDialog);
                    updateTree(dicomProcessor.getPatientData());
                    SwingUtilities.invokeLater(() -> loadingDialog.dispose());
                    cardLayout.show(leftSidePanel, "TREE");
                } catch (Exception e) {
                    e.printStackTrace();
                    SwingUtilities.invokeLater(() -> {
                        loadingDialog.dispose();
                        JOptionPane.showMessageDialog(frame, "Error loading directory: " + e.getMessage(),
                                "Error", JOptionPane.ERROR_MESSAGE);
                    });
                }
            }).start();
        }
    }

    public static JDialog createLoadingDialog(String message) {
        JDialog dialog = new JDialog(frame, "", true);
        dialog.setSize(300, 100);
        dialog.setLocationRelativeTo(frame);
        dialog.setLayout(new BorderLayout());

        JLabel label = new JLabel(message, JLabel.CENTER);
        JProgressBar progressBar = new JProgressBar();
        progressBar.setStringPainted(true);
        progressBar.setValue(0);

        dialog.add(label, BorderLayout.NORTH);
        dialog.add(progressBar, BorderLayout.CENTER);

        SwingUtilities.invokeLater(() -> dialog.setVisible(true));
        return dialog;
    }

    private void updateTree(Map<String, Map<String, Map<String, Map<String, SeriesData>>>> patientData) {
        root.removeAllChildren();
        for (String patient : patientData.keySet()) {
            DefaultMutableTreeNode patientNode = new DefaultMutableTreeNode(patient);
            root.add(patientNode);
            for (String study : patientData.get(patient).keySet()) {
                DefaultMutableTreeNode studyNode = new DefaultMutableTreeNode(study);
                patientNode.add(studyNode);
                for (String protocol : patientData.get(patient).get(study).keySet()) {
                    DefaultMutableTreeNode protocolNode = new DefaultMutableTreeNode(protocol);
                    studyNode.add(protocolNode);
                    for (SeriesData series : patientData.get(patient).get(study).get(protocol).values()) {
                        DefaultMutableTreeNode seriesNode = new DefaultMutableTreeNode(series);
                        protocolNode.add(seriesNode);
                        series.sortFileMetadataByInstanceNumber();
                    }
                }
            }
        }
        treeModel.reload();
    }

    public static void main(final String... args) throws Exception {
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        ij.command().run(DetectabilidadePlugin.class, true);
    }
}
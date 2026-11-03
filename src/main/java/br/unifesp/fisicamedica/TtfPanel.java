package br.unifesp.fisicamedica;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.jfree.chart.ChartPanel;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;

public class TtfPanel {

    private static JTextField startSliceField;
    private static JTextField endSliceField;
    
    private static JComboBox<String> detectionModeCombo;
    private static JCheckBox avgSlicesCheck;
    
    private static ChartPanel chartPanel;
    private static JComboBox<String> curveTypeCombo;
    
    private static JPanel taskChecksPanel;
    private static List<JCheckBox> insertCheckBoxes;
    
    private static Segmenter.DetectorResult lastDetectionResult = null;
    private static ImageStack lastStack = null;
    private static Calibration lastCal = null;
    
    private static TtfCalc.mtfResult accumulatedTtfResult = new TtfCalc.mtfResult();
    private static List<Integer> accumulatedTaskIDs = new ArrayList<>();
    private static int nextTaskId = 0; 

    public static JPanel create() {
        JPanel root = new JPanel(new GridBagLayout());
        GridBagConstraints g = new GridBagConstraints();
        g.fill = GridBagConstraints.BOTH;
        g.insets = new Insets(8, 8, 8, 8);

        JPanel params = new JPanel(new GridBagLayout());
        params.setPreferredSize(new Dimension(350, 320));
        params.setBorder(BorderFactory.createTitledBorder("TTF Parameters"));
        
        GridBagConstraints pGBC = new GridBagConstraints();
        pGBC.insets = new Insets(5, 5, 5, 5);
        pGBC.fill = GridBagConstraints.HORIZONTAL;

        int row = 0;
        
        pGBC.gridx = 0; pGBC.gridy = row; pGBC.gridwidth = 4;
        params.add(new JLabel("Slice selection:"), pGBC);
        row++;

        pGBC.gridwidth = 1;
        pGBC.gridx = 0; pGBC.gridy = row;
        params.add(new JLabel("Start slice:"), pGBC);
        startSliceField = new JTextField("1", 3); 
        pGBC.gridx = 1; pGBC.gridy = row;
        params.add(startSliceField, pGBC);

        pGBC.gridx = 2; pGBC.gridy = row;
        params.add(new JLabel("End slice:"), pGBC);
        endSliceField = new JTextField("1", 3); 
        pGBC.gridx = 3; pGBC.gridy = row;
        params.add(endSliceField, pGBC);
        row++;

        pGBC.gridx = 0; pGBC.gridy = row; pGBC.gridwidth = 4;
        params.add(new JLabel("Detection Method:"), pGBC);
        row++;
        
        String[] detectionOptions = {"Automatic", "ROI"};
        detectionModeCombo = new JComboBox<>(detectionOptions);
        pGBC.gridx = 0; pGBC.gridy = row; pGBC.gridwidth = 4;
        params.add(detectionModeCombo, pGBC);
        row++;

        avgSlicesCheck = new JCheckBox("Use average slice", true); 
        pGBC.gridx = 0; pGBC.gridy = row; pGBC.gridwidth = 4;
        params.add(avgSlicesCheck, pGBC);
        row++;

        JButton detectButton = new JButton("Find Edges");
        pGBC.gridx = 0; pGBC.gridy = row; pGBC.gridwidth = 4;
        params.add(detectButton, pGBC);
        row++;
        
        JButton calcButton = new JButton("Calculate TTF");
        pGBC.gridx = 0; pGBC.gridy = row; pGBC.gridwidth = 4;
        params.add(calcButton, pGBC);
        row++;

        pGBC.gridx = 0; pGBC.gridy = row; pGBC.gridwidth = 4;
        pGBC.weighty = 1.0; 
        params.add(Box.createVerticalGlue(), pGBC);

        g.gridx = 0; g.gridy = 0; g.weightx = 0.3; g.weighty = 1.0;
        root.add(params, g);
        
        JPanel results = new JPanel(new GridBagLayout());
        results.setBorder(BorderFactory.createTitledBorder("TTF Results"));
        GridBagConstraints rGBC = new GridBagConstraints();
        rGBC.insets = new Insets(5, 5, 5, 5);

        JPanel buttonsPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 10, 0));
        String[] curveOptions = {"TTF", "LSF", "ESF"};
        curveTypeCombo = new JComboBox<>(curveOptions);
        
        JButton showTableButton = new JButton("Show Data Table");
        JButton clearAllButton = new JButton("Clear All");
        
        buttonsPanel.add(new JLabel("Show Curve:"));
        buttonsPanel.add(curveTypeCombo);
        buttonsPanel.add(showTableButton);
        buttonsPanel.add(clearAllButton);
        
        taskChecksPanel = new JPanel(new GridBagLayout());
        insertCheckBoxes = new ArrayList<>();
        
        rGBC.gridx = 0; rGBC.gridy = 0; rGBC.fill = GridBagConstraints.HORIZONTAL; rGBC.weightx = 1.0;
        results.add(buttonsPanel, rGBC);
        
        rGBC.gridy = 1;
        results.add(taskChecksPanel, rGBC);
        
        chartPanel = TtfGraphManager.createChartPanel();
        rGBC.gridy = 2; rGBC.fill = GridBagConstraints.BOTH; rGBC.weighty = 1.0; 
        results.add(chartPanel, rGBC);

        g.gridx = 1; g.gridy = 0; g.weightx = 0.7; g.weighty = 1.0;
        root.add(results, g);

        detectButton.addActionListener(e -> detectInserts());
        calcButton.addActionListener(e -> calculateTtf());
        showTableButton.addActionListener(e -> openTableWindow());
        clearAllButton.addActionListener(e -> clearAllTasks());
        curveTypeCombo.addActionListener(e -> updateGraph());
        
        return root;
    }

    private static void detectInserts() {
        try {
            int startSlice = Integer.parseInt(startSliceField.getText().trim());
            int endSlice = Integer.parseInt(endSliceField.getText().trim());

            if (endSlice < startSlice) {
                int temp = endSlice;
                endSlice = startSlice;
                startSlice = temp;
            }

            SeriesData sd = ImageOpener.getActiveSeriesData();
            if (sd == null) {
                IJ.log("No active series. Use the file browser to open images");
                return;
            }
            
            List<File> files = sd.getFiles();
            ImageStack stack = new ImageStack();
            Calibration cal = null;

            for (int i = startSlice - 1; i < endSlice; i++) {
                File file = files.get(i);
                ImagePlus tempImp = new ImagePlus(file.getAbsolutePath());
                if (tempImp.getWidth() == 0) continue;
                if (cal == null) cal = tempImp.getCalibration();
                stack.addSlice(tempImp.getProcessor());
            }

            String detectionMode = getDetectionMode();
            boolean isAuto = detectionMode.equals("Automatic");
            boolean avgSlices = avgSlicesCheck.isSelected();
            
            Roi userRoi = null;
            if (!isAuto) {
                ImagePlus activeImage = IJ.getImage();
                if (activeImage != null) {
                    userRoi = activeImage.getRoi();
                }
                if (userRoi == null) {
                    IJ.error("Please draw a ROI on the image first, or select 'Automatic' detection.");
                    return;
                }
            }
            
            lastDetectionResult = Segmenter.edgeFinder(stack, cal, userRoi, isAuto, avgSlices);
            
            lastStack = stack; 
            lastCal = cal;

            if (lastDetectionResult == null || lastDetectionResult.rois.isEmpty()) {
                JOptionPane.showMessageDialog(null, "Could not segment edges.", "", JOptionPane.WARNING_MESSAGE);
            } else {
                RoiManager manager = RoiManager.getInstance();
                if(manager==null) manager = new RoiManager();
                
                for(List<Roi> list: lastDetectionResult.rois) {
                    for(int i=0; i<list.size();i++) {
                        Roi r = list.get(i);
                        r.setPosition(startSlice+i);
                        manager.addRoi(r);
                    }
                }
            } 
            
        } catch (NumberFormatException ex) {
            IJ.error("Invalid slice parameters.");
        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }    
    

    private static void calculateTtf() {
        if (lastDetectionResult == null || lastDetectionResult.rois.isEmpty()) {
            return;
        }

        try {
            TtfCalc.mtfResult newResult = TtfCalc.analisarStackComInserts(lastStack, lastCal, lastDetectionResult, 128);

            if (newResult == null || newResult.ttfs.isEmpty()) {
                IJ.error("No TTF calculated for this region.");
                return;
            }

            for (Integer oldId : newResult.ttfs.keySet()) {
                int newId = nextTaskId++;
                accumulatedTtfResult.ttfs.put(newId, newResult.ttfs.get(oldId));
                accumulatedTtfResult.freqs.put(newId, newResult.freqs.get(oldId));
                accumulatedTtfResult.lsfs.put(newId, newResult.lsfs.get(oldId));
                accumulatedTtfResult.esfs.put(newId, newResult.esfs.get(oldId));
                accumulatedTtfResult.contrasts.put(newId, newResult.contrasts.get(oldId));
                accumulatedTtfResult.I0.put(newId, newResult.I0.get(oldId));
                accumulatedTtfResult.I1.put(newId, newResult.I1.get(oldId));
                
                accumulatedTaskIDs.add(newId);
            }

            rebuildTaskPanel();

        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }

    private static void removeTask(int id) {
        accumulatedTtfResult.ttfs.remove(id);
        accumulatedTtfResult.freqs.remove(id);
        accumulatedTtfResult.lsfs.remove(id);
        accumulatedTtfResult.esfs.remove(id);
        accumulatedTtfResult.contrasts.remove(id);
        accumulatedTtfResult.I0.remove(id);
        accumulatedTtfResult.I1.remove(id);
        accumulatedTaskIDs.remove(Integer.valueOf(id));
        
        rebuildTaskPanel();
    }

    private static void clearAllTasks() {
        accumulatedTtfResult = new TtfCalc.mtfResult();
        accumulatedTaskIDs.clear();
        nextTaskId = 0;
        
        rebuildTaskPanel();
    }
    
    private static void rebuildTaskPanel() {
        taskChecksPanel.removeAll();
        insertCheckBoxes.clear();
        
        GridBagConstraints cbGbc = new GridBagConstraints();
        cbGbc.insets = new Insets(2, 5, 2, 10);
        cbGbc.anchor = GridBagConstraints.WEST;
        int row = 0;
        
        DecimalFormat df = new DecimalFormat("0.00");

        for (Integer id : accumulatedTaskIDs) {
            JCheckBox cb = new JCheckBox("Task " + (id + 1), true);
            cb.addActionListener(e -> updateGraph());
            insertCheckBoxes.add(cb);
            
            cbGbc.gridx = 0; cbGbc.gridy = row; cbGbc.weightx = 0;
            taskChecksPanel.add(cb, cbGbc);
            
            Double contrast = accumulatedTtfResult.contrasts.get(id);
            Double I0 = accumulatedTtfResult.I0.get(id);
            Double I1 = accumulatedTtfResult.I1.get(id);

            String strContrast = "N/A";
            String strI0 = "N/A";
            String strI1 = "N/A";
            
            if (contrast != null) {
                strContrast = df.format(contrast) + " HU";
                strI0 = df.format(I0) + " HU";
                strI1 = df.format(I1) + " HU";
            }
            
            JLabel statsLabel = new JLabel("I0: "+strI0+" | I1: "+strI1+" | Contrast: " + strContrast);
            statsLabel.setFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 12));
            
            cbGbc.gridx = 1; cbGbc.gridy = row; cbGbc.weightx = 1.0;
            taskChecksPanel.add(statsLabel, cbGbc);
            
            // --- BOTÃO DE REMOÇÃO ATUALIZADO ---
            JButton removeBtn = new JButton("Clear");
            removeBtn.setToolTipText("Remove this task from the list");
            // Margens reduzidas e fonte menor para não estufar a linha
            removeBtn.setMargin(new Insets(2, 8, 2, 8));
            removeBtn.setFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 11));
            removeBtn.addActionListener(e -> removeTask(id));

            cbGbc.gridx = 2; cbGbc.gridy = row; cbGbc.weightx = 0;
            taskChecksPanel.add(removeBtn, cbGbc);
            
            row++;
        }
        
        taskChecksPanel.revalidate();
        taskChecksPanel.repaint();
        updateGraph();
    }
    
    private static void updateGraph() {
        if (accumulatedTtfResult == null || insertCheckBoxes.isEmpty()) {
            TtfGraphManager.clearChart();
            return;
        }
        
        List<Integer> selectedIDs = new ArrayList<>();
        for (JCheckBox cb : insertCheckBoxes) {
            if (cb.isSelected()) {
                String text = cb.getText().replace("Task ", "");
                selectedIDs.add(Integer.parseInt(text) - 1);
            }
        }
        
        String selectedCurve = (String) curveTypeCombo.getSelectedItem();
        TtfGraphManager.updateChart(selectedIDs, selectedCurve, accumulatedTtfResult);
    }

    private static void openTableWindow() {
        if (accumulatedTtfResult == null || accumulatedTtfResult.ttfs == null || accumulatedTaskIDs == null || accumulatedTaskIDs.isEmpty()) {
            return;
        }

        String selectedCurve = (String) curveTypeCombo.getSelectedItem();
        List<String> columnNames = new ArrayList<>();
        
        boolean isTtf = selectedCurve.equals("TTF");
        columnNames.add(isTtf ? "Freq (cycles/mm)" : "Position (Points)");
        
        java.util.List<double[]> dataList = new java.util.ArrayList<>();
        
        Integer firstId = accumulatedTaskIDs.get(0);
        
        if (isTtf) {
            dataList.add(accumulatedTtfResult.freqs.get(firstId));
        } else {
            // Cria um eixo X sintético para as curvas do domínio espacial (ESF/LSF) centrado em zero
            double[] referenceArray = selectedCurve.equals("ESF") ? 
                                      accumulatedTtfResult.esfs.get(firstId) : 
                                      accumulatedTtfResult.lsfs.get(firstId);
            double[] xAxis = new double[referenceArray.length];
            int center = referenceArray.length / 2;
            for (int i = 0; i < xAxis.length; i++) {
                xAxis[i] = i - center; 
            }
            dataList.add(xAxis);
        }

        for (Integer id : accumulatedTaskIDs) {
            double[] curve = null;
            if (selectedCurve.equals("TTF")) curve = accumulatedTtfResult.ttfs.get(id);
            else if (selectedCurve.equals("LSF")) curve = accumulatedTtfResult.lsfs.get(id);
            else if (selectedCurve.equals("ESF")) curve = accumulatedTtfResult.esfs.get(id);
            
            Double contrast = accumulatedTtfResult.contrasts.get(id);
            
            if (curve != null) {
                String colName = "Task " + (id + 1) + " " + selectedCurve;
                if (contrast != null) {
                    colName += " (\u0394HU " + String.format("%.0f", contrast) + ")";
                }
                columnNames.add(colName);
                dataList.add(curve);
            }
        }

        String[] columns = columnNames.toArray(new String[0]);
        DataTableFrame tableFrame = new DataTableFrame("1D " + selectedCurve + " Data", columns, dataList);
        tableFrame.setVisible(true);
    }
    
    public static String getDetectionMode() {
        return (String) detectionModeCombo.getSelectedItem();
    }
    
    public static Calibration getLastCal() { return lastCal; }
    public static List<Integer> getLastInsertIDs() { return accumulatedTaskIDs; }
    public static TtfCalc.mtfResult getLastTtfResult() { return accumulatedTtfResult; }
}
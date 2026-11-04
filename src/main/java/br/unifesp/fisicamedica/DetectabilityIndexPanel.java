package br.unifesp.fisicamedica;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.text.DecimalFormat;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;

import org.apache.commons.math3.special.Erf;

import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

public class DetectabilityIndexPanel {

    private static JComboBox<String> TaskSelector;
    private static JTextField taskContrastField;
    private static JTextField taskRadiusField;
    
    private static JTextField viewDistanceField;
    private static JTextField displaySizeField;
    private static JTextField fovField;
    
    private static JTextField dPrimeResultField;
    private static JTextField AucResultField;
    
    private static JButton btnShowImages;
    
    private static double[][] lastNps2D = null;
    private static double[] lastNpsFreq = null;
    private static double[] lastNpsData = null;
    private static double[] lastTtfFreq = null; 
    private static double[] lastTtfData = null; 
    private static double lastTaskContrast = 0;
    private static double lastTaskRadius = 0;
    
    private static final DecimalFormat DF = new DecimalFormat("0.####");

    public static JPanel create() {

        JPanel root = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(8, 8, 8, 8);
        gbc.fill = GridBagConstraints.HORIZONTAL;

        JPanel params = new JPanel(new GridBagLayout());
        params.setBorder(BorderFactory.createTitledBorder("d' Parameters"));
        GridBagConstraints pGBC = new GridBagConstraints();
        pGBC.insets = new Insets(5, 5, 5, 5);
        pGBC.fill = GridBagConstraints.HORIZONTAL;

        int row = 0;

        pGBC.gridx = 0; pGBC.gridy = row;
        params.add(new JLabel("Target Insert:"), pGBC);
        
        TaskSelector = new JComboBox<>();
        pGBC.gridx = 1; pGBC.gridy = row++;
        params.add(TaskSelector, pGBC);

        JButton calcBtn = new JButton("Calculate d'");
        pGBC.gridx = 0; pGBC.gridy = row++; pGBC.gridwidth = 2;
        params.add(calcBtn, pGBC);

        pGBC.gridwidth = 1;
        pGBC.gridx = 0; pGBC.gridy = row;
        params.add(new JLabel("Task Contrast (HU):"), pGBC);
        taskContrastField = new JTextField("7", 8);
        pGBC.gridx = 1; pGBC.gridy = row;
        params.add(taskContrastField, pGBC);
        row++;
        
        pGBC.gridx = 0; pGBC.gridy = row;
        params.add(new JLabel("Task Radius (mm):"), pGBC);
        taskRadiusField = new JTextField("6", 8);
        pGBC.gridx = 1; pGBC.gridy = row;
        params.add(taskRadiusField, pGBC);
        row++;
        
        pGBC.gridx = 0; pGBC.gridy = row;
        params.add(new JLabel("Viewing Distance (mm):"), pGBC);
        viewDistanceField = new JTextField("500", 8);
        pGBC.gridx = 1; pGBC.gridy = row;
        params.add(viewDistanceField, pGBC);
        row++;
        
        pGBC.gridx = 0; pGBC.gridy = row;
        params.add(new JLabel("Display Size (mm):"), pGBC);
        displaySizeField = new JTextField("305", 8);
        pGBC.gridx = 1; pGBC.gridy = row;
        params.add(displaySizeField, pGBC);
        row++;

        pGBC.gridx = 0; pGBC.gridy = row;
        params.add(new JLabel("Display FOV (mm):"), pGBC);
        fovField = new JTextField("250", 8);
        pGBC.gridx = 1; pGBC.gridy = row;
        params.add(fovField, pGBC);
        row++;
        
        updateDataSources();        
        
        JPanel results = new JPanel(new GridBagLayout());
        results.setBorder(BorderFactory.createTitledBorder("Result"));
        GridBagConstraints rGBC = new GridBagConstraints();
        rGBC.insets = new Insets(5, 5, 5, 5);
        
        rGBC.gridx = 0; rGBC.gridy = 0;
        results.add(new JLabel("Detectability Index (d'):"), rGBC);
        dPrimeResultField = new JTextField("-", 12);
        dPrimeResultField.setEditable(false);
        rGBC.gridx = 1; rGBC.gridy = 0;
        results.add(dPrimeResultField, rGBC);
        
        rGBC.gridx = 0; rGBC.gridy = 1;
        results.add(new JLabel("AUC:"), rGBC);
        AucResultField = new JTextField("-", 12);
        AucResultField.setEditable(false);
        rGBC.gridx = 1; rGBC.gridy = 1;
        results.add(AucResultField, rGBC);

        btnShowImages = new JButton("Show Detectability Task");
        btnShowImages.setEnabled(false); 
        rGBC.gridx = 0; rGBC.gridy = 2; rGBC.gridwidth = 2;
        rGBC.fill = GridBagConstraints.HORIZONTAL;
        results.add(btnShowImages, rGBC);

        gbc.gridx = 0; gbc.gridy = 0; gbc.weightx = 1;
        params.setPreferredSize(new Dimension(350, 320));
        root.add(params, gbc);
        
        gbc.gridx = 0; gbc.gridy = 1;
        results.setPreferredSize(new Dimension(350, 160)); 
        root.add(results, gbc);

        TaskSelector.addPopupMenuListener(new PopupMenuListener() {
            @Override
            public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
                updateDataSources();
            }
            @Override
            public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {}
            @Override
            public void popupMenuCanceled(PopupMenuEvent e) {}
        });

        calcBtn.addActionListener(e -> calculate(root));
        
        btnShowImages.addActionListener(e -> {
            if (lastNps2D != null && lastTtfFreq != null && lastTtfData != null) {
               Calibration realCal = TtfPanel.getLastCal();
               double realPixelSize = (realCal != null && realCal.pixelWidth > 0) ? realCal.pixelWidth : 0.5;
               exibirStackNoImageJ(realPixelSize, 512, lastNps2D, lastTtfFreq, lastTtfData, lastTaskContrast, lastTaskRadius);
            }
        });
        

        return root;
    }

    private static void updateDataSources() {
        List<Integer> insertIDs = TtfPanel.getLastInsertIDs();

        if (insertIDs == null || insertIDs.isEmpty()) {
            TaskSelector.setModel(new DefaultComboBoxModel<>());
            return;
        }

        Object currentSelection = TaskSelector.getSelectedItem();
        
        DefaultComboBoxModel<String> newModel = new DefaultComboBoxModel<>();
        for (Integer id : insertIDs) {
            newModel.addElement("TTF Task " + (id + 1));
        }
        TaskSelector.setModel(newModel);
        
        if (currentSelection != null && newModel.getIndexOf(currentSelection) >= 0) {
            TaskSelector.setSelectedItem(currentSelection);
        }
    }
    
    private static void calculate(JPanel parent) {
        btnShowImages.setEnabled(false);
        
        try {
            double[] npsFreq = NPSpanel.getLastRadialFreq();
            double[] npsData = NPSpanel.getLastRadialNPS();
            double[][] nps2D = NPSpanel.getLastNPS2D();
            
            if (npsFreq == null || npsData == null || nps2D == null) {
                JOptionPane.showMessageDialog(parent, "Error: No NPS calculated.", "Erro", JOptionPane.ERROR_MESSAGE);
                return;
            }

            TtfCalc.mtfResult ttfRes = TtfPanel.getLastTtfResult();
            String selectedStr = (String) TaskSelector.getSelectedItem();
            
            if (selectedStr == null || ttfRes == null || ttfRes.ttfs == null) {
                 JOptionPane.showMessageDialog(parent, "No valid insert selected or no TTF data available.", "Erro", JOptionPane.ERROR_MESSAGE);
                return;
            }

            int selectedID = Integer.parseInt(selectedStr.replace("TTF Task ", "")) - 1;

            if (!ttfRes.ttfs.containsKey(selectedID)) {
                 JOptionPane.showMessageDialog(parent, "No valid TTF data for this insert.", "Erro", JOptionPane.ERROR_MESSAGE);
                 return;
            }

            double[] ttfData = ttfRes.ttfs.get(selectedID);
            double[] ttfFreq = ttfRes.freqs.get(selectedID);
            
            if (ttfData == null || ttfFreq == null) {
                JOptionPane.showMessageDialog(parent, "Nenhum dado de TTF encontrado para este inserto.", "Erro", JOptionPane.ERROR_MESSAGE);
                return;
            }
            
            double taskContrast = Double.parseDouble(taskContrastField.getText().trim());
            double taskRadius = Double.parseDouble(taskRadiusField.getText().trim());
            double viewDistance = Double.parseDouble(viewDistanceField.getText().trim());
            double displaySize = Double.parseDouble(displaySizeField.getText().trim());
            double fov = Double.parseDouble(fovField.getText().trim());

            DetectabilityIndexCalculator calculator = new DetectabilityIndexCalculator(
                ttfData, ttfFreq, 
                npsData, npsFreq, 
                taskRadius, 
                taskContrast, 
                viewDistance, 
                displaySize,
                fov
            );
            
            double dPrime = calculator.compute();

            dPrimeResultField.setText(DF.format(dPrime));
            double auc = 0.5 + 0.5 * Erf.erf(dPrime / 2.0);
            AucResultField.setText(DF.format(auc));

            lastNps2D = nps2D;
            lastNpsFreq = npsFreq;
            lastNpsData = npsData;
            lastTtfFreq = ttfFreq;
            lastTtfData = ttfData;
            lastTaskContrast = taskContrast;
            lastTaskRadius = taskRadius;
            
            btnShowImages.setEnabled(true);

        } catch (NumberFormatException ex) {
            JOptionPane.showMessageDialog(parent, "Error in numerical input", "Error", JOptionPane.ERROR_MESSAGE);
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(parent, "Error: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
            ex.printStackTrace();
        }
    }
    
    private static double[][] expandNpsToSize(double[][] smallNps, int targetSize) {
        int h = smallNps.length;
        int w = smallNps[0].length;
        FloatProcessor fp = new FloatProcessor(w, h);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                fp.setf(x, y, (float) smallNps[y][x]);
            }
        }
        fp.setInterpolationMethod(ImageProcessor.BICUBIC);
        FloatProcessor resized = (FloatProcessor) fp.resize(targetSize, targetSize);
        
        double[][] bigNps = new double[targetSize][targetSize];
        for (int y = 0; y < targetSize; y++) {
            for (int x = 0; x < targetSize; x++) {
                bigNps[y][x] = resized.getf(x, y);
            }
        }
        return bigNps;
    }

    private static void exibirStackNoImageJ(double pixelSize, int matrizSize, double[][] nps2D, double[] ttfFreq, double[] ttfData, double contrast, double radius) {
       
        FloatProcessor ipTemplate = new FloatProcessor(matrizSize, matrizSize);
        double rPixels = radius / pixelSize;
        int center = matrizSize / 2;
        int subSamples = 4; 
        
        for (int y = 0; y < matrizSize; y++) {
            for (int x = 0; x < matrizSize; x++) {
                int insideCount = 0;
                
                for (int sy = 0; sy < subSamples; sy++) {
                    for (int sx = 0; sx < subSamples; sx++) {
                        double subX = (x - center) + (sx + 0.5) / subSamples - 0.5;
                        double subY = (y - center) + (sy + 0.5) / subSamples - 0.5;
                        
                        if (Math.sqrt(subX * subX + subY * subY) <= rPixels) {
                            insideCount++;
                        }
                    }
                }
                
                double coverage = (double) insideCount / (subSamples * subSamples);
                ipTemplate.putPixelValue(x, y, contrast * coverage);
            }
        }

        double[][] mtf2D = new double[matrizSize][matrizSize];
        double df2D = 1.0 / (matrizSize * pixelSize);
        
        for (int y = 0; y < matrizSize; y++) {
            for (int x = 0; x < matrizSize; x++) {
                double fx = (x - center) * df2D;
                double fy = (y - center) * df2D;
                double f = Math.sqrt(fx * fx + fy * fy);
                mtf2D[y][x] = interpolate1D(ttfFreq, ttfData, f);
            }
        }

        ImageProcessor ipMtf = MTFModel.applyMTF(ipTemplate, mtf2D);
        NPSModel npsModel = new NPSModel();
        
        double[][] expandedNps2D = expandNpsToSize(nps2D, matrizSize);
        ImageProcessor ipRuido = npsModel.texturedNoiseIP(expandedNps2D, pixelSize);
        
        FloatProcessor ipFinal = (FloatProcessor) ipMtf.duplicate();
        for (int i = 0; i < ipFinal.getPixelCount(); i++) {
            ipFinal.setf(i, ipFinal.getf(i) + ipRuido.getf(i));
        }

        ImageStack stack = new ImageStack(matrizSize, matrizSize);
        stack.addSlice("1 - Task", ipTemplate);
        stack.addSlice("2 - Noise", ipRuido);
        stack.addSlice("3 - Task*TTF + NPS", ipFinal);

        ImagePlus imp = new ImagePlus("Detectability Task", stack);
        
        double w = Math.max(contrast * 4, 30);
        double l = contrast / 2.0;
        imp.setDisplayRange(l - w/2.0, l + w/2.0);
        
        imp.show();
    }

    private static double interpolate1D(double[] freq, double[] data, double targetF) {
        if (freq == null || freq.length == 0) return 0.0;
        if (targetF >= freq[freq.length - 1]) return 0.0;
        if (targetF <= freq[0]) return data[0];
        
        for (int j = 0; j < freq.length - 1; j++) {
            if (targetF >= freq[j] && targetF <= freq[j + 1]) {
                double f0 = freq[j];
                double f1 = freq[j + 1];
                double m0 = data[j];
                double m1 = data[j + 1];
                return m0 + (targetF - f0) * (m1 - m0) / (f1 - f0);
            }
        }
        return 0.0;
    }
}
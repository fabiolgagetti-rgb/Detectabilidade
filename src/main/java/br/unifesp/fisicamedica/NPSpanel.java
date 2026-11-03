package br.unifesp.fisicamedica;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Rectangle;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
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

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class NPSpanel {

    private static JTextField startSliceField;
    private static JTextField endSliceField;
    private static JTextField recordPixelSizeField;
    private static JTextField recordNumberField;
    private static JTextField recordRadiusField;
    private static JCheckBox detrendCheck;

    private static JComboBox<String> positioningCombo;
    private static JButton btnAddRoi;
    private static JButton btnClearRois;
    private static JLabel lblRoiCount;

    private static JTextField varianceField;
    private static JTextField integratedField;
    private static JTextField peakFreqField;
    private static JTextField meanFreqField;

    private static double[][] lastNPS2D;         
    private static double[] lastRadialFreq;      
    private static double[] lastRadialNPS;       
    private static final DecimalFormat DF = new DecimalFormat("0.######");

    private static List<Roi> userRois = new ArrayList<>();

    public static JPanel create() {
        JPanel root = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(6, 6, 6, 6);
        gbc.fill = GridBagConstraints.BOTH;

        JPanel params = new JPanel(new GridBagLayout());
        params.setBorder(BorderFactory.createTitledBorder("NPS Parameters"));

        int row = 0;
        
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 1;
        params.add(new JLabel("Slice selection:"), gbc);
        row++;

        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 1;
        params.add(new JLabel("Start slice:"), gbc);
        startSliceField = new JTextField("1", 3);
        gbc.gridx = 1; gbc.gridy = row; gbc.gridwidth = 1;
        params.add(startSliceField, gbc);

        gbc.gridx = 2; gbc.gridy = row; gbc.gridwidth = 1;
        params.add(new JLabel("End slice:"), gbc);
        endSliceField = new JTextField("1", 3);
        gbc.gridx = 3; gbc.gridy = row; gbc.gridwidth = 1;
        params.add(endSliceField, gbc);
        row++;
        
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 2;
        params.add(new JLabel("ROIs positioning:"), gbc);
        
        positioningCombo = new JComboBox<>(new String[]{"Automatic", "Manual"});
        gbc.gridx = 2; gbc.gridy = row; gbc.gridwidth = 2;
        params.add(positioningCombo, gbc);
        row++;

        btnAddRoi = new JButton("Add Current ROI");
        btnClearRois = new JButton("Clear all");
        lblRoiCount = new JLabel("ROIs: 0");
        
        btnAddRoi.setEnabled(false);
        btnClearRois.setEnabled(false);

        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 2;
        params.add(btnAddRoi, gbc);
        gbc.gridx = 2; gbc.gridy = row; gbc.gridwidth = 1;
        params.add(btnClearRois, gbc);
        gbc.gridx = 3; gbc.gridy = row; gbc.gridwidth = 1;
        params.add(lblRoiCount, gbc);
        row++; 

        gbc.gridx = 0; gbc.gridy = row;
        params.add(new JLabel("ROI size (px):"), gbc);
        recordPixelSizeField = new JTextField("64", 6);
        gbc.gridx = 1; gbc.gridy = row;
        params.add(recordPixelSizeField, gbc);

        gbc.gridx = 2; gbc.gridy = row;
        params.add(new JLabel("N. ROIs:"), gbc);
        recordNumberField = new JTextField("12", 6);
        gbc.gridx = 3; gbc.gridy = row;
        params.add(recordNumberField, gbc);
        row++;
        
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 2;
        params.add(new JLabel("Distribution Radius (px):"), gbc);
        recordRadiusField = new JTextField("64", 6);
        gbc.gridx = 2; gbc.gridy = row; gbc.gridwidth = 2;
        params.add(recordRadiusField, gbc);
        row++;

        detrendCheck = new JCheckBox("Detrend ROIs with 2nd order polynomial fit", true);
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 4;
        params.add(detrendCheck, gbc);
        row++;

        JButton calcBtn = new JButton("Calculate NPS");
        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 4;
        params.add(calcBtn, gbc);
        row++;

        gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 4;
        params.add(Box.createVerticalGlue(), gbc);

        positioningCombo.addActionListener(e -> {
            boolean isManual = positioningCombo.getSelectedItem().equals("Manual");
            toggleManualControls(isManual);
        });

        btnAddRoi.addActionListener(e -> {
            ImagePlus imp = IJ.getImage();
            if (imp == null) {
                IJ.showMessage("No open images.");
                return;
            }
            Roi roi = imp.getRoi();
            if (roi == null) {
                IJ.showMessage("Draw a ROI first");
                return;
            }
            userRois.add((Roi) roi.clone());
            lblRoiCount.setText("Manual ROIs: " + userRois.size());
            
            RoiManager rm = RoiManager.getInstance();
            if (rm == null) rm = new RoiManager();
            rm.addRoi(roi);
        });

        btnClearRois.addActionListener(e -> {
            userRois.clear();
            lblRoiCount.setText("ROIs: 0");
            RoiManager rm = RoiManager.getInstance();
            if (rm != null) rm.reset();
        });

        JPanel results = new JPanel(new GridBagLayout());
        results.setBorder(BorderFactory.createTitledBorder("Results"));

        int rrow = 0;
        gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 6, 4, 6); 
        gbc.fill = GridBagConstraints.HORIZONTAL;

        gbc.gridx = 0; gbc.gridy = rrow;
        results.add(new JLabel("Variance:"), gbc);
        varianceField = new JTextField("-", 5);
        varianceField.setEditable(false);
        gbc.gridx = 1; gbc.gridy = rrow;
        results.add(varianceField, gbc);
        rrow++;

        gbc.gridx = 0; gbc.gridy = rrow;
        results.add(new JLabel("Integrated NPS:"), gbc);
        integratedField = new JTextField("-", 5);
        integratedField.setEditable(false);
        gbc.gridx = 1; gbc.gridy = rrow;
        results.add(integratedField, gbc);
        rrow++;
        
        gbc.gridx = 0; gbc.gridy = rrow;
        results.add(new JLabel("Peak Freq (cy/mm):"), gbc);
        peakFreqField = new JTextField("-", 5);
        peakFreqField.setEditable(false);
        gbc.gridx = 1; gbc.gridy = rrow;
        results.add(peakFreqField, gbc);
        rrow++;
        
        gbc.gridx = 0; gbc.gridy = rrow;
        results.add(new JLabel("Mean Freq (cy/mm):"), gbc);
        meanFreqField = new JTextField("-", 5);
        meanFreqField.setEditable(false);
        gbc.gridx = 1; gbc.gridy = rrow;
        results.add(meanFreqField, gbc);
        rrow++;

        JButton show2DBtn = new JButton("Show 2D NPS");
        JButton showDataBtn = new JButton("Show 1D NPS Table");
        
        gbc.gridx = 0; gbc.gridy = rrow; gbc.gridwidth = 1;
        results.add(show2DBtn, gbc);
        
        gbc.gridx = 1; gbc.gridy = rrow; gbc.gridwidth = 1;
        results.add(showDataBtn, gbc);
        rrow++;

        gbc.gridx = 0; gbc.gridy = rrow; 
        gbc.gridwidth = 2; 
        gbc.weightx = 1.0; 
        gbc.weighty = 1.0; 
        gbc.fill = GridBagConstraints.BOTH;
        results.add(NpsGraphManager.createChartPanel(), gbc);

        GridBagConstraints g = new GridBagConstraints();
        g.insets = new Insets(8, 8, 8, 8);
        g.fill = GridBagConstraints.BOTH;

        g.gridx = 0; g.gridy = 0; 
        g.weightx = 0.3; 
        g.weighty = 1.0;
        root.add(params, g);

        g.gridx = 1; g.gridy = 0; 
        g.weightx = 0.7;
        g.weighty = 1.0;
        root.add(results, g);

        calcBtn.addActionListener(e -> doCalculate());
        show2DBtn.addActionListener(e -> {
            if (lastNPS2D != null) {
                NPScalc.showPowerSpectrum(lastNPS2D, "NPS");
            } else {
                JOptionPane.showMessageDialog(root, "Run calculation first.");
            }
        });
        showDataBtn.addActionListener(e -> {
            if (lastRadialFreq == null || lastRadialNPS == null) {
                JOptionPane.showMessageDialog(null, "Run NPS calculation first.", "No Data", JOptionPane.WARNING_MESSAGE);
                return;
            }
            String[] columns = {"Freq (cycles/mm)", "NPS (HU²·mm²)"};
            List<double[]> data = Arrays.asList(lastRadialFreq, lastRadialNPS);
            DataTableFrame tableFrame = new DataTableFrame("1D Radial NPS Data", columns, data);
            tableFrame.setVisible(true);
        });        
        
        return root;
    }

    private static void toggleManualControls(boolean isManual) {
        btnAddRoi.setEnabled(isManual);
        btnClearRois.setEnabled(isManual);
        recordNumberField.setEnabled(!isManual);
        recordPixelSizeField.setEnabled(!isManual);
        recordRadiusField.setEnabled(!isManual);
        
        startSliceField.setEnabled(true);
        endSliceField.setEnabled(true);
    }

    private static void doCalculate() {
        try {
            SeriesData sd = ImageOpener.getActiveSeriesData();
            if (sd == null) {
                JOptionPane.showMessageDialog(null, "No active series. Use the file browser to open images.");
                return;
            }

            int total = sd.getFiles().size();
            int start = Integer.parseInt(startSliceField.getText().trim());
            int end   = Integer.parseInt(endSliceField.getText().trim());
            
            if (start < 1 || end > total || start > end) {
                JOptionPane.showMessageDialog(null, "Invalid slice range.");
                return;
            }

            boolean isAutoMode = positioningCombo.getSelectedItem().equals("Automatic");
            int padSize = 0;

            if (isAutoMode) {
                padSize = Integer.parseInt(recordPixelSizeField.getText().trim());
            } else {
                if (userRois.isEmpty()) {
                    JOptionPane.showMessageDialog(null, "No manual ROIs defined. Add ROIs first.");
                    return;
                }
                for (Roi userRoi : userRois) {
                    Rectangle bounds = userRoi.getBounds();
                    int maxDim = Math.max(bounds.width, bounds.height);
                    if (maxDim > padSize) {
                        padSize = maxDim;
                    }
                }
            }
            
            List<File> files = sd.getFiles();
            ImageStack stack = new ImageStack();
            Calibration cal = null;
            
            for (int i = start - 1; i < end; i++) {
                File file = files.get(i);
                ImagePlus tempImp = new ImagePlus(file.getAbsolutePath());
                if (tempImp.getWidth() == 0) continue;
                if (cal == null) cal = tempImp.getCalibration();
                stack.addSlice(tempImp.getProcessor());
            }

            if (stack.getSize() == 0) {
                JOptionPane.showMessageDialog(null, "Could not load images.");
                return;
            }
            
            double[] varianceAccumulator = {0.0}; 
            
            double[][] sumNPS2D = new double[padSize][padSize];
            int validRecordsTotal = 0;

            if (isAutoMode) {
                int recordPx = padSize;
                int recordNum = Integer.parseInt(recordNumberField.getText().trim());
                double recordRadiusMM = Double.parseDouble(recordRadiusField.getText().trim()); // Parse do novo raio
                boolean addToManager = true;

                RoiManager manager = RoiManager.getInstance();
                if (addToManager) {
                    if (manager == null) manager = new RoiManager();
                    manager.reset();
                }

                for (int z = start; z <= end; z++) {
                    int stackIndex = z - start + 1;
                    ImageProcessor ip = stack.getProcessor(stackIndex);
                    PolygonRoi phantom = Segmenter.segmenter(ip, cal, -200, 5000, 0.2, false, 0, 20);
                    
                    
                    List<Roi> records = Segmenter.createCircularRecordRois(
                            ip, cal, phantom, recordNum, recordPx, recordRadiusMM, addToManager, z);

                    for (Roi r : records) {
                        validRecordsTotal = processRoiForNPS(ip, cal, r, padSize, varianceAccumulator, sumNPS2D, validRecordsTotal);
                    }
                }
                if (addToManager && manager != null) manager.runCommand("Show All");

            } else {
                for (int z = start; z <= end; z++) {
                    int stackIndex = z - start + 1;
                    ImageProcessor ip = stack.getProcessor(stackIndex);
                    for (Roi userRoi : userRois) {
                        validRecordsTotal = processRoiForNPS(ip, cal, userRoi, padSize, varianceAccumulator, sumNPS2D, validRecordsTotal);
                    }
                }
            }

            if (validRecordsTotal == 0) {
                JOptionPane.showMessageDialog(null, "No valid record ROIs were generated or processed.");
                return;
            }

            for (int y = 0; y < padSize; y++) {
                for (int x = 0; x < padSize; x++) {
                    sumNPS2D[y][x] /= validRecordsTotal;
                }
            }

            double[] radialFreq = new double[padSize / 2];
            double[] radialNPS = NPScalc.radialProfile(sumNPS2D, cal, radialFreq);
            double df_spatial = 1.0 / (padSize * cal.pixelWidth);
            double freqPico = encontrarFrequenciaDePico(radialFreq, radialNPS, 4);
            
            int indicePico = (int) Math.round(freqPico / df_spatial);
            int numPontosParaAjuste = Math.max(1, indicePico); 

            WeightedObservedPoints obs = new WeightedObservedPoints();
            for (int i = 1; i <= numPontosParaAjuste; i++) {
                if(i < radialFreq.length) { 
                    obs.add(radialFreq[i], radialNPS[i]);
                }
            }

            PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2);
            double[] coeff = fitter.fit(obs.toList());
            double centerNPS = Math.max(0, coeff[0]);

            radialNPS[0] = centerNPS;
            sumNPS2D[padSize/2][padSize/2] = centerNPS;

            double integratedNPS = 0.0;
            double sumFreqTimesNPS = 0.0;
            double sumNPS_1D = 0.0;
            
            for (int y = 0; y < padSize; y++) {
                for (int x = 0; x < padSize; x++) {
                    integratedNPS += sumNPS2D[y][x];
                }
            }
            integratedNPS *= (1.0 / (padSize * cal.pixelWidth)) * (1.0 / (padSize * cal.pixelHeight));

            for (int i = 0; i < radialFreq.length; i++) {
                sumFreqTimesNPS += (radialFreq[i] * radialNPS[i]);
                sumNPS_1D += radialNPS[i];
            }
            double meanFreq = (sumNPS_1D > 0) ? (sumFreqTimesNPS / sumNPS_1D) : 0.0;

            double varianceMean = varianceAccumulator[0] / validRecordsTotal;

            lastNPS2D = sumNPS2D;
            lastRadialFreq = radialFreq;
            lastRadialNPS = radialNPS;
           
            NpsGraphManager.updateChart(radialFreq, radialNPS);

            varianceField.setText(DF.format(varianceMean));
            integratedField.setText(DF.format(integratedNPS));
            peakFreqField.setText(DF.format(freqPico));
            meanFreqField.setText(DF.format(meanFreq));

        } catch (NumberFormatException ex) {
            IJ.error("Invalid parameters.");
        } catch (Exception ex) {
            IJ.handleException(ex);
        }
    }

    private static int processRoiForNPS(ImageProcessor ip, Calibration cal, Roi r, int padSize, double[] varAccum, double[][] sumNPS2D, int validRecordsTotal) {
        ip.setRoi(r);
        ImageStatistics stats = ImageStatistics.getStatistics(ip);
        
        ImageProcessor record = ip.crop();
        
        int recordW = record.getWidth();
        int recordH = record.getHeight();
        int N = recordW * recordH;
        
        varAccum[0] += Math.pow(stats.stdDev, 2) * (double)(N - 1) / N;
        
        if (detrendCheck.isSelected()) {
            float[][] roiData = record.getFloatArray(); 
            double[][] detrendedRoiData = Detrending.detrendWith2ndOrderPolynomial(roiData);
            
            float[][] detrendedFloatData = new float[recordW][recordH]; 
            
            for (int x = 0; x < recordW; x++) {
                for (int y = 0; y < recordH; y++) {
                    detrendedFloatData[x][y] = (float) detrendedRoiData[x][y];
                }
            }
            record = new ij.process.FloatProcessor(detrendedFloatData);
        }

        double[][] PowerSpectrum = NPScalc.calculate2DNPS(record, cal, padSize);
        ip.resetRoi();

        for (int y = 0; y < padSize; y++) {
            for (int x = 0; x < padSize; x++) {
                sumNPS2D[y][x] += PowerSpectrum[y][x];
            }
        }
        return validRecordsTotal + 1;
    }

    public static double encontrarFrequenciaDePico(double[] freqRadial, double[] npsRadial1D, int larguraJanela) {
        int indicePicoAproximado = -1;
        double maxNPS = 0;
        for (int i = 1; i < npsRadial1D.length; i++) {
            if (npsRadial1D[i] > maxNPS) {
                maxNPS = npsRadial1D[i];
                indicePicoAproximado = i;
            }
        }

        if (indicePicoAproximado == -1) {
            return -1;
        }

        int inicioJanela = Math.max(1, indicePicoAproximado - larguraJanela);
        int fimJanela = Math.min(npsRadial1D.length - 1, indicePicoAproximado + larguraJanela);

        if ((fimJanela - inicioJanela + 1) < 3) {
            return freqRadial[indicePicoAproximado];
        }
        
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for (int i = inicioJanela; i <= fimJanela; i++) {
            obs.add(freqRadial[i], npsRadial1D[i]);
        }

        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2); 
        double[] coeff = fitter.fit(obs.toList()); 

        double a = coeff[2];
        double b = coeff[1];
        
        double frequenciaDePico = -b / (2 * a);
        
        if (a >= 0 || frequenciaDePico < freqRadial[inicioJanela] || frequenciaDePico > freqRadial[fimJanela]) {
            return freqRadial[indicePicoAproximado];
        }
        
        return frequenciaDePico;
    }

    public static double[] getLastRadialFreq() {
        return lastRadialFreq;
    }
    
    public static double[][] getLastNPS2D() {
        return lastNPS2D;
    }

    public static double[] getLastRadialNPS() {
        return lastRadialNPS;
    }
}
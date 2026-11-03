package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.FloatProcessor;

import java.util.ArrayList;
import java.util.List;

public class NpsValidation {

    public static void main(String[] args) {
    	System.out.println("Kernel\tNumROIs\tMean_PENps(%)\tSD_PENps(%)\tFailures");

        int matriz = 64;
        double pixelSize = 250.0 / 512.0;
        double fov = matriz * pixelSize;
        
        double sdConstante = 50.0;
        double variance = sdConstante * sdConstante;
        
        Calibration cal = new Calibration();
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.setUnit("mm");
        
        int[] numRoisList = {1, 5, 10, 15, 20, 40, 60};
        String[] kernels = {"SOFT", "LUNG"};
        int numMeasurements = 30;

        for (String kernel : kernels) {
            
            double[] npsParams = NPSModel.getNPSCoefficients(kernel);
            NPSModel npsModel = new NPSModel();
            
            double df = 1.0 / (matriz * pixelSize);
            double nyquist = 1.0 / (2.0 * pixelSize);
            
            int numBins = matriz / 2;
            double[] npsModeloRadial = new double[numBins];
            for (int i = 0; i < numBins; i++) {
                double freq = i * df;
                if (freq <= nyquist) {
                    npsModeloRadial[i] = npsModel.getNPSValue(freq, npsParams, variance);
                }
            }

            for (int numRois : numRoisList) {
                List<Double> peNpsList = new ArrayList<>();
                int failures = 0;

                for (int i = 0; i < numMeasurements; i++) {
                    try {
                        PhantomImageCreator creator = new PhantomImageCreator();
                        ImagePlus imp = creator.criarImagemCTPhantom(
                                numRois,            // nSlicesHomogeneo = N ROIs
                                0,                  // nSlicesInserts = 0 (Apenas ruído homogêneo)
                                matriz,             // matriz = 64
                                fov,                // fov = 31.25
                                0.0,                // bgHU = 0
                                0.0,                // objectMeanHU = 0
                                500,                // diametroMM = 0
                                0,                  // trendAmpHU = 0
                                new double[0],      // insertHU vazio
                                0.0,                // insertRaioMM = 0
                                0.0,                // insertDistanciaMM = 0
                                kernel,             // Kernel alvo
                                variance            // Ruído alvo
                        );
                        
                        ImageStack stack = imp.getStack();
                        List<double[][]> all2DNps = new ArrayList<>();
                        
                        for (int s = 1; s <= numRois; s++) {
                            FloatProcessor record = (FloatProcessor) stack.getProcessor(s).convertToFloat();
                            float[][] roiData = record.getFloatArray();
                            
                            double[][] detrendedRoiData = Detrending.detrendWith2ndOrderPolynomial(roiData);
                            
                            float[][] detrendedFloatData = new float[record.getHeight()][record.getWidth()];
                            for(int i1 = 0; i1 < record.getHeight(); i1++) {
                                for(int j1 = 0; j1 < record.getWidth(); j1++) {
                                    detrendedFloatData[i1][j1] = (float) detrendedRoiData[i1][j1];
                                }
                            }
                            record = new FloatProcessor(detrendedFloatData);

                            double[][] powerSpectrum = NPScalc.calculate2DNPS(record, cal, matriz);
                            all2DNps.add(powerSpectrum);
                        }
                        
                        double[][] avg2DNps = media2D(all2DNps, matriz);
                        double[] radialFreq = new double[matriz / 2];

                        double[] npsMedidaRadial = NPScalc.radialProfile(avg2DNps, cal, radialFreq);
                        
                        double pe = calcularPENps(npsMedidaRadial, npsModeloRadial, numBins);
                        if (!Double.isNaN(pe)) {
                            peNpsList.add(pe);
                        } else {
                            failures++;
                        }
                        
                    } catch (Exception e) {
                        failures++;
                    }
                }

                double meanPE = 0;
                for (double p : peNpsList) meanPE += p;
                if (!peNpsList.isEmpty()) meanPE /= peNpsList.size();
                
                double sumSq = 0;
                for (double p : peNpsList) sumSq += Math.pow(p - meanPE, 2);
                double stdPE = peNpsList.isEmpty() ? 0 : Math.sqrt(sumSq / peNpsList.size());
                
                System.out.printf("%s\t%d\t%.4f\t%.4f\t%d\n", kernel, numRois, meanPE, stdPE, failures);
            }
        }
        System.out.println("-------------------------------------------------------------------------");
        System.out.println("Validação Concluída.");
    }

    /**
     * Calcula a média dos espectros 2D.
     */
    private static double[][] media2D(List<double[][]> list2D, int size) {
        double[][] avg = new double[size][size];
        int n = list2D.size();
        for (double[][] ps : list2D) {
            for (int y = 0; y < size; y++) {
                for (int x = 0; x < size; x++) {
                    avg[y][x] += ps[y][x];
                }
            }
        }
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                avg[y][x] /= n;
            }
        }
        return avg;
    }

    /**
     * Calcula o Erro Percentual Integral (PENps)
     */
    private static double calcularPENps(double[] measured, double[] model, int numBins) {
        double sumDiff = 0.0;
        double sumModel = 0.0;
        
        int limit = Math.min(Math.min(measured.length, model.length), numBins);
        
        for (int i = 0; i < limit; i++) {
            sumDiff += Math.abs(measured[i] - model[i]);
            sumModel += model[i];
        }
        
        if (sumModel == 0) return Double.NaN;
        return (sumDiff / sumModel) * 100.0;
    }
}
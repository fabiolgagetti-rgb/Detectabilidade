package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;

import java.util.ArrayList;
import java.util.List;

import br.unifesp.fisicamedica.Segmenter.DetectorResult;

public class TtfValidation {

    public static void main(String[] args) {
    	
    	System.out.println("Kernel\tNoise(SD)\tMean_PEttf(%)\tSD_PEttf(%)\tFalhas");

        int nSlicesHomogeneo = 0;
        int nSlicesInserts = 1;
        int matriz = 512;
        double fov = 250.0;
        double bgHU = -1000.0;
        double objectMeanHU = 0.0; 
        double diametroMM = 200.0; 
        double trendAmpHU = 0.0;
        double[] insertHU = {100};
        double insertRaioMM = 12.5; 
        double insertDistanciaMM = 0.0; 
        
        double pixelSize = fov / matriz; 
        double nyquist = 1.0 / (2.0 * pixelSize); 
        
        Calibration cal = new Calibration();
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.setUnit("mm");
        
        int roiSize = (int ) (insertRaioMM*3/pixelSize);
        int centro = matriz / 2;
        int roiX = centro - (roiSize / 2);
        int roiY = centro - (roiSize / 2);
        Roi userRoi = new Roi(roiX, roiY, roiSize, roiSize);
        
        String[] kernels = {"SOFT", "Lung"}; 
        int[] noiseLevelsSD = {0, 1, 10, 20, 30, 45};
        int numMeasurements = 30;
        int fftPoints = 64; 
        
        double fNy = 1/(2*pixelSize);
        
        double df = 2*fNy/fftPoints;

        for (String kernel : kernels) {

            double[] modelMtf = new double[fftPoints];
            try {
                double[][] gaussians = MTFModel.getMTFparameters(kernel);
                
                for (int i = 0; i < fftPoints; i++) {
                    double freq = i * df;
                    modelMtf[i] = MTFModel.getMTFValue(freq, gaussians);
                }
                
                double dcValue = modelMtf[0];
                if (dcValue > 0) {
                    for (int i = 0; i < fftPoints; i++) {
                        modelMtf[i] /= dcValue;
                    }
                }
                
            } catch (Exception e) {
                System.out.println("ERRO " + kernel);
                continue; 
            }

            for (int sd : noiseLevelsSD) {
                double variance = sd * sd; 
                
                List<Double> peList = new ArrayList<>();
                int failures = 0;
                
                for (int i = 0; i < numMeasurements; i++) {
                    try {
                        PhantomImageCreator creator = new PhantomImageCreator();
                        ImagePlus imp = creator.criarImagemCTPhantom(
                            nSlicesHomogeneo, nSlicesInserts, matriz, fov, bgHU, objectMeanHU, 
                            diametroMM, trendAmpHU, insertHU, insertRaioMM, insertDistanciaMM, 
                            kernel, variance
                        );
                        
                        DetectorResult dR = Segmenter.edgeFinder(imp.getStack(), cal, userRoi, false, true);
                        
                        if (dR == null) {
                            failures++;i--;
                            continue;
                        }
                        
                        TtfCalc.mtfResult mtfRes = TtfCalc.analisarStackComInserts(imp.getStack(), cal, dR, fftPoints);
                        
                        if (mtfRes == null || mtfRes.ttfs == null || mtfRes.ttfs.isEmpty()) {
                            failures++;
                            continue;
                        }

                        double[] measuredMtf = mtfRes.ttfs.values().iterator().next();
                        
                        double pe = calcularPEttf(measuredMtf, modelMtf, df, nyquist);
                        if (!Double.isNaN(pe)) {
                            peList.add(pe);
                        } else {
                            failures++;
                        }
                        
                    } catch (Exception e) {
                        failures++;
                    }
                }
                
                double meanPE = 0;
                for (double p : peList) {
                    meanPE += p;
                }
                if (!peList.isEmpty()) {
                    meanPE /= peList.size();
                }
                
                double sumSq = 0;
                for (double p : peList) {
                    sumSq += Math.pow(p - meanPE, 2);
                }
                double stdPE = peList.isEmpty() ? 0 : Math.sqrt(sumSq / peList.size());
                
                System.out.printf("%s\t%d\t%.4f\t%.4f\t%d\n", kernel, sd, meanPE, stdPE, failures);
            }
        }
    }

    /**
     * Calcula o Erro Percentual (PEttf) comparando a MTF medida e a MTF modelo (Bujila).
     */
    private static double calcularPEttf(double[] measuredMtf, double[] model, double df, double nyquist) {
        double sumDiff = 0.0;
        double sumModel = 0.0;
        
        int limit = Math.min(measuredMtf.length, model.length);
        
        for (int i = 0; i < limit; i++) {
            double freq = i * df;
            
            if (freq > nyquist) {
                break;
            }
            
            sumDiff += Math.abs(measuredMtf[i] - model[i]);
            sumModel += model[i];
        }
        
        if (sumModel == 0) return Double.NaN; 
        
        return (sumDiff / sumModel) * 100.0;
    }
}
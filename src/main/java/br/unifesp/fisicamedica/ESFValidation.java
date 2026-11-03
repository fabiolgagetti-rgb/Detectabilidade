package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.List;

import br.unifesp.fisicamedica.Segmenter.DetectorResult;
import br.unifesp.fisicamedica.TtfCalc.mtfResult;

public class ESFValidation {

    public static void main(String[] args) {
        
        System.out.println("Kernel\tNoise(SD)\tRMSE_ESF(HU)\tBias_ESF(HU)\tFalhas");

        int nSlicesHomogeneo = 0;
        int nSlicesInserts = 1;
        int matriz = 512;
        double fov = 250.0;
        double bgHU = -1000.0;
        double objectMeanHU = 0.0; 
        double diametroMM = 200.0; 
        double trendAmpHU = 0.0;
        double[] insertHU = {100}; // Contraste do inserto
        double insertRaioMM = 12.5; 
        double insertDistanciaMM = 0.0; 
        
        double pixelSize = fov / matriz; 
        
        Calibration cal = new Calibration();
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.setUnit("mm");
        
        int roiSize = (int) (insertRaioMM * 3 / pixelSize);
        double gtCentroPx = matriz / 2.0; 
        int roiX = (int) gtCentroPx - (roiSize / 2);
        int roiY = (int) gtCentroPx - (roiSize / 2);
        Roi userRoi = new Roi(roiX, roiY, roiSize, roiSize);
        
        String[] kernels = {"SOFT", "Lung"}; 
        int[] noiseLevelsSD = {1, 10, 30, 45}; 
        int numMeasurements = 30;
        int fftPoints = 64; 
        
        // Tamanho exato do array da ESF calculado no seu TtfCalc
        int profileLength = (int) (10 / cal.pixelWidth);
        profileLength *= 2;
        if (profileLength % 2 == 0) profileLength++;

        for (String kernel : kernels) {
            
            // 1. GERAR A ESF GROUND TRUTH ABSOLUTA (Direto da Matriz 1D)
            double[] esfIdeal = null;
            try {
                PhantomImageCreator creator = new PhantomImageCreator();
                ImagePlus impZero = creator.criarImagemCTPhantom(
                    nSlicesHomogeneo, nSlicesInserts, matriz, fov, bgHU, objectMeanHU, 
                    diametroMM, trendAmpHU, insertHU, insertRaioMM, insertDistanciaMM, 
                    kernel, 0.0 // Ruído ZERO
                );
                
                ImageProcessor ipZero = impZero.getStack().getProcessor(1);
                
                // Extrai a ESF exatamente às 12 horas
                esfIdeal = gerarEsfIdealDireta(ipZero, gtCentroPx, insertRaioMM / pixelSize, profileLength, cal);
                
            } catch (Exception e) {
                System.out.println("Erro ao gerar Ground Truth Absoluto para " + kernel);
                continue;
            }

            // 2. TESTAR COM RUÍDO E O PIPELINE COMPLETO
            for (int sd : noiseLevelsSD) {
                double variance = sd * sd; 
                
                List<Double> rmseList = new ArrayList<>();
                List<Double> biasList = new ArrayList<>();
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
                            failures++; i--; continue;
                        }
                        
                        mtfResult mtfRes = TtfCalc.analisarStackComInserts(imp.getStack(), cal, dR, fftPoints);
                        if (mtfRes == null || mtfRes.esfs == null || mtfRes.esfs.isEmpty()) {
                            failures++; continue;
                        }

                        double[] measuredEsf = mtfRes.esfs.values().iterator().next();
                        
                        // Validar a ESF medida contra a ESF matemática absoluta
                        double[] metrics = calcularErroESF(measuredEsf, esfIdeal);
                        if (!Double.isNaN(metrics[0])) {
                            rmseList.add(metrics[0]); // RMSE
                            biasList.add(metrics[1]); // BIAS
                        } else {
                            failures++;
                        }
                        
                    } catch (Exception e) {
                        failures++;
                    }
                }
                
                double meanRmse = calculateMean(rmseList);
                double meanBias = calculateMean(biasList);
                System.out.printf("%s\t%d\t%.4f\t%.4f\t%d\n", kernel, sd, meanRmse, meanBias, failures);
            }
        }
    }

    /**
     * Extrai a ESF de referência matematicamente perfeita isolando a coluna central 
     * (amostragem bruta) e aplicando interpolação linear estritamente em 1D 
     * para evitar a suavização cruzada da interpolação bilinear 2D.
     */
    private static double[] gerarEsfIdealDireta(ImageProcessor ip, double centroPx, double raioPx, int profileLength, Calibration cal) {
        int half = profileLength / 2;
        double step = 0.5; // Espaçamento subpixel do TtfCalc
        double[] esf = new double[profileLength];

        // Isola a coluna central exata (amostragem bruta em X)
        int x = (int) Math.round(centroPx);
        
        // A borda "12 horas" fica em y = centro - raio (lembrando que Y cresce para baixo)
        double edgeY = centroPx - raioPx;

        for (int j = -half; j <= half; j++) {
            // j negativo -> dentro do inserto (y maior, voltando pro centro)
            // j positivo -> fora do inserto (y menor, subindo na tela)
            double sampleY = edgeY + (j * step);
            
            // 1. Amostragem bruta na coluna isolada (Sem interferência 2D)
            int y0 = (int) Math.floor(sampleY);
            int y1 = y0 + 1;
            
            // Limites de segurança da matriz
            if (y0 < 0) y0 = 0; 
            if (y0 >= ip.getHeight()) y0 = ip.getHeight() - 1;
            if (y1 < 0) y1 = 0; 
            if (y1 >= ip.getHeight()) y1 = ip.getHeight() - 1;
            
            double val0 = ip.getf(x, y0);
            double val1 = ip.getf(x, y1);
            
            // 2. Interpolação linear puramente 1D no eixo Y
            double dy = sampleY - Math.floor(sampleY); // Parte fracionária
            double pixelValue = val0 * (1.0 - dy) + val1 * dy;
            
            esf[j + half] = cal.getCValue(pixelValue);
        }
        
        return esf;
    }

    private static double[] calcularErroESF(double[] measured, double[] ideal) {
        if (measured.length != ideal.length) return new double[]{Double.NaN, Double.NaN};
        
        double sumSqError = 0.0;
        double sumBias = 0.0;
        int n = measured.length;
        
        for (int i = 0; i < n; i++) {
            double error = measured[i] - ideal[i];
            sumBias += error;
            sumSqError += (error * error);
        }
        
        double rmse = Math.sqrt(sumSqError / n);
        double bias = sumBias / n;
        
        return new double[]{rmse, bias};
    }

    private static double calculateMean(List<Double> values) {
        if (values == null || values.isEmpty()) return Double.NaN;
        double sum = 0;
        for (double v : values) sum += v;
        return sum / values.size();
    }
}
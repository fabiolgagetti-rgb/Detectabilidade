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
     * ESF de referência (ground truth).
     *
     * IMPORTANTE: para que o erro medido (measured - ideal) reflita apenas o
     * efeito de ruido + detector, a ESF ideal precisa ser amostrada com a MESMA
     * geometria e o MESMO kernel de interpolacao que a ESF medida em
     * TtfCalc.calculateESF. Antes isso nao acontecia:
     *   - a medida anda ao longo da NORMAL local e usa interpolacao Lanczos-3;
     *   - a ideal andava so no eixo Y e usava interpolacao linear 1D.
     * Essa diferenca de metodo (Lanczos vs linear; normal vs eixo Y) contaminava
     * o RMSE/Bias com "mismatch de metodo", nao com erro de deteccao.
     *
     * Aqui a borda esta as "12 horas" de um circulo centrado em (centroPx,
     * centroPx), entao a normal e exatamente vertical. Amostramos ao longo de Y
     * com passo 0.5 px, ancorados na borda analitica edgeY = centroPx - raioPx,
     * usando o MESMO Lanczos-3 de TtfCalc (getLanczosInterpolatedPixel). A
     * imagem passada deve ser a versao SEM ruido (MTF aplicada) -- essa e a ESF
     * verdadeira borrada pelo kernel.
     *
     * Convencao de indice identica a calculateESF: indice 0 = fundo (fora do
     * inserto, Y menor), ultimo indice = objeto (dentro, Y maior).
     */
    private static double[] gerarEsfIdealDireta(ImageProcessor ip, double centroPx, double raioPx, int profileLength, Calibration cal) {
        int half = profileLength / 2;
        double step = 0.5; // mesmo passo subpixel do TtfCalc
        double[] esf = new double[profileLength];

        // Coluna central (12 horas). A normal aponta de fora (Y menor) para
        // dentro (Y maior), i.e. (nx, ny) = (0, +1).
        double x0 = centroPx;
        double edgeY = centroPx - raioPx;   // borda analitica

        for (int j = -half; j <= half; j++) {
            double x = x0;
            double y = edgeY + (j * step);  // normal vertical, mesmo sentido da medida

            // Mesmo kernel de interpolacao da ESF medida (Lanczos-3).
            double pixelValue = TtfCalc.getLanczosInterpolatedPixel(ip, x, y);
            esf[j + half] = cal.getCValue(pixelValue);
        }

        return esf;
    }

    /**
     * NOTA (registro/registration): a ESF medida e ancorada nos pontos de borda
     * DETECTADOS, enquanto a ideal e ancorada na borda ANALITICA
     * (centroPx - raioPx). Um deslocamento subpixel entre as duas ancoras
     * aparece como BIAS em HU, amplificado pela inclinacao da ESF. Portanto o
     * campo Bias_ESF mede, em parte, erro de REGISTRO e nao so erro de
     * intensidade. Para isolar a fidelidade de forma da ESF, considere alinhar
     * measured e ideal pelo ponto de 50% (ou pelo pico da correlacao cruzada)
     * antes de calcular o RMSE, e reportar o offset de alinhamento separadamente.
     */
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
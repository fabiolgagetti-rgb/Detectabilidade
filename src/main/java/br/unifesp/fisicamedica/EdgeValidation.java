package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.process.FloatPolygon;
import java.awt.Point;
import java.util.ArrayList;
import java.util.List;
import br.unifesp.fisicamedica.Segmenter.DetectorResult;

public class EdgeValidation {

    public static void main(String[] args) {
        
        System.out.println("Kernel\tNoise(SD)\tCenterShift(mm)\tRadialBias(mm)\tRadialRMSE(mm)\tFalhas");

        int nSlicesHomogeneo = 0;
        int nSlicesInserts = 1;
        int matriz = 512;
        double fov = 250.0;
        double bgHU = -1000.0;
        double objectMeanHU = 0.0; 
        double diametroMM = 200.0; 
        double trendAmpHU = 0.0;
        double[] insertHU = {100}; // Contraste do inserto
        double insertRaioMM = 12.5; // Raio Verdadeiro (Ground Truth)
        double insertDistanciaMM = 0.0; // Fica exatamente no centro
        
        double pixelSize = fov / matriz; 
        
        Calibration cal = new Calibration();
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.setUnit("mm");
        
        int roiSize = (int) (insertRaioMM * 3 / pixelSize);
        double gtCentroPx = matriz / 2.0; // Centro Verdadeiro (Ground Truth) em pixels
        
        int roiX = (int) gtCentroPx - (roiSize / 2);
        int roiY = (int) gtCentroPx - (roiSize / 2);
        Roi userRoi = new Roi(roiX, roiY, roiSize, roiSize);
        
        String[] kernels = {"SOFT", "Lung"}; 
        int[] noiseLevelsSD = {0, 1, 10, 20, 30, 40};
        int numMeasurements = 30;

        for (String kernel : kernels) {
            for (int sd : noiseLevelsSD) {
                double variance = sd * sd; 
                
                List<Double> centerShiftList = new ArrayList<>();
                List<Double> radialBiasList = new ArrayList<>();
                List<Double> radialRmseList = new ArrayList<>();
                int failures = 0;
                
                for (int i = 0; i < numMeasurements; i++) {
                    try {
                        PhantomImageCreator creator = new PhantomImageCreator();
                        ImagePlus imp = creator.criarImagemCTPhantom(
                            nSlicesHomogeneo, nSlicesInserts, matriz, fov, bgHU, objectMeanHU, 
                            diametroMM, trendAmpHU, insertHU, insertRaioMM, insertDistanciaMM, 
                            kernel, variance
                        );
                        
                        // Passamos false, true (assumindo que o último é para detectar insertos)
                        DetectorResult dR = Segmenter.edgeFinder(imp.getStack(), cal, userRoi, false, true);
                        
                        if (dR == null || dR.rois.isEmpty() || dR.rois.get(0).isEmpty()) {
                            failures++;
                            continue;
                        }
                        
                        // Pega o centro detectado e o polígono da borda subpixel
                        Point detectedCenter = dR.Pos.get(0);
                        Roi borderRoi = dR.rois.get(0).get(0);
                        FloatPolygon boundary = borderRoi.getFloatPolygon();
                        
                        if (boundary == null || boundary.npoints == 0) {
                            failures++;
                            continue;
                        }

                        // 1. Calcula o Erro do Centro em mm
                        double dx = (detectedCenter.x - gtCentroPx) * pixelSize;
                        double dy = (detectedCenter.y - gtCentroPx) * pixelSize;
                        double centerShiftMM = Math.sqrt(dx * dx + dy * dy);
                        
                        // 2 e 3. Calcula o Viés e o RMSE da borda detectada em mm
                        double[] radialMetrics = calcularMetricasRadiais(boundary, gtCentroPx, gtCentroPx, insertRaioMM, pixelSize);
                        
                        centerShiftList.add(centerShiftMM);
                        radialBiasList.add(radialMetrics[0]);
                        radialRmseList.add(radialMetrics[1]);
                        
                    } catch (Exception e) {
                        failures++;
                    }
                }
                
                // Agregação dos resultados (Média de cada métrica)
                double meanCenterShift = calculateMean(centerShiftList);
                double meanRadialBias = calculateMean(radialBiasList);
                double meanRadialRmse = calculateMean(radialRmseList);
                
                // Print formatado para ir direto para o Excel/Origin
                System.out.printf("%s\t%d\t%.4f\t%.4f\t%.4f\t%d\n", 
                                  kernel, sd, meanCenterShift, meanRadialBias, meanRadialRmse, failures);
            }
        }
    }

    /**
     * Calcula o viés (bias) e o erro quadrático médio (RMSE) da borda detectada
     * comparada com o círculo ideal.
     * Retorna: [0] = Bias (mm), [1] = RMSE (mm)
     */
    private static double[] calcularMetricasRadiais(FloatPolygon boundary, double gtCenterX, double gtCenterY, double gtRadiusMM, double pixelSize) {
        double sumError = 0.0;
        double sumSqError = 0.0;
        int n = boundary.npoints;
        
        for (int k = 0; k < n; k++) {
            // Pega a coordenada subpixel do ponto da borda
            float px = boundary.xpoints[k];
            float py = boundary.ypoints[k];
            
            // Distância do ponto até o centro verdadeiro (em pixels)
            double distPx = Math.sqrt(Math.pow(px - gtCenterX, 2) + Math.pow(py - gtCenterY, 2));
            
            // Converte para MM e acha o erro em relação ao raio perfeito
            double distMM = distPx * pixelSize;
            double errorMM = distMM - gtRadiusMM; 
            
            sumError += errorMM;
            sumSqError += (errorMM * errorMM);
        }
        
        double bias = sumError / n; // Se positivo, a borda foi detectada maior que o real
        double rmse = Math.sqrt(sumSqError / n); // Espalhamento geral da borda (rugosidade)
        
        return new double[]{bias, rmse};
    }

    private static double calculateMean(List<Double> values) {
        if (values == null || values.isEmpty()) return Double.NaN;
        double sum = 0;
        for (double v : values) sum += v;
        return sum / values.size();
    }
}
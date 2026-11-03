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
        
        // RMSE_raw/Bias_raw: sem alinhar (contem o erro de registro).
        // RMSE_align: so a fidelidade de forma/intensidade, apos remover o
        //             deslocamento lateral.
        // Offset_px:  vies de posicionamento da borda (detectada vs analitica).
        // FailReasons: quebra das falhas por motivo (ver EdgeDiagnostics).
        EdgeDiagnostics.enable();
        System.out.println("Kernel\tNoise(SD)\tRMSE_raw(HU)\tBias_raw(HU)\tRMSE_align(HU)\tOffset(px)\tFalhas\tFailReasons");

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
                List<Double> rmseAlignList = new ArrayList<>();
                List<Double> offsetList = new ArrayList<>();
                int failures = 0;

                EdgeDiagnostics.reset();  // zera a contagem de motivos para esta condicao

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
                            rmseList.add(metrics[0]);      // RMSE_raw
                            biasList.add(metrics[1]);      // Bias_raw
                            if (!Double.isNaN(metrics[2])) rmseAlignList.add(metrics[2]); // RMSE_align
                            if (!Double.isNaN(metrics[3])) offsetList.add(metrics[3]);    // Offset_px
                        } else {
                            failures++;
                        }
                        
                    } catch (Exception e) {
                        failures++;
                    }
                }
                
                double meanRmse = calculateMean(rmseList);
                double meanBias = calculateMean(biasList);
                double meanRmseAlign = calculateMean(rmseAlignList);
                double meanOffset = calculateMean(offsetList);
                System.out.printf("%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%s\n",
                        kernel, sd, meanRmse, meanBias, meanRmseAlign, meanOffset,
                        failures, EdgeDiagnostics.summary());
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
     * Erros da ESF, separando REGISTRO de FIDELIDADE DE FORMA.
     *
     * Por que alinhar? A ESF ideal e ancorada na borda ANALITICA (passo j=0 cai
     * exatamente na borda verdadeira), entao seu ponto de 50% fica no indice
     * central. Ja a ESF medida e ancorada no ponto de borda DETECTADO
     * (TtfCalc usa x0,y0 do poligono). O detector e subpixel mas tem um
     * deslocamento residual -- em parte aleatorio (RMS ~0.06..0.5 px) e em parte
     * SISTEMATICO (o pre-blur Gaussiano empurra a borda de um circulo levemente
     * para dentro; restauracao/cleanWand/snap de 3px tambem deslocam). Logo o
     * ponto de 50% da ESF medida NAO cai no indice central: ele fica deslocado
     * pelo erro de posicionamento da borda.
     *
     * Num trecho ingreme a ESF muda ~40-50 HU por passo de 0.5 px, entao um
     * deslocamento de 0.2 px vira dezenas de HU de erro aparente, com sinal
     * consistente -> infla RMSE e, principalmente, o Bias, sem que as
     * INTENSIDADES reconstruidas estejam erradas.
     *
     * Solucao: medir o offset entre as duas curvas pelo cruzamento de 50% do
     * contraste, alinhar a medida (reamostrando com o MESMO Lanczos-3) e so
     * entao calcular o RMSE. Reportamos:
     *   [0] rmseRaw   - RMSE sem alinhar (como antes)
     *   [1] biasRaw   - Bias sem alinhar (como antes)
     *   [2] rmseAlign - RMSE depois de remover o deslocamento lateral
     *                   (fidelidade de FORMA/intensidade, ~independente da
     *                    posicao da borda)
     *   [3] offsetPx  - deslocamento de registro em pixels (positivo = borda
     *                   detectada deslocada no sentido +indice em relacao a
     *                   analitica); mede o VIES DE POSICIONAMENTO do detector.
     */
    private static double[] calcularErroESF(double[] measured, double[] ideal) {
        if (measured.length != ideal.length)
            return new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN};

        int n = measured.length;
        double step = 0.5; // px por amostra, igual ao TtfCalc

        // --- metricas cruas (sem alinhar) ---
        double sumSq = 0.0, sumBias = 0.0;
        for (int i = 0; i < n; i++) {
            double e = measured[i] - ideal[i];
            sumBias += e;
            sumSq += e * e;
        }
        double rmseRaw = Math.sqrt(sumSq / n);
        double biasRaw = sumBias / n;

        // --- offset de registro pelo cruzamento de 50% do contraste ---
        double idxMeas = halfRiseIndex(measured);
        double idxIdeal = halfRiseIndex(ideal);
        double rmseAlign = Double.NaN, offsetPx = Double.NaN;
        if (!Double.isNaN(idxMeas) && !Double.isNaN(idxIdeal)) {
            double shift = idxMeas - idxIdeal;          // em amostras (passos)
            offsetPx = shift * step;                    // em pixels

            // Reamostra a medida deslocada de -shift para casar com a ideal,
            // usando Lanczos-3 (mesmo kernel da ESF). Compara apenas a regiao
            // interior valida (margem = raio do Lanczos) para evitar borda.
            int a = 3;
            double sq = 0; int cnt = 0;
            for (int i = 0; i < n; i++) {
                double src = i + shift;                 // posicao na medida que corresponde ao indice i da ideal
                if (src < a || src > n - 1 - a) continue; // fora do suporte seguro do kernel
                double m = lanczosResample(measured, src, a);
                double e = m - ideal[i];
                sq += e * e; cnt++;
            }
            if (cnt > 0) rmseAlign = Math.sqrt(sq / cnt);
        }

        return new double[]{rmseRaw, biasRaw, rmseAlign, offsetPx};
    }

    /**
     * Indice fracionario onde a ESF cruza 50% entre seu nivel de fundo (media
     * dos 3 primeiros pontos) e de objeto (media dos 3 ultimos). Interpolacao
     * linear entre as duas amostras que cercam o cruzamento. Robusto para ESF
     * monotonica (SOFT) e tambem para a nao-monotonica do Lung, pois usa o
     * primeiro cruzamento do nivel medio na subida.
     */
    private static double halfRiseIndex(double[] esf) {
        int n = esf.length;
        if (n < 7) return Double.NaN;
        double bg = (esf[0] + esf[1] + esf[2]) / 3.0;
        double obj = (esf[n - 1] + esf[n - 2] + esf[n - 3]) / 3.0;
        double mid = 0.5 * (bg + obj);
        boolean rising = obj >= bg;
        for (int i = 0; i < n - 1; i++) {
            double a = esf[i], b = esf[i + 1];
            boolean crosses = rising ? (a <= mid && b >= mid) : (a >= mid && b <= mid);
            if (crosses) {
                double denom = (b - a);
                double frac = denom != 0 ? (mid - a) / denom : 0.0;
                return i + frac;
            }
        }
        return Double.NaN;
    }

    /** Reamostragem Lanczos-3 de um sinal 1D numa posicao fracionaria. */
    private static double lanczosResample(double[] s, double pos, int a) {
        int n = s.length;
        int base = (int) Math.floor(pos);
        double sum = 0, w = 0;
        for (int k = base - a + 1; k <= base + a; k++) {
            int idx = Math.max(0, Math.min(n - 1, k));
            double d = pos - k;
            double wt = lanczos1d(d, a);
            sum += s[idx] * wt;
            w += wt;
        }
        return w != 0 ? sum / w : 0;
    }

    private static double lanczos1d(double x, int a) {
        if (x == 0.0) return 1.0;
        if (Math.abs(x) >= a) return 0.0;
        double piX = Math.PI * x;
        return (a * Math.sin(piX) * Math.sin(piX / a)) / (piX * piX);
    }

    private static double calculateMean(List<Double> values) {
        if (values == null || values.isEmpty()) return Double.NaN;
        double sum = 0;
        for (double v : values) sum += v;
        return sum / values.size();
    }
}
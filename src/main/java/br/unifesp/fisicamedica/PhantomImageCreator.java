package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class PhantomImageCreator {

    /**
     * Cria uma imagem de phantom virtual, contendo uma seção homogênea e uma seção com inserts.
     * * @param nSlicesHomogeneo Número de fatias para a seção homogênea
     * @param nSlicesInserts Número de fatias para a seção com inserts
     * @param matriz Tamanho da matriz de pixels
     * @param fov Tamanho do campo de visão em mm
     * @param bgHU Valor HU do background
     * @param objectMeanHU Valor HU do phantom
     * @param diametroMM Diâmetro do phantom em mm
     * @param trendAmpHU Amplitude da não uniformidade do sinal, do centro à borda
     * @param insertHU Array com os valores HU dos inserts
     * @param insertRaioMM Raio de cada insert em mm
     * @param insertDistanciaMM Distância radial do centro do phantom até os inserts em mm
     * @param kernelName Nome do kernel para recuperar os coeficientes de MTF e NPS
     * @param variance Variância do ruído
     * @return O ImagePlus final
     */
    public ImagePlus criarImagemCTPhantom(
            int nSlicesHomogeneo,
            int nSlicesInserts,
            int matriz,
            double fov,
            double bgHU,
            double objectMeanHU,
            double diametroMM,
            double trendAmpHU,
            double[] insertHU,
            double insertRaioMM,
            double insertDistanciaMM,
            String kernelName,
            double variance
    ) {
        ImageStack stack = new ImageStack(matriz, matriz);
        
        double bg = bgHU + 8192;
        double objectMean = objectMeanHU + 8192;
        
        double pixelSize = fov / matriz;
        int centroX = matriz / 2;
        int centroY = matriz / 2;
        int raioPhantomPixels = (int) Math.round((diametroMM / 2.0) / pixelSize);

        Calibration cal = new Calibration();
        cal.setUnit("mm");
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.pixelDepth = 1.0;
        cal.setValueUnit("HU");
        cal.setFunction(Calibration.STRAIGHT_LINE, new double[]{-8192, 1}, "mm");

        MTFModel mtfModel = new MTFModel();
        double[][] mtfParams = MTFModel.getMTFparameters(kernelName);
        double[][] mtf2D = mtfModel.gerarMTF2D(matriz, matriz, cal, mtfParams);

        NPSModel npsModel = new NPSModel();
        double[] npsParams = NPSModel.getNPSCoefficients(kernelName);
        double[][] nps2D = npsModel.gerarNPS2D(matriz, matriz, cal, npsParams, variance);

        FloatProcessor ipIdealHomogeneo = new FloatProcessor(matriz, matriz);
        FloatProcessor ipIdealInserts = new FloatProcessor(matriz, matriz);

        for (int y = 0; y < matriz; y++) {
            for (int x = 0; x < matriz; x++) {
                int dx = x - centroX;
                int dy = y - centroY;
                double rSquared = dx * dx + dy * dy;
                
                if (rSquared <= raioPhantomPixels * raioPhantomPixels) {
                    double rSquaredNorm = rSquared / (double)(raioPhantomPixels * raioPhantomPixels);
                    double trend = trendAmpHU * (rSquaredNorm - 0.5);
                    float val = (float) (objectMean + trend);
                    
                    ipIdealHomogeneo.putPixelValue(x, y, val);
                    ipIdealInserts.putPixelValue(x, y, val);
                } else {
                    ipIdealHomogeneo.putPixelValue(x, y, bg);
                    ipIdealInserts.putPixelValue(x, y, bg);
                }
            }
        }

        if (insertHU != null && insertHU.length > 0 && nSlicesInserts > 0) {
            int nInserts = insertHU.length;
            int raioInsertPixels = (int) Math.round(insertRaioMM / pixelSize);
            int distInsertPixels = (int) Math.round(insertDistanciaMM / pixelSize);

            for (int i = 0; i < nInserts; i++) {
                double angle = 2 * Math.PI * i / nInserts;
                int cx = centroX + (int) Math.round(distInsertPixels * Math.cos(angle));
                int cy = centroY + (int) Math.round(distInsertPixels * Math.sin(angle));
                double valInsert = insertHU[i] + 8192;

                for (int y = cy - raioInsertPixels; y <= cy + raioInsertPixels; y++) {
                    for (int x = cx - raioInsertPixels; x <= cx + raioInsertPixels; x++) {
                        if (x >= 0 && y >= 0 && x < matriz && y < matriz) {
                            int dx = x - cx;
                            int dy = y - cy;
                            if (dx * dx + dy * dy <= raioInsertPixels * raioInsertPixels) {
                                ipIdealInserts.putPixelValue(x, y, valInsert);
                            }
                        }
                    }
                }
            }
        }

        ImageProcessor ipMtfHomogeneo = null;
        if (nSlicesHomogeneo > 0) {
            ipMtfHomogeneo = MTFModel.applyMTF(ipIdealHomogeneo, mtf2D);
        }
        
        ImageProcessor ipMtfInserts = null;
        if (nSlicesInserts > 0) {
            ipMtfInserts = MTFModel.applyMTF(ipIdealInserts, mtf2D);
        }
        
        for (int i = 0; i < nSlicesHomogeneo; i++) {
            ImageProcessor ipRuido = npsModel.texturedNoiseIP(nps2D, pixelSize);
            FloatProcessor sliceFim = adicionarRuido(ipMtfHomogeneo, ipRuido);
            stack.addSlice(criarFatiaShortUnscaled(sliceFim));
        }

        for (int i = 0; i < nSlicesInserts; i++) {
            ImageProcessor ipRuido = npsModel.texturedNoiseIP(nps2D, pixelSize);
            FloatProcessor sliceFim = adicionarRuido(ipMtfInserts, ipRuido);
            stack.addSlice(criarFatiaShortUnscaled(sliceFim));
        }

        ImagePlus imp = new ImagePlus("Phantom_" + kernelName, stack);
        imp.setCalibration(cal);
        return imp;
    }

    /**
     * Creates an anthropomorphic spine/marrow test phantom and runs it through the
     * same MTF + textured-noise pipeline as the circular-insert phantom.
     *
     * Geometry (concentric about the image centre), per the reference object:
     *   - marrow: solid circle, diameter {@code marrowDiamMM} (default 20.5 mm), {@code marrowHU} (200).
     *   - vertebra: ring enclosing the marrow, inner diameter = marrow diameter,
     *     outer diameter {@code ringOuterDiamMM} (24 mm), {@code boneHU} (600).
     *   - wings: a wide isosceles triangle below the marrow, base {@code triBaseMM} (50 mm)
     *     centred horizontally with its base {@code triBaseBelowMM} (17.7 mm) below the
     *     marrow centre, and base angle {@code triBaseAngleDeg} (23 deg) between each
     *     slanted side and the base. Same HU as bone.
     *
     * All structures are rendered with NxN sub-pixel supersampling so the edges
     * (especially the triangle's diagonals) are anti-aliased, giving a clean ESF.
     *
     * @return single-slice ImagePlus (one inserts slice), calibrated in HU.
     */
    public ImagePlus criarImagemCTPhantomSpine(
            int matriz, double fov, double bgHU, double marrowHU, double boneHU,
            double marrowDiamMM, double ringOuterDiamMM,
            double triBaseMM, double triBaseBelowMM, double triBaseAngleDeg,
            String kernelName, double variance) {

        ImageStack stack = new ImageStack(matriz, matriz);
        double bg = bgHU + 8192;
        double marrow = marrowHU + 8192;
        double bone = boneHU + 8192;

        double pixelSize = fov / matriz;
        double cxPx = matriz / 2.0;
        double cyPx = matriz / 2.0;

        Calibration cal = new Calibration();
        cal.setUnit("mm");
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.pixelDepth = 1.0;
        cal.setValueUnit("HU");
        cal.setFunction(Calibration.STRAIGHT_LINE, new double[]{-8192, 1}, "mm");

        MTFModel mtfModel = new MTFModel();
        double[][] mtfParams = MTFModel.getMTFparameters(kernelName);
        double[][] mtf2D = mtfModel.gerarMTF2D(matriz, matriz, cal, mtfParams);
        NPSModel npsModel = new NPSModel();
        double[] npsParams = NPSModel.getNPSCoefficients(kernelName);
        double[][] nps2D = npsModel.gerarNPS2D(matriz, matriz, cal, npsParams, variance);

        // ---- geometry in pixels (y increases downward; "below" = +y) ----
        double marrowR = (marrowDiamMM / 2.0) / pixelSize;
        double ringOuterR = (ringOuterDiamMM / 2.0) / pixelSize;
        double halfBase = (triBaseMM / 2.0) / pixelSize;
        double baseY = cyPx + (triBaseBelowMM / pixelSize);           // base line (bottom of wings)
        double apexRise = halfBase * Math.tan(Math.toRadians(triBaseAngleDeg));
        double apexY = baseY - apexRise;                              // top vertex of the triangle
        // triangle vertices: left/right base corners and top apex (centred on x)
        double axL = cxPx - halfBase, axR = cxPx + halfBase;
        double apex = cxPx;

        FloatProcessor ideal = new FloatProcessor(matriz, matriz);
        final int SS = 4;                 // sub-pixel supersampling per axis
        final double inv = 1.0 / SS;

        for (int y = 0; y < matriz; y++) {
            for (int x = 0; x < matriz; x++) {
                int marrowHits = 0, boneHits = 0, total = SS * SS;
                for (int sy = 0; sy < SS; sy++) {
                    for (int sx = 0; sx < SS; sx++) {
                        double px = x + (sx + 0.5) * inv - 0.5;
                        double py = y + (sy + 0.5) * inv - 0.5;
                        double dx = px - cxPx, dy = py - cyPx;
                        double r = Math.hypot(dx, dy);
                        boolean inTri = pointInTriangle(px, py, axL, baseY, axR, baseY, apex, apexY);
                        if (r <= marrowR) {
                            marrowHits++;                       // marrow interior
                        } else if (r <= ringOuterR || inTri) {
                            boneHits++;                         // vertebra ring or wings
                        }
                    }
                }
                double val;
                if (marrowHits == 0 && boneHits == 0) {
                    val = bg;
                } else {
                    // coverage-weighted blend (bone takes precedence where it overlaps)
                    double fBone = boneHits / (double) total;
                    double fMarrow = marrowHits / (double) total;
                    double fBg = 1.0 - fBone - fMarrow;
                    if (fBg < 0) fBg = 0;
                    val = fBone * bone + fMarrow * marrow + fBg * bg;
                }
                ideal.putPixelValue(x, y, val);
            }
        }

        ImageProcessor ipMtf = MTFModel.applyMTF(ideal, mtf2D);
        ImageProcessor ipRuido = npsModel.texturedNoiseIP(nps2D, pixelSize);
        FloatProcessor sliceFim = adicionarRuido(ipMtf, ipRuido);
        stack.addSlice(criarFatiaShortUnscaled(sliceFim));

        ImagePlus imp = new ImagePlus("SpinePhantom_" + kernelName, stack);
        imp.setCalibration(cal);
        return imp;
    }

    /** Barycentric point-in-triangle test. */
    private static boolean pointInTriangle(double px, double py,
            double x1, double y1, double x2, double y2, double x3, double y3) {
        double d = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
        if (d == 0) return false;
        double a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / d;
        double b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / d;
        double c = 1 - a - b;
        return a >= 0 && b >= 0 && c >= 0;
    }

    /**
     * Soma uma imagem determinística (fundo) a uma imagem estocástica (ruído)
     */
    private static FloatProcessor adicionarRuido(ImageProcessor base, ImageProcessor ruido) {
        FloatProcessor result = (FloatProcessor) base.duplicate().convertToFloatProcessor();
        float[] pBase = (float[]) result.getPixels();
        float[] pRuido = (float[]) ruido.getPixels();
        
        for (int i = 0; i < pBase.length; i++) {
            pBase[i] += pRuido[i];
        }
        return result;
    }

    /**
     * Converte um FloatProcessor para ShortProcessor rigorosamente sem aplicar auto-scaling.
     * Garante que os valores brutos calibrados (+8192) sejam preservados.
     */
    private static ij.process.ShortProcessor criarFatiaShortUnscaled(FloatProcessor fp) {
        float[] pixelsFloat = (float[]) fp.getPixels();
        short[] pixelsShort16 = new short[pixelsFloat.length];
        
        for (int i = 0; i < pixelsFloat.length; i++) {
            int val = Math.round(pixelsFloat[i]);
            if (val < 0) val = 0;
            if (val > 65535) val = 65535;
            pixelsShort16[i] = (short) val;
        }
        return new ij.process.ShortProcessor(fp.getWidth(), fp.getHeight(), pixelsShort16, null);
    }
}
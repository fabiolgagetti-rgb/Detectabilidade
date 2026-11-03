package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.List;

import br.unifesp.fisicamedica.Segmenter.DetectorResult;

/**
 * Staged TTF error-attribution harness.
 *
 * Builds the TTF estimation pipeline one layer at a time so each stage adds
 * exactly ONE new source of error. The marginal change between consecutive
 * stages (printed as d(Si-Sj)) shows which part of the pipeline contributes
 * most to the TTF error.
 *
 *   Stage 0  MODEL LSF  -> window -> FFT
 *            Isolates the LSF-window + FFT/zero-pad stage only.
 *            Analytic ideal LSF; no ESF, no derivative, no interpolation, no noise.
 *
 *   Stage 1  MODEL ESF  -> first-difference + window (derivate) -> FFT
 *            Adds the finite-difference derivative on top of Stage 0.
 *            Analytic ideal ESF; still no interpolation, no noise.
 *
 *   Stage 2  REAL ESF from a PERFECT edge -> average -> derivate -> FFT
 *            Adds Lanczos interpolation + profile sampling/averaging along the
 *            TRUE (known) circular boundary. Noise enters here, because
 *            calculateESF samples the actual rendered (textured-noise) image.
 *
 *   Stage 3  REAL ESF from a DETECTED edge -> average -> derivate -> FFT
 *            Full pipeline: adds edge detection / sub-pixel registration
 *            (Segmenter.edgeFinder) on top of Stage 2.
 *
 * Noise only affects Stages 2 and 3 (they sample the noisy image). Stages 0 and
 * 1 are analytic, noise-free references: their PE is constant across SD rows.
 *
 * Metric: PE = 100 * sum|measured-model| / sum(model) over 0..pixel-Nyquist,
 * identical to TtfValidation, so all stages are directly comparable.
 *
 * NOTE: Stage 0's windowing replicates derivate()'s Tukey window (which couples
 * first-difference + window and so cannot be reused for a difference-free LSF).
 * Keep it in sync with TtfCalc.derivate() if that window changes.
 */
public class TtfContributionTest {

    public static void main(String[] args) {
        int matriz = 512;
        double fov = 250.0;
        double pixelSize = fov / matriz;
        double bgHU = -1000.0;
        double objectMeanHU = 0.0;
        double diametroMM = 200.0;
        double[] insertHU = {100};
        double insertRaioMM = 12.5;

        Calibration cal = new Calibration();
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.setUnit("mm");

        int roiSize = (int) (insertRaioMM * 3 / pixelSize);
        int centro = matriz / 2;
        Roi userRoi = new Roi(centro - roiSize / 2, centro - roiSize / 2, roiSize, roiSize);

        int fftPoints = 64;
        double nyquistPixel = 1.0 / (2.0 * pixelSize);
        double df = 2.0 * nyquistPixel / fftPoints;      // model-grid spacing
        double dxEsf = 0.5 * pixelSize;                  // supersampled ESF step (mm)

        String[] kernels = {"SOFT", "LUNG"};
        int[] noiseLevelsSD = {0, 5, 10, 20, 30};
        int reps = Integer.getInteger("opt.reps", 20);

        System.out.println("# TtfContributionTest  window=" + System.getProperty("ttf.window", "current")
                + "  deconv=" + Boolean.getBoolean("ttf.deconv"));
        System.out.println("# PE = integral |measured-model| / integral model (%), up to pixel Nyquist");
        System.out.println("Kernel\tSD\tS0_modelLSF\tS1_modelESF\tS2_perfectEdge\tS3_detectedEdge"
                + "\t|\td(S1-S0)\td(S2-S1)\td(S3-S2)\tFails");

        for (String kernel : kernels) {
            double[][] mtfParams = MTFModel.getMTFparameters(kernel);

            // model TTF on the comparison grid (ground truth)
            double[] model = new double[fftPoints];
            for (int i = 0; i < fftPoints; i++) model[i] = MTFModel.getMTFValue(i * df, mtfParams);
            double dcm = model[0];
            if (dcm > 0) for (int i = 0; i < fftPoints; i++) model[i] /= dcm;

            // analytic ideal LSF/ESF on the supersampled grid (length matches real ESF)
            int esfLen = 2 * ((int) (10 / pixelSize)) * 2 + 1; // mirror calculateESF length
            double[] idealLsf = buildIdealLsf(mtfParams, dxEsf, nyquistPixel, esfLen);
            double[] idealEsf = integrate(idealLsf);

            // Stage 0: ideal LSF -> window -> FFT  (no difference)
            double[] s0mtf = TtfCalc.calculateTtf(tukeyLikeDerivateWindow(idealLsf), 2 * fftPoints);
            double s0 = pe(s0mtf, model, df, nyquistPixel);

            // Stage 1: ideal ESF -> derivate (difference+window) -> FFT
            double[] s1mtf = TtfCalc.calculateTtf(TtfCalc.derivate(idealEsf), 2 * fftPoints);
            double s1 = pe(s1mtf, model, df, nyquistPixel);

            for (int sd : noiseLevelsSD) {
                double variance = sd * sd;
                List<Double> s2list = new ArrayList<>();
                List<Double> s3list = new ArrayList<>();
                int fails = 0;

                for (int r = 0; r < reps; r++) {
                    try {
                        ImagePlus imp = new PhantomImageCreator().criarImagemCTPhantom(
                                0, 1, matriz, fov, bgHU, objectMeanHU, diametroMM, 0.0,
                                insertHU, insertRaioMM, 0.0, kernel, variance);
                        ImageStack stack = imp.getStack();
                        ImageProcessor ip = stack.getProcessor(1);

                        // Stage 2: perfect circular edge, real ESF sampling
                        FloatPolygon perfect = perfectCircle(centro, centro, insertRaioMM / pixelSize);
                        double s2 = ttfFromBoundary(ip, cal, perfect, fftPoints, df, nyquistPixel, model);
                        if (!Double.isNaN(s2)) s2list.add(s2);

                        // Stage 3: detected edge, full pipeline
                        DetectorResult dR = Segmenter.edgeFinder(stack, cal, userRoi, false, true);
                        if (dR == null) { fails++; r--; continue; }
                        TtfCalc.mtfResult res = TtfCalc.analisarStackComInserts(stack, cal, dR, fftPoints);
                        if (res == null || res.ttfs == null || res.ttfs.isEmpty()) { fails++; continue; }
                        double[] s3mtf = res.ttfs.values().iterator().next();
                        double s3 = pe(s3mtf, model, df, nyquistPixel);
                        if (!Double.isNaN(s3)) s3list.add(s3);
                    } catch (Exception e) {
                        fails++;
                    }
                }

                double s2 = mean(s2list);
                double s3 = mean(s3list);
                System.out.printf("%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t|\t%+.3f\t%+.3f\t%+.3f\t%d%n",
                        kernel, sd, s0, s1, s2, s3,
                        (s1 - s0), (s2 - s1), (s3 - s2), fails);
            }
        }
        System.exit(0);
    }

    // ------------------------------------------------------------------
    // Stage building blocks
    // ------------------------------------------------------------------

    /** Ideal LSF on the supersampled grid via inverse cosine transform of the model MTF. */
    private static double[] buildIdealLsf(double[][] mtfParams, double dx, double fNy, int N) {
        double[] lsf = new double[N];
        int M = 600;
        double dcModel = MTFModel.getMTFValue(0.0, mtfParams);
        for (int i = 0; i < N; i++) {
            double x = (i - N / 2) * dx;
            double s = 0.0, prevF = 0.0, prevG = 1.0; // normalized MTF(0)=1
            for (int k = 1; k < M; k++) {
                double f = fNy * k / (M - 1);
                double mtf = MTFModel.getMTFValue(f, mtfParams) / dcModel;
                double g = mtf * Math.cos(2 * Math.PI * f * x);
                s += 0.5 * (prevG + g) * (f - prevF);
                prevF = f; prevG = g;
            }
            lsf[i] = s;
        }
        return lsf;
    }

    /** Cumulative sum: ideal LSF -> ideal ESF on the same grid. */
    private static double[] integrate(double[] lsf) {
        double[] esf = new double[lsf.length];
        double acc = 0.0;
        for (int i = 0; i < lsf.length; i++) { acc += lsf[i]; esf[i] = acc; }
        return esf;
    }

    /**
     * Replicates the Tukey window applied inside TtfCalc.derivate(), but applied
     * DIRECTLY to a supplied LSF (no first-difference). Used only for Stage 0 so
     * it isolates the window+FFT from the finite-difference. Mirror of the window
     * branch in derivate(); keep in sync if that changes.
     */
    private static double[] tukeyLikeDerivateWindow(double[] lsf) {
        int sizeLsf = lsf.length;
        double max = 0;
        for (double v : lsf) if (Math.abs(v) > max) max = Math.abs(v);
        int left = 0, right = 0; boolean peak = false;
        for (int i = 0; i < sizeLsf; i++) {
            if (Math.abs(lsf[i]) > max / 2 && !peak) { left = i; peak = true; }
            if (Math.abs(lsf[i]) < max / 2 && peak)  { right = i; break; }
        }
        double FWHM = right - left;
        double alpha = 1 - (FWHM * 4.0 / sizeLsf);
        double steepness = 5;
        double[] out = new double[sizeLsf];
        if (alpha <= 0) {
            System.arraycopy(lsf, 0, out, 0, sizeLsf);
        } else if (alpha >= 1.0) {
            for (int i = 0; i < sizeLsf; i++) {
                double x = (double) i / (sizeLsf - 1);
                out[i] = lsf[i] * 0.5 * (1 - Math.cos(2 * Math.PI * x));
            }
        } else {
            double halfAlpha = alpha / 2.0;
            for (int i = 0; i < sizeLsf; i++) {
                double x = (double) i / (sizeLsf - 1);
                double w;
                if (x < halfAlpha) {
                    double t = Math.pow(x / halfAlpha, steepness);
                    w = 0.5 * (1 + Math.cos(Math.PI * (t - 1)));
                } else if (x <= 1.0 - halfAlpha) {
                    w = 1.0;
                } else {
                    double t = Math.pow((1.0 - x) / halfAlpha, steepness);
                    w = 0.5 * (1 + Math.cos(Math.PI * (t - 1)));
                }
                out[i] = lsf[i] * w;
            }
        }
        return out;
    }

    /** TTF from a given boundary polygon using the real ESF sampling + averaging. */
    private static double ttfFromBoundary(ImageProcessor ip, Calibration cal, FloatPolygon boundary,
                                          int fftPoints, double df, double fNy, double[] model) {
        int profLen = (int) (10 / cal.pixelWidth);
        // expectedContrasts=null -> no contrast gating; every profile is used,
        // isolating interpolation+averaging from the contrast-matching logic.
        List<List<double[]>> grouped = TtfCalc.calculateESF(cal, ip, boundary, profLen, false, null);
        if (grouped == null || grouped.isEmpty() || grouped.get(0).isEmpty()) return Double.NaN;
        double[] avgEsf = CurveAverage.averageEsfs(grouped.get(0));
        if (avgEsf == null) return Double.NaN;
        double[] mtf = TtfCalc.calculateTtf(TtfCalc.derivate(avgEsf), 2 * fftPoints);
        return pe(mtf, model, df, fNy);
    }

    /** Sub-pixel polygon tracing the true insert circle (perfect edge for Stage 2). */
    private static FloatPolygon perfectCircle(double cx, double cy, double rPix) {
        FloatPolygon fp = new FloatPolygon();
        int n = Math.max(360, (int) (2 * Math.PI * rPix * 2));
        for (int i = 0; i < n; i++) {
            double a = 2 * Math.PI * i / n;
            fp.addPoint((float) (cx + rPix * Math.cos(a)), (float) (cy + rPix * Math.sin(a)));
        }
        return fp;
    }

    // ------------------------------------------------------------------
    // Metric / stats
    // ------------------------------------------------------------------

    private static double pe(double[] measured, double[] model, double df, double nyquist) {
        double sumDiff = 0.0, sumModel = 0.0;
        int limit = Math.min(measured.length, model.length);
        for (int i = 0; i < limit; i++) {
            double f = i * df;
            if (f > nyquist) break;
            sumDiff += Math.abs(measured[i] - model[i]);
            sumModel += model[i];
        }
        return sumModel == 0 ? Double.NaN : 100.0 * sumDiff / sumModel;
    }

    private static double mean(List<Double> v) {
        if (v.isEmpty()) return Double.NaN;
        double s = 0; for (double d : v) s += d; return s / v.size();
    }
}

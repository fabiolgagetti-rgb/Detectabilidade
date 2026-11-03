package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.process.FloatProcessor;

import java.util.ArrayList;
import java.util.List;

import br.unifesp.fisicamedica.Segmenter.DetectorResult;

/**
 * End-to-end detectability validation.
 *
 * Mirrors the component validations (TtfValidation, NpsValidation): it drives
 * the FULL measured pipeline and compares the resulting detectability index d'
 * against the THEORETICAL d' computed from the analytic MTF/NPS models (the
 * ground truth), across kernels and noise levels.
 *
 * Both d' values go through the SAME DetectabilityIndexCalculator with the SAME
 * task/observer parameters, so the only difference is measured vs model
 * MTF/NPS. This surfaces any residual pipeline issue that only shows up once
 * TTF and NPS are combined into a task-based figure of merit (units, frequency
 * grids, normalization, interpolation limits, ...).
 *
 * Measured inputs per trial:
 *   - measured TTF: TtfCalc.analisarStackComInserts on an insert slice
 *     (same path as TtfValidation / the TtfPanel).
 *   - measured NPS: NPScalc.calculate2DNPS + radialProfile on homogeneous ROIs
 *     (same path as NpsValidation, including 2nd-order detrending).
 *
 * Report: Mean/SD of PE_d' = 100*(d'_measured - d'_theory)/d'_theory, and the
 * mean measured/theoretical d' for context, plus failure counts.
 */
public class DetectabilityValidation {

    public static void main(String[] args) {
        // ---- phantom / geometry (match TtfValidation's insert geometry) ----
        int matriz = 512;
        double fov = 250.0;
        double pixelSize = fov / matriz;
        double bgHU = -1000.0;
        double objectMeanHU = 0.0;
        double diametroMM = 200.0;
        double trendAmpHU = 0.0;
        double[] insertHU = {100};
        double insertRaioMM = 12.5;
        double insertDistanciaMM = 0.0;

        Calibration cal = new Calibration();
        cal.pixelWidth = pixelSize;
        cal.pixelHeight = pixelSize;
        cal.setUnit("mm");

        int roiSize = (int) (insertRaioMM * 3 / pixelSize);
        int centro = matriz / 2;
        Roi userRoi = new Roi(centro - roiSize / 2, centro - roiSize / 2, roiSize, roiSize);

        // ---- NPS measurement geometry (match NpsValidation) ----
        int npsMatriz = 64;
        double npsFov = npsMatriz * pixelSize;
        Calibration npsCal = new Calibration();
        npsCal.pixelWidth = pixelSize;
        npsCal.pixelHeight = pixelSize;
        npsCal.setUnit("mm");
        int numNpsRois = 20;

        // ---- task / observer parameters (identical for measured and theory) ----
        double taskRadius = 3.0;      // mm
        double taskContrast = 15.0;   // HU
        double viewDistance = 500.0;  // mm
        double displaySize = 250.0;   // mm
        double dFov = 250.0;          // mm (display FOV)

        int fftPoints = 64;
        double df = 2.0 * (1.0 / (2.0 * pixelSize)) / fftPoints;  // TTF model grid spacing

        String[] kernels = {"SOFT", "LUNG"};
        int[] noiseLevelsSD = {5, 10, 20, 30};
        int numMeasurements = 20;
        boolean deconv = Boolean.getBoolean("ttf.deconv");

        System.out.println("# DetectabilityValidation  deconv=" + deconv
                + "  window=" + System.getProperty("ttf.window", "current")
                + "  task: C=" + taskContrast + "HU r=" + taskRadius + "mm");
        System.out.println("Kernel\tSD\td'_theory\td'_meas\tMean_PEd(%)\tSD_PEd(%)\tFails");

        for (String kernel : kernels) {
            double[][] mtfParams = MTFModel.getMTFparameters(kernel);
            double[] npsParams = NPSModel.getNPSCoefficients(kernel);

            for (int sd : noiseLevelsSD) {
                double variance = sd * sd;

                // ---------------- theoretical d' (ground truth) ----------------
                // Model MTF on the TTF grid; model NPS on the NPS radial grid.
                double[] modelTtfFreq = new double[fftPoints];
                double[] modelTtf = new double[fftPoints];
                for (int i = 0; i < fftPoints; i++) {
                    modelTtfFreq[i] = i * df;
                    modelTtf[i] = MTFModel.getMTFValue(i * df, mtfParams);
                }
                double dc = modelTtf[0];
                if (dc > 0) for (int i = 0; i < fftPoints; i++) modelTtf[i] /= dc;

                double npsDf = 1.0 / (npsMatriz * pixelSize);
                double npsNyq = 1.0 / (2.0 * pixelSize);
                int numBins = npsMatriz / 2;
                double[] modelNpsFreq = new double[numBins];
                double[] modelNps = new double[numBins];
                for (int i = 0; i < numBins; i++) {
                    double f = i * npsDf;
                    modelNpsFreq[i] = f;
                    modelNps[i] = (f <= npsNyq) ? NPSModel.getNPSValue(f, npsParams, variance) : 0.0;
                }

                double dPrimeTheory = new DetectabilityIndexCalculator(
                        modelTtf, modelTtfFreq, modelNps, modelNpsFreq,
                        taskRadius, taskContrast, viewDistance, displaySize, dFov).compute();

                // ---------------- measured d' (full pipeline) ----------------
                List<Double> peList = new ArrayList<>();
                List<Double> measList = new ArrayList<>();
                int failures = 0;

                for (int m = 0; m < numMeasurements; m++) {
                    try {
                        // (1) measured TTF from an insert phantom
                        ImagePlus impTtf = new PhantomImageCreator().criarImagemCTPhantom(
                                0, 1, matriz, fov, bgHU, objectMeanHU, diametroMM, trendAmpHU,
                                insertHU, insertRaioMM, insertDistanciaMM, kernel, variance);
                        DetectorResult dR = Segmenter.edgeFinder(impTtf.getStack(), cal, userRoi, false, true);
                        if (dR == null) { failures++; m--; continue; }
                        TtfCalc.mtfResult ttfRes =
                                TtfCalc.analisarStackComInserts(impTtf.getStack(), cal, dR, fftPoints);
                        if (ttfRes == null || ttfRes.ttfs == null || ttfRes.ttfs.isEmpty()) { failures++; continue; }
                        Integer id = ttfRes.ttfs.keySet().iterator().next();
                        double[] ttfData = ttfRes.ttfs.get(id);
                        double[] ttfFreq = ttfRes.freqs.get(id);

                        // (2) measured NPS from homogeneous ROIs (NpsValidation path)
                        ImagePlus impNps = new PhantomImageCreator().criarImagemCTPhantom(
                                numNpsRois, 0, npsMatriz, npsFov, 0.0, 0.0, 500, 0.0,
                                new double[0], 0.0, 0.0, kernel, variance);
                        ImageStack nStack = impNps.getStack();
                        List<double[][]> all2DNps = new ArrayList<>();
                        for (int s = 1; s <= numNpsRois; s++) {
                            FloatProcessor rec = (FloatProcessor) nStack.getProcessor(s).convertToFloat();
                            float[][] roiData = rec.getFloatArray();
                            double[][] detr = Detrending.detrendWith2ndOrderPolynomial(roiData);
                            float[][] detrF = new float[rec.getHeight()][rec.getWidth()];
                            for (int yy = 0; yy < rec.getHeight(); yy++)
                                for (int xx = 0; xx < rec.getWidth(); xx++)
                                    detrF[yy][xx] = (float) detr[yy][xx];
                            rec = new FloatProcessor(detrF);
                            all2DNps.add(NPScalc.calculate2DNPS(rec, npsCal, npsMatriz));
                        }
                        double[][] avg2D = mean2D(all2DNps, npsMatriz);
                        double[] npsFreq = new double[npsMatriz / 2];
                        double[] npsData = NPScalc.radialProfile(avg2D, npsCal, npsFreq);

                        // (3) measured d' through the same calculator / same task.
                        // Diagnostic isolation: optionally substitute the model TTF or
                        // model NPS to attribute any d' deficit to one factor.
                        //   -Dval.useModelTtf=true -> measured NPS only
                        //   -Dval.useModelNps=true -> measured TTF only
                        double[] useTtf = ttfData, useTtfF = ttfFreq, useNps = npsData, useNpsF = npsFreq;
                        if (Boolean.getBoolean("val.useModelTtf")) { useTtf = modelTtf; useTtfF = modelTtfFreq; }
                        if (Boolean.getBoolean("val.useModelNps")) { useNps = modelNps; useNpsF = modelNpsFreq; }
                        double dPrimeMeas = new DetectabilityIndexCalculator(
                                useTtf, useTtfF, useNps, useNpsF,
                                taskRadius, taskContrast, viewDistance, displaySize, dFov).compute();

                        if (Double.isNaN(dPrimeMeas) || Double.isInfinite(dPrimeMeas)) { failures++; continue; }
                        measList.add(dPrimeMeas);
                        if (dPrimeTheory != 0.0)
                            peList.add(100.0 * (dPrimeMeas - dPrimeTheory) / dPrimeTheory);
                    } catch (Exception e) {
                        failures++;
                    }
                }

                double meanPE = mean(peList);
                double sdPE = sd(peList, meanPE);
                double meanMeas = mean(measList);

                System.out.printf("%s\t%d\t%.4f\t%.4f\t%.3f\t%.3f\t%d%n",
                        kernel, sd, dPrimeTheory, meanMeas, meanPE, sdPE, failures);
            }
        }
        System.exit(0);
    }

    private static double[][] mean2D(List<double[][]> list, int size) {
        double[][] avg = new double[size][size];
        int n = list.size();
        for (double[][] ps : list)
            for (int y = 0; y < size; y++)
                for (int x = 0; x < size; x++)
                    avg[y][x] += ps[y][x];
        if (n > 0)
            for (int y = 0; y < size; y++)
                for (int x = 0; x < size; x++)
                    avg[y][x] /= n;
        return avg;
    }

    private static double mean(List<Double> v) {
        if (v.isEmpty()) return 0;
        double s = 0; for (double d : v) s += d; return s / v.size();
    }

    private static double sd(List<Double> v, double mean) {
        if (v.isEmpty()) return 0;
        double s = 0; for (double d : v) s += (d - mean) * (d - mean);
        return Math.sqrt(s / v.size());
    }
}

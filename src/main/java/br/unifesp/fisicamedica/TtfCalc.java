package br.unifesp.fisicamedica;

import java.awt.Point;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jtransforms.fft.DoubleFFT_1D;

import br.unifesp.fisicamedica.Segmenter.DetectorResult;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;

public class TtfCalc {
    static int phantomWidth = 0;
    static int phantomHeight = 0;
    
    public static class mtfResult {
        public Map<Integer, double[]> ttfs = new HashMap<>();
        public Map<Integer, double[]> freqs = new HashMap<>();
        public Map<Integer, double[]> lsfs = new HashMap<>();
        public Map<Integer, double[]> esfs = new HashMap<>();
        public Map<Integer, Double> contrasts = new HashMap<>();
        public Map<Integer, Double> I0 = new HashMap<>();
        public Map<Integer, Double> I1 = new HashMap<>();
        public Map<Integer, Point> positions = new HashMap<>();
    }
    
    public static mtfResult analisarStackComInserts(ImageStack stack, Calibration cal, DetectorResult dR, int fftPoints) {
        int nSlices = stack.getSize();
        List<Point> points = dR.Pos;
        mtfResult result = new mtfResult();
        
        for (int i = 0; i < points.size(); i++) {
            if (dR.rois.get(i) == null) continue;
            
            boolean invert = false;
            if (dR.HU.get(i) < dR.Threshold.get(i)) invert = true;
            
            double targetContrast = 0.0;
            double I0 = 0.0;
            double I1 = 0.0;

            if (!dR.modeContrast.isEmpty() && i < dR.modeContrast.size() && !dR.modeContrast.get(i).isEmpty()) {
                targetContrast = dR.modeContrast.get(i).get(0);
                I0 = dR.meanI0.get(i).get(0);
                I1 = dR.meanI1.get(i).get(0);

            }
            
            List<Double> singleExpectedContrastList = new ArrayList<>();
            singleExpectedContrastList.add(targetContrast);
            
            List<double[]> esfProfilesAccumulated = new ArrayList<>();
            
            for (int j = 1; j <= nSlices; j++) {
                ImageProcessor ip = stack.getProcessor(j);
                Roi roi = dR.rois.get(i).get(j - 1);
                FloatPolygon boundary = roi.getFloatPolygon();

                double[][] pointNormals = null;
                if (dR.edgeNormals != null && i < dR.edgeNormals.size()
                        && dR.edgeNormals.get(i) != null && (j - 1) < dR.edgeNormals.get(i).size()) {
                    pointNormals = dR.edgeNormals.get(i).get(j - 1);
                }

                List<List<double[]>> groupedEsfs = calculateESF(cal, ip, boundary, (int) (10 / cal.pixelWidth), invert, singleExpectedContrastList, pointNormals);
                
                if (groupedEsfs != null && !groupedEsfs.isEmpty() && groupedEsfs.get(0) != null) {
                    esfProfilesAccumulated.addAll(groupedEsfs.get(0));
                }
            }
            
            if (esfProfilesAccumulated.isEmpty()) continue; 
            
            double[] avgEsf = CurveAverage.averageEsfs(esfProfilesAccumulated);
            if (avgEsf == null) continue;

            double[] lsf = derivate(avgEsf);
            double[] mtf = calculateTtf(lsf, 2 * fftPoints);

            // ----------------------------------------------------------------
            // TODO(Problem A): deconvolve the measurement-chain transfer
            // functions out of `mtf` to remove the systematic low bias
            // (~5% on edge-enhancing kernels, worst in the mid-band).
            //
            // WHY THIS IS MATHEMATICALLY JUSTIFIED
            // ------------------------------------
            // The measured TTF is the true system MTF multiplied by the
            // transfer functions of the operators we apply when estimating it.
            // Each operator below is linear and shift-invariant along the
            // supersampled ESF axis, so in the Fourier domain they are exact
            // multiplicative factors H_i(f). The pipeline computes
            //
            //     ESF --(1) Lanczos resample at step dx --(2) first difference
            //         --(3) Tukey window--> LSF,   MTF = |FFT(LSF)| / FFT[0].
            //
            // By the convolution theorem the measured transform factorises as
            //
            //     MTF_meas(f) = MTF_true(f) * H_lanczos(f) * H_diff(f) * (window bias)
            //
            // where the first two factors are KNOWN in closed form:
            //
            //  (2) First difference lsf[i] = esf[i+1]-esf[i] is a convolution
            //      with the kernel {+1,-1} at spacing dx. Its transfer function
            //      is the forward-difference operator
            //          D(f)       = (e^{+j2*pi*f*dx} - 1)/dx            (complex),
            //          |D(f)|      = (2/dx) * |sin(pi*f*dx)|.
            //      The IDEAL derivative we intend to approximate has transfer
            //      |D_ideal(f)| = 2*pi*f. Hence the finite-difference imposes a
            //      frequency-dependent gain
            //          H_diff(f)  = |sin(pi*f*dx)| / (pi*f*dx)  = sinc(f*dx),
            //      i.e. a sinc droop: ~1% at f=0.3 lp/mm, ~10% near Nyquist for
            //      dx = 0.5*pixelWidth. (The linear-phase factor only recentres
            //      the LSF and drops out of the magnitude |FFT|.)
            //
            //  (1) Lanczos-3 interpolation is convolution with the Lanczos
            //      kernel L_a(t) = sinc(t)*sinc(t/a), a=3, sampled at the ESF
            //      step. Its transfer function H_lanczos(f) is the (known)
            //      Fourier transform of that kernel evaluated on [0, f_Nyquist];
            //      it is ~1 at low f and rolls off toward Nyquist, accounting
            //      for the high-frequency part of the observed droop.
            //
            // THE CORRECTION
            // --------------
            // Because the factors are multiplicative and non-zero over the band
            // of interest (0 .. pixel Nyquist, well below the first sinc zero at
            // f = 1/dx), the unbiased estimate is a straight division:
            //
            //     MTF_corrected(f) = MTF_meas(f) / (H_diff(f) * H_lanczos(f)).
            //
            // This is a deconvolution by pointwise division in the frequency
            // domain. It is safe here (no noise amplification blow-up) because
            // both factors stay comfortably above ~0.9 across the reported band
            // and we never divide near a transfer-function zero. The window bias
            // is handled separately by the side-lobe-aware window in derivate();
            // only the two analytic, kernel-independent factors are divided out
            // here. Implement as: build H_diff(f) and H_lanczos(f) on the same
            // frequency grid as `frequency` below, then divide `mtf` by their
            // product before storing it.
            // ----------------------------------------------------------------

            double samplingStepPx = 0.5; 
            double dx_mm = samplingStepPx * cal.pixelWidth; 
            double fNyquist = 1.0 / (2.0 * dx_mm); 
            int length = mtf.length;
             
            double[] frequency = new double[length];
            for (int i1 = 0; i1 < length; i1++) {
                frequency[i1] = (fNyquist * i1) / (length - 1);
            }
            
            int insertId = i;
            
            result.esfs.put(insertId, avgEsf);
            result.lsfs.put(insertId, lsf);
            result.ttfs.put(insertId, mtf);
            result.freqs.put(insertId, frequency);
            result.contrasts.put(insertId, targetContrast);
            result.I0.put(insertId, I0);
            result.I1.put(insertId, I1);

            result.positions.put(insertId, points.get(i));
        }
        return result;
    }

    public static List<List<double[]>> calculateESF(Calibration cal, ImageProcessor ip, FloatPolygon boundary, int profileLength, boolean invert, List<Double> expectedContrasts) {
        return calculateESF(cal, ip, boundary, profileLength, invert, expectedContrasts, null);
    }

    /**
     * ESF sampling with an optional per-point normal array (Option A).
     *
     * When pointNormals is non-null, aligned 1:1 with the boundary polygon, and
     * -Dttf.paeNormals=true, the profile at each boundary point is sampled along
     * that point's own analytic (PAE) normal instead of a least-squares fit over
     * neighbouring polygon points. The neighbour-fit is sensitive to gaps in the
     * polygon ("holes" left when a wand point has no nearby PAE edge): a hole makes
     * distant vertices adjacent and skews the fitted tangent/normal. A per-point
     * analytic normal is independent of polygon topology and so is immune to holes.
     * Falls back to the neighbour-fit when normals are unavailable for a point.
     */
    public static List<List<double[]>> calculateESF(Calibration cal, ImageProcessor ip, FloatPolygon boundary, int profileLength, boolean invert, List<Double> expectedContrasts, double[][] pointNormals) {
        if (boundary == null) throw new IllegalArgumentException("null ROI");
        // Option C (ttf.paePrimary) builds the polygon from PAE points and implies
        // using their analytic per-point normals; Option A (ttf.paeNormals) enables
        // the same normal use on the legacy wand-anchored polygon.
        boolean usePaeNormals = (Boolean.getBoolean("ttf.paeNormals") || Boolean.getBoolean("ttf.paePrimary")
                || Boolean.getBoolean("ttf.paeFill"))
                && pointNormals != null && pointNormals.length == boundary.npoints;

        // Optional smoothing of the PAE normals along the edge ordering. PAE normals
        // are local and hole-immune but per-point noisy; averaging each with its
        // neighbours (a +/- window, in polygon index order) recovers the noise
        // reduction the neighbour-fit provided, WITHOUT reintroducing hole
        // sensitivity (a dropped point just means one fewer term in the average).
        // -Dttf.paeNormalSmooth = half-window in points (0 = raw PAE normals).
        double[][] normals = pointNormals;
        if (usePaeNormals) {
            int halfWin = Integer.getInteger("ttf.paeNormalSmooth", 3);
            if (halfWin > 0) {
                int n = pointNormals.length;
                normals = new double[n][2];
                for (int p = 0; p < n; p++) {
                    double sx = 0, sy = 0;
                    for (int w = -halfWin; w <= halfWin; w++) {
                        int q = ((p + w) % n + n) % n;      // wrap around closed contour
                        sx += pointNormals[q][0];
                        sy += pointNormals[q][1];
                    }
                    double len = Math.hypot(sx, sy);
                    if (len == 0) { normals[p][0] = pointNormals[p][0]; normals[p][1] = pointNormals[p][1]; }
                    else { normals[p][0] = sx / len; normals[p][1] = sy / len; }
                }
            }
        }
        profileLength *= 2;
        if (profileLength % 2 == 0) profileLength++;

        int half = profileLength / 2;
        double step = 0.5;
        // Contrast-matching tolerance. Kept as an ABSOLUTE HU floor (rejects profiles
        // corrupted by internal structures, e.g. bone marrow inside a spine insert,
        // which shift the endpoint contrast) PLUS a RELATIVE band that scales with the
        // expected contrast. The relative term bridges the systematic gap between the
        // segmenter's contrast estimate and the ESF endpoint contrast on edge-enhancing
        // kernels, without widening the window for low-contrast structures.
        double toleranceHU = Double.parseDouble(System.getProperty("ttf.tolHU", "15"));
        double relTol      = Double.parseDouble(System.getProperty("ttf.relTol", "0.20"));
        
        List<List<double[]>> groupedProfiles = new ArrayList<>();
        boolean checkContrast = expectedContrasts != null && !expectedContrasts.isEmpty() && expectedContrasts.get(0) > 0;
        
        if (!checkContrast) {
            groupedProfiles.add(new ArrayList<>()); 
        } else {
            for (int i = 0; i < expectedContrasts.size(); i++) {
                groupedProfiles.add(new ArrayList<>());
            }
        }

        for (int i = 0; i < boundary.npoints - 1; i++) {
            float x0 = boundary.xpoints[i];
            float y0 = boundary.ypoints[i];

            double nx, ny;

            if (usePaeNormals) {
                // Option A: use this point's own analytic PAE normal; independent of
                // polygon neighbours, so holes in the polygon cannot skew it.
                double pnx = normals[i][0];
                double pny = normals[i][1];
                double plen = Math.hypot(pnx, pny);
                if (plen == 0) continue;   // degenerate normal; skip (one fewer profile)
                nx = pnx / plen;
                ny = pny / plen;
                if (invert) { nx = -nx; ny = -ny; }
                double[] profile = new double[profileLength];
                for (int j = -half; j <= half; j++) {
                    double x = x0 + j * step * nx;
                    double y = y0 + j * step * ny;
                    double val = getLanczosInterpolatedPixel(ip, x, y);
                    profile[j + half] = cal.getCValue(val);
                }
                addProfileWithContrastCheck(profile, profileLength, checkContrast,
                        expectedContrasts, toleranceHU, relTol, groupedProfiles);
                continue;
            }

            double searchRadiusMM = 5.0;
            double searchRadiusPx = searchRadiusMM / cal.pixelWidth;

            int startIdx = i;
            double distBack = 0;
            while (startIdx > 0 && distBack < searchRadiusPx) {
                distBack += Math.hypot(
                    boundary.xpoints[startIdx] - boundary.xpoints[startIdx - 1], 
                    boundary.ypoints[startIdx] - boundary.ypoints[startIdx - 1]
                );
                startIdx--;
            }

            int endIdx = i;
            double distFwd = 0;
            while (endIdx < boundary.npoints - 1 && distFwd < searchRadiusPx) {
                distFwd += Math.hypot(
                    boundary.xpoints[endIdx + 1] - boundary.xpoints[endIdx], 
                    boundary.ypoints[endIdx + 1] - boundary.ypoints[endIdx]
                );
                endIdx++;
            }
           
            if (distBack < searchRadiusPx || distFwd < searchRadiusPx) {
                continue;
            }

            int count = endIdx - startIdx + 1;

            double sumX = 0, sumY = 0;
            for (int w = startIdx; w <= endIdx; w++) {
                sumX += boundary.xpoints[w];
                sumY += boundary.ypoints[w];
            }
            double meanX = sumX / count;
            double meanY = sumY / count;

            double Sxx = 0, Syy = 0, Sxy = 0;
            for (int w = startIdx; w <= endIdx; w++) {
                double diffX = boundary.xpoints[w] - meanX;
                double diffY = boundary.ypoints[w] - meanY;
                Sxx += diffX * diffX;
                Syy += diffY * diffY;
                Sxy += diffX * diffY;
            }

            double dx, dy;
            if (Sxx > Syy) {
                double slope = Sxy / Sxx; 
                dx = 1.0;
                dy = slope;
            } else {
                if (Syy == 0) continue;
                double slope = Sxy / Syy; 
                dx = slope;
                dy = 1.0;
            }

            double originalDx = boundary.xpoints[endIdx] - boundary.xpoints[startIdx];
            double originalDy = boundary.ypoints[endIdx] - boundary.ypoints[startIdx];

            if ((dx * originalDx + dy * originalDy) < 0) {
                dx = -dx;
                dy = -dy;
            }
            double length = Math.hypot(dx, dy);
            if (length == 0) continue;

            dx /= length;
            dy /= length;
            
            nx = -dy;
            ny = dx;

            if (invert) { nx = -nx; ny = -ny; }

            double[] profile = new double[profileLength];

            for (int j = -half; j <= half; j++) {
                double x = x0 + j * step * nx;
                double y = y0 + j * step * ny;
                double val = getLanczosInterpolatedPixel(ip, x, y);
                profile[j + half] = cal.getCValue(val);
            }
            
            if (checkContrast) { 
                int avgPts = 3;
                double valBgd = 0;
                double valObject = 0;
                for(int p = 0; p < avgPts; p++) {
                    // profile[] is already calibrated (cal.getCValue applied when sampled);
                    // do NOT calibrate a second time.
                    valBgd += profile[p];
                    valObject += profile[profileLength - 1 - p];
                }
                valBgd /= avgPts;
                valObject /= avgPts;
                
                double edgeContrast = Math.abs(valObject - valBgd);
                
                int bestMatchIdx = -1;
                double minDiff = Double.MAX_VALUE;
                
                for (int c = 0; c < expectedContrasts.size(); c++) {
                    double expected = expectedContrasts.get(c);
                    double diff = Math.abs(edgeContrast - expected);
                    double tol = Math.max(toleranceHU, relTol * Math.abs(expected));
                    if (diff <= tol && diff < minDiff) {
                        minDiff = diff;
                        bestMatchIdx = c;
                    }
                }

                if (bestMatchIdx != -1) {
                    groupedProfiles.get(bestMatchIdx).add(profile);
                }
            } else {
                groupedProfiles.get(0).add(profile);
            }
        }
        
        return groupedProfiles;
    }

    /**
     * Applies the contrast-match gate and files a profile into the right group.
     * Shared by the neighbour-fit and PAE-normal paths so both enforce the same
     * internal-structure rejection (absolute HU floor + relative band).
     */
    private static void addProfileWithContrastCheck(double[] profile, int profileLength,
            boolean checkContrast, List<Double> expectedContrasts,
            double toleranceHU, double relTol, List<List<double[]>> groupedProfiles) {
        if (checkContrast) {
            int avgPts = 3;
            double valBgd = 0, valObject = 0;
            for (int p = 0; p < avgPts; p++) {
                valBgd += profile[p];
                valObject += profile[profileLength - 1 - p];
            }
            valBgd /= avgPts;
            valObject /= avgPts;
            double edgeContrast = Math.abs(valObject - valBgd);

            int bestMatchIdx = -1;
            double minDiff = Double.MAX_VALUE;
            for (int c = 0; c < expectedContrasts.size(); c++) {
                double expected = expectedContrasts.get(c);
                double diff = Math.abs(edgeContrast - expected);
                double tol = Math.max(toleranceHU, relTol * Math.abs(expected));
                if (diff <= tol && diff < minDiff) {
                    minDiff = diff;
                    bestMatchIdx = c;
                }
            }
            if (bestMatchIdx != -1) groupedProfiles.get(bestMatchIdx).add(profile);
        } else {
            groupedProfiles.get(0).add(profile);
        }
    }

    public static double getLanczosInterpolatedPixel(ImageProcessor ip, double x, double y) {
        int a = 3; 
        
        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        
        double sum = 0.0;
        double weightSum = 0.0;
        int width = ip.getWidth();
        int height = ip.getHeight();
        
        for (int i = x0 - a + 1; i <= x0 + a; i++) {
            for (int j = y0 - a + 1; j <= y0 + a; j++) {
                
                double dx = x - i;
                double dy = y - j;
                
                double weight = lanczos(dx, a) * lanczos(dy, a);
                
                int px = Math.max(0, Math.min(i, width - 1));
                int py = Math.max(0, Math.min(j, height - 1));
                
                double val = ip.getf(px, py);
                
                sum += val * weight;
                weightSum += weight;
            }
        }
        
        if (weightSum == 0) return 0;
        return sum / weightSum;
    }

    private static double lanczos(double x, int a) {
        if (x == 0.0) return 1.0;
        if (Math.abs(x) >= a) return 0.0;
        
        double piX = Math.PI * x;
        return (a * Math.sin(piX) * Math.sin(piX / a)) / (piX * piX);
    }
    
    public static double[] calculateTtf(double[] lsf, int fftPoints) {
        int sizeLsf = lsf.length;
        int n = Math.max(fftPoints, sizeLsf);
        
        int fftHalf = n / 2;
        double[] mtfSum = new double[fftHalf + 1];
        DoubleFFT_1D fft = new DoubleFFT_1D(n);
                
        double[] paddedLsf = new double[n];
        int padLeft = (n - sizeLsf) / 2;
        
        System.arraycopy(lsf, 0, paddedLsf, padLeft, sizeLsf);
        
        double[] fftData = new double[2 * n];
        System.arraycopy(paddedLsf, 0, fftData, 0, n);
        fft.realForwardFull(fftData);

        for (int k = 0; k <= fftHalf; k++) {
            double re = fftData[2 * k];
            double im = fftData[2 * k + 1];
            mtfSum[k] += Math.sqrt(re * re + im * im);
        }
        
        double f0 = mtfSum[0];
        if (!(f0 > 0)) {
            f0 = 1.0;
        }
        for (int k = 0; k <= fftHalf; k++) {
            mtfSum[k] /= f0;
        }
        
        return mtfSum;
    }

    // --- Windowing strategy selection (opt-in; default preserves legacy behavior) ---
    // Select with -Dttf.window=current|sidelobe . Tunables (sidelobe) via sysprops below.
    private static String winStrategy() {
        return System.getProperty("ttf.window", "current");
    }
    private static double prop(String key, double def) {
        String v = System.getProperty(key);
        if (v == null) return def;
        try { return Double.parseDouble(v); } catch (NumberFormatException e) { return def; }
    }

    public static double[] derivate(double[] esf) {
        int n = esf.length;
        double[] lsf = new double[n];

        for (int i = 0; i < n - 1; i++) {
            lsf[i] = esf[i + 1] - esf[i];
        }
        lsf[n - 1] = 0;

        if ("sidelobe".equalsIgnoreCase(winStrategy())) {
            return windowSidelobe(lsf);
        }
        return windowCurrent(lsf);
    }

    /**
     * Legacy window: Tukey whose flat-top is derived from the central-peak FWHM
     * (alpha = 1 - 4*FWHM/N). Suppresses noise well but tapers the LSF side-lobes
     * that carry edge-enhancement overshoot (e.g. Lung kernel).
     */
    private static double[] windowCurrent(double[] lsf) {
        int sizeLsf = lsf.length;
        double max = 0;
        for (int k = 0; k < lsf.length; k++) {
            if (Math.abs(lsf[k]) > max) max = Math.abs(lsf[k]);
        }

        int left = 0, right = 0;
        boolean peak = false;
        for (int i = 0; i < sizeLsf; i++) {
            if (Math.abs(lsf[i]) > max / 2 && !peak) { left = i; peak = true; }
            if (Math.abs(lsf[i]) < max / 2 && peak)  { right = i; break; }
        }

        double FWHM = right - left;
        double alpha = 1 - (FWHM * 4 / sizeLsf);
        double steepness = 5;
        return applyTukey(lsf, alpha, steepness);
    }

    /**
     * Side-lobe-aware, noise-adaptive Tukey window.
     *
     * Rationale: an edge-enhancing kernel produces an LSF with significant
     * negative side-lobes whose radius is governed by the kernel (not by the
     * central-peak FWHM). The legacy FWHM-based taper clips those lobes and
     * suppresses the low/mid-frequency MTF overshoot. This strategy instead:
     *   1. Estimates the tail noise sigma from the flat outer region of the LSF.
     *   2. Scans outward from the LSF peak and finds the farthest sample whose
     *      magnitude exceeds k*sigma -> the signal radius R (keeps side-lobes
     *      when they rise above noise; clamps tight on SOFT / high noise).
     *   3. Builds a flat top of half-width R+guard and tapers only beyond it,
     *      with taper steepness scaled by the estimated noise.
     *
     * Tunables (sysprops): ttf.k, ttf.guard, ttf.minAlpha, ttf.steepMin, ttf.steepMax.
     */
    private static double[] windowSidelobe(double[] lsf) {
        int size = lsf.length;
        double k        = prop("ttf.k", 3.0);
        int guard       = (int) prop("ttf.guard", 3);
        double minAlpha = prop("ttf.minAlpha", 0.10);
        double steepMin = prop("ttf.steepMin", 2.0);
        double steepMax = prop("ttf.steepMax", 6.0);

        // 1) peak index and peak magnitude
        int peakIdx = 0;
        double peakMag = 0;
        for (int i = 0; i < size; i++) {
            if (Math.abs(lsf[i]) > peakMag) { peakMag = Math.abs(lsf[i]); peakIdx = i; }
        }

        // 2) tail noise sigma from outer 20% (both ends), far from the edge structure
        int tail = Math.max(3, size / 5);
        double mean = 0; int cnt = 0;
        for (int i = 0; i < tail; i++)            { mean += lsf[i]; cnt++; }
        for (int i = size - tail; i < size; i++)  { mean += lsf[i]; cnt++; }
        mean /= cnt;
        double var = 0;
        for (int i = 0; i < tail; i++)            { double d = lsf[i] - mean; var += d * d; }
        for (int i = size - tail; i < size; i++)  { double d = lsf[i] - mean; var += d * d; }
        double sigma = Math.sqrt(var / cnt);

        // 3) signal radius: farthest sample from the peak above k*sigma
        double thresh = Math.max(k * sigma, peakMag * 1e-6);
        int radius = 0;
        for (int r = 1; r < size; r++) {
            int li = peakIdx - r, ri = peakIdx + r;
            boolean sig = false;
            if (li >= 0   && Math.abs(lsf[li]) > thresh) sig = true;
            if (ri < size && Math.abs(lsf[ri]) > thresh) sig = true;
            if (sig) radius = r;
        }

        // 4) flat-top half-width -> alpha
        int flatHalf = Math.min(radius + guard, size / 2);
        double flatWidth = 2.0 * flatHalf;
        double alpha = 1.0 - flatWidth / size;
        if (alpha < minAlpha) alpha = minAlpha;
        if (alpha > 1.0) alpha = 1.0;

        // 5) steepness scales with relative noise (more noise -> steeper roll-off)
        double snr = peakMag > 0 ? sigma / peakMag : 0.0;          // ~0 low noise .. larger high noise
        double steepness = steepMin + (steepMax - steepMin) * Math.min(1.0, snr / 0.1);

        return applyTukey(lsf, alpha, steepness);
    }

    /** Asymmetric-power Tukey window applied symmetrically over the LSF support. */
    private static double[] applyTukey(double[] lsf, double alpha, double steepness) {
        int sizeLsf = lsf.length;
        double[] windowedLsf = new double[sizeLsf];
        if (alpha <= 0) {
            System.arraycopy(lsf, 0, windowedLsf, 0, sizeLsf);
        } else if (alpha >= 1.0) {
            for (int i = 0; i < sizeLsf; i++) {
                double x = (double) i / (sizeLsf - 1);
                double w = 0.5 * (1 - Math.cos(2 * Math.PI * x));
                windowedLsf[i] = lsf[i] * w;
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
                windowedLsf[i] = lsf[i] * w;
            }
        }
        return windowedLsf;
    }
}

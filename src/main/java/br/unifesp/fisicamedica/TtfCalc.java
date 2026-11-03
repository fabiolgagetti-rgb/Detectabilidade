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

                List<List<double[]>> groupedEsfs = calculateESF(cal, ip, boundary, (int) (10 / cal.pixelWidth), invert, singleExpectedContrastList);
                
                if (groupedEsfs != null && !groupedEsfs.isEmpty() && groupedEsfs.get(0) != null) {
                    esfProfilesAccumulated.addAll(groupedEsfs.get(0));
                }
            }
            
            if (esfProfilesAccumulated.isEmpty()) continue; 
            
            double[] avgEsf = CurveAverage.averageEsfs(esfProfilesAccumulated);
            if (avgEsf == null) continue;

            double[] lsf = derivate(avgEsf);
            double[] mtf = calculateTtf(lsf, 2 * fftPoints);
            
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
        if (boundary == null) throw new IllegalArgumentException("null ROI");
        profileLength *= 2;
        if (profileLength % 2 == 0) profileLength++;

        int half = profileLength / 2;
        double step = 0.5;
        double toleranceHU = 15.0; 
        
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
            
            double nx = -dy;
            double ny = dx;

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
                    valBgd += cal.getCValue(profile[p]);
                    valObject += cal.getCValue(profile[profileLength - 1 - p]);
                }
                valBgd /= avgPts;
                valObject /= avgPts;
                
                double edgeContrast = Math.abs(valObject - valBgd);
                
                int bestMatchIdx = -1;
                double minDiff = Double.MAX_VALUE;
                
                for (int c = 0; c < expectedContrasts.size(); c++) {
                    double expected = expectedContrasts.get(c);
                    double diff = Math.abs(edgeContrast - expected);
                    if (diff <= toleranceHU && diff < minDiff) {
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

    public static double[] derivate(double[] esf) {
        int n = esf.length;
        double[] lsf = new double[n];
        
        for(int i = 0; i < n - 1; i++) {
            lsf[i] = esf[i + 1] - esf[i];
        }
        lsf[n - 1] = 0;

        int sizeLsf = n;
        double max = 0;
        for(int k = 0; k < lsf.length; k++) {
            if(Math.abs(lsf[k]) > max) {
                max = Math.abs(lsf[k]);
            }
        }
        
        int left = 0, right = 0;
        boolean peak = false;
        
        for(int i = 0; i < sizeLsf; i++) {
            if(Math.abs(lsf[i]) > max / 2 && !peak) {
                left = i; peak = true;
            }
            if(Math.abs(lsf[i]) < max / 2 && peak) {
                right = i; break;
            }
        }
        
        double FWHM = right - left;
        double alpha = 1 - (FWHM * 4 / sizeLsf);
        double steepness = 5;

        double[] windowedLsf = new double[sizeLsf];
        if (alpha <= 0) {
            for (int i = 0; i < sizeLsf; i++) windowedLsf[i] = lsf[i];
        } else if (alpha >= 1.0) {
            for (int i = 0; i < sizeLsf; i++) {
                double x = (double)i / (sizeLsf - 1);
                double w = 0.5 * (1 - Math.cos(2 * Math.PI * x));
                windowedLsf[i] = lsf[i] * w;
            }
        } else {
            for (int i = 0; i < sizeLsf; i++) {
                double x = (double) i / (sizeLsf - 1);
                double w;
                double halfAlpha = alpha / 2.0;

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
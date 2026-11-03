package br.unifesp.fisicamedica;

import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.jtransforms.fft.DoubleFFT_2D;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CrossCorrelation {

    public static class CCResult {
        public Point2D.Double center;
        public double contrast;
        
        public CCResult(Point2D.Double center, double contrast) {
            this.center = center;
            this.contrast = contrast;
        }
    }

    public static CCResult findCenter_CC(ImageProcessor imOrig, int scale, double rad, ByteProcessor mask, double tForeground, double tBackground) {
        int w = imOrig.getWidth();
        int h = imOrig.getHeight();
        
        double cx = w / 2.0;
        double cy = h / 2.0;
        
        List<Double> rValsList = new ArrayList<>();
        List<Double> pixelValsList = new ArrayList<>();
        
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (mask.getPixel(x, y) > 0) {
                    // REMOVIDO o +1 que causava o deslocamento para cima e para a esquerda
                    double valX = x - cx; 
                    double valY = y - cy; 
                    double r = Math.hypot(valX, valY);
                    
                    rValsList.add(r);
                    pixelValsList.add((double) imOrig.getf(x, y));
                }
            }
        }
        
        double[] rArray = rValsList.stream().mapToDouble(Double::doubleValue).toArray();
        double[] rSorted = rArray.clone();
        Arrays.sort(rSorted);
        
        double qBkg = getQuantile(rSorted, tBackground);
        double qObj = getQuantile(rSorted, tForeground);
        
        double sumBkg = 0, sumObj = 0;
        int countBkg = 0, countObj = 0;
        
        for (int i = 0; i < rArray.length; i++) {
            if (rArray[i] > qBkg) {
                sumBkg += pixelValsList.get(i);
                countBkg++;
            }
            if (rArray[i] < qObj) {
                sumObj += pixelValsList.get(i);
                countObj++;
            }
        }
        
        double bkg = countBkg > 0 ? sumBkg / countBkg : 0;
        double obj = countObj > 0 ? sumObj / countObj : 0;
        double contrast = obj - bkg;
        
        FloatProcessor imProc = (FloatProcessor) imOrig.convertToFloatProcessor();
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (mask.getPixel(x, y) == 0) {
                    imProc.setf(x, y, 0);
                } else {
                    double val = imProc.getf(x, y) - bkg;
                    if (contrast < 0) {
                        val = val * -1;
                    }
                    imProc.setf(x, y, (float) val);
                }
            }
        }
        
        int newW = w * scale;
        int newH = h * scale;
        imProc.setInterpolationMethod(ImageProcessor.BILINEAR);
        FloatProcessor scaledIm = (FloatProcessor) imProc.resize(newW, newH);
        
        double[] F = new double[newW * newH * 2];
        double[] H = new double[newW * newH * 2];
        
        double sumIm = 0;
        for (int i = 0; i < newW * newH; i++) sumIm += scaledIm.getf(i);
        double MV = sumIm / (newW * newH);
        
        double sumCircle = 0;
        
        for (int y = 0; y < newH; y++) {
            for (int x = 0; x < newW; x++) {
                int idx = 2 * (y * newW + x);
                
                double interpX = (x / (double) scale) - cx; 
                double interpY = (y / (double) scale) - cy; 
                double rInterp = Math.hypot(interpX, interpY);
                
                double circleVal = (rInterp <= rad) ? 1.0 : 0.0;
                sumCircle += circleVal;
                
                F[idx] = scaledIm.getf(x, y) - MV; 
                H[idx] = circleVal;                
            }
        }
        
        if (sumCircle > 0) {
            for (int i = 0; i < H.length; i += 2) {
                H[i] /= sumCircle;
            }
        }
        
        DoubleFFT_2D fft = new DoubleFFT_2D(newH, newW);
        fft.complexForward(F);
        fft.complexForward(H);
        
        F[0] = 0; 
        F[1] = 0; 
        
        for (int i = 0; i < F.length; i += 2) {
            double r1 = F[i], i1 = F[i + 1];
            double r2 = H[i], i2 = H[i + 1];
            F[i] = r1 * r2 - i1 * i2;
            F[i + 1] = r1 * i2 + i1 * r2;
        }
        
        fft.complexInverse(F, true);
        
        double maxPeak = -Double.MAX_VALUE;
        int maxR = 0, maxC = 0;
        
        for (int r = 0; r < newH; r++) {
            for (int c = 0; c < newW; c++) {
                int shiftedR = (r + newH / 2) % newH;
                int shiftedC = (c + newW / 2) % newW;
                
                double val = F[2 * (r * newW + c)] + MV;
                
                if (val > maxPeak) {
                    maxPeak = val;
                    maxR = shiftedR;
                    maxC = shiftedC;
                }
            }
        }
        
        double centerX = (maxC / (double) scale);
        double centerY = (maxR / (double) scale);
        
        return new CCResult(new Point2D.Double(centerX, centerY), contrast);
    }

    private static double getQuantile(double[] sortedData, double p) {
        if (sortedData == null || sortedData.length == 0) return 0;
        if (p <= 0) return sortedData[0];
        if (p >= 1) return sortedData[sortedData.length - 1];
        
        double pos = (sortedData.length - 1) * p;
        int i = (int) pos;
        double frac = pos - i;
        
        return sortedData[i] + frac * (sortedData[i + 1] - sortedData[i]);
    }
}
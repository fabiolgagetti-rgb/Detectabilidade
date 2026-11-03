package br.unifesp.fisicamedica;

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.jtransforms.fft.DoubleFFT_2D;

import ij.ImagePlus;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class NPScalc {

    public static double[][] calculate2DNPS(ImageProcessor ip, Calibration cal, int npsSize) {
        double mean = 0;
        ImageStatistics stats = ip.getStatistics();
        mean = stats.mean;

        int nx = ip.getWidth();
        int ny = ip.getHeight();

        double dx = cal.pixelWidth;
        double dy = cal.pixelHeight;

        double[][] temp = new double[npsSize][npsSize * 2];

        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                temp[y][x] = ip.getPixelValue(x, y) - mean;
            }
        }

        DoubleFFT_2D fft = new DoubleFFT_2D(npsSize, npsSize);
        fft.realForwardFull(temp);

        double[][] powerSpectrum = new double[npsSize][npsSize];
        double normFactor = (dx * dy) / (nx * ny);

        for (int y = 0; y < npsSize; y++) {
            for (int x = 0; x < npsSize; x++) {
                double re = temp[y][2 * x];
                double im = temp[y][2 * x + 1];
                powerSpectrum[y][x] = (re * re + im * im) * normFactor;
            }
        }

        double[][] shifted = fftshift(powerSpectrum);
        shifted[npsSize / 2][npsSize / 2] = 0.0;

        return shifted;
    }

    public static double getIntegratedNPS(double[][] NPS2D) {
        double integratedNPS = 0;
        for (int y = 0; y < NPS2D.length; y++) {
            for (int x = 0; x < NPS2D[0].length; x++) {
                integratedNPS += NPS2D[y][x];
            }
        }
        return integratedNPS;
    }

    public static double[] radialProfile(double[][] powerSpectrum, Calibration cal, double[] radialFreq) {
        int w = powerSpectrum[0].length;
        int h = powerSpectrum.length;

        double pixelSizeMM = cal.pixelWidth;
        double freqStepX = 1.0 / (w * pixelSizeMM);
        double freqStepY = 1.0 / (h * pixelSizeMM);

        int centerX = w / 2;
        int centerY = h / 2;

        int maxBins = (int) Math.ceil((w / 2.0) * Math.sqrt(2));
        double[] tempRadialNPS = new double[maxBins];
        double[] tempRadialFreq = new double[maxBins];
        int[] count = new int[maxBins];

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double fx = (x - centerX) * freqStepX;
                double fy = (y - centerY) * freqStepY;
                double fr = Math.sqrt(fx * fx + fy * fy);

                int bin = (int) Math.round(fr / freqStepX);

                if (bin < maxBins) {
                    tempRadialFreq[bin] += fr;
                    tempRadialNPS[bin] += powerSpectrum[y][x];
                    count[bin]++;
                }
            }
        }

        for (int i = 0; i < maxBins; i++) {
            if (count[i] > 0) {
                tempRadialFreq[i] /= count[i];
                tempRadialNPS[i] /= count[i];
            }
        }

        int outSize = radialFreq.length; 
        double[] finalRadialNPS = new double[outSize];
        
        for (int i = 0; i < outSize; i++) {
            radialFreq[i] = tempRadialFreq[i];
            finalRadialNPS[i] = tempRadialNPS[i];
        }

        return finalRadialNPS;
    }

    public static double[][] fftshift(double[][] input) {
        int h = input.length;
        int w = input[0].length;
        double[][] output = new double[h][w];

        int halfH = (int) Math.ceil(h / 2.0);
        int halfW = (int) Math.ceil(w / 2.0);

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int shiftedY = (y < halfH) ? (y + h / 2) : (y - halfH);
                int shiftedX = (x < halfW) ? (x + w / 2) : (x - halfW);
                output[shiftedY][shiftedX] = input[y][x];
            }
        }
        return output;
    }

    public static void showPowerSpectrum(double[][] spectrum, String title) {
        int height = spectrum.length;
        int width = spectrum[0].length;

        float[] pixels = new float[width * height];

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                pixels[y * width + x] = (float) spectrum[y][x];
            }
        }

        FloatProcessor fp = new FloatProcessor(width, height, pixels);
        ImagePlus imp = new ImagePlus(title, fp);
        imp.show();
    }

    public static double estimateCenterPixel(double[][] spectrum) {
        int centerY = spectrum.length / 2;
        int centerX = spectrum[0].length / 2;
        double media = 0;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                if ((i == 0 && j == 0) || Math.abs(i) + Math.abs(j) == 2) continue;
                media += spectrum[centerY + i][centerX + j];
            }
        }
        return media / 4.0;
    }

    public static double findPeakFreq(double[] freqRadial, double[] npsRadial1D, int larguraJanela) {
        int indicePicoAproximado = -1;
        double maxNPS = 0;
        for (int i = 1; i < npsRadial1D.length; i++) {
            if (npsRadial1D[i] > maxNPS) {
                maxNPS = npsRadial1D[i];
                indicePicoAproximado = i;
            }
        }

        if (indicePicoAproximado == -1) {
            return -1;
        }

        int inicioJanela = Math.max(1, indicePicoAproximado - larguraJanela);
        int fimJanela = Math.min(npsRadial1D.length - 1, indicePicoAproximado + larguraJanela);

        if ((fimJanela - inicioJanela + 1) < 3) {
            return freqRadial[indicePicoAproximado];
        }

        WeightedObservedPoints obs = new WeightedObservedPoints();
        for (int i = inicioJanela; i <= fimJanela; i++) {
            obs.add(freqRadial[i], npsRadial1D[i]);
        }

        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2);
        double[] coeff = fitter.fit(obs.toList());

        double a = coeff[2];
        double b = coeff[1];

        double frequenciaDePico = -b / (2 * a);

        if (a >= 0 || frequenciaDePico < freqRadial[inicioJanela] || frequenciaDePico > freqRadial[fimJanela]) {
            return freqRadial[indicePicoAproximado];
        }

        return frequenciaDePico;
    }
}
package br.unifesp.fisicamedica;

import org.jtransforms.fft.DoubleFFT_2D;

import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class MTFModel{
	

    public static double getMTFValue(double frequencyLpMm, double[][] gaussians) {
    	
        double u = frequencyLpMm * 10.0; 
        
        double mtf = 0.0;
        
        for (double[] parameter : gaussians) {
            mtf += parameter[0] * Math.exp(-Math.pow((u - parameter[1]) / parameter[2],2));
        }
        
        return mtf;
    }
    
    public double[][] gerarMTF2D(int width, int height, Calibration cal, double[][] gaussians) {
        double[][] mtf2D = new double[height][width];
        
        double pixelSizeMm = cal.pixelWidth; 
        double freqStepX = 1.0 / (width * pixelSizeMm);
        double freqStepY = 1.0 / (height * pixelSizeMm);
        
        int centerX = width / 2;
        int centerY = height / 2;

        double dcValue = getMTFValue(0.0, gaussians); 

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double fx = (x - centerX) * freqStepX;
                double fy = (y - centerY) * freqStepY;
                
                double fr = Math.sqrt(fx * fx + fy * fy);
                mtf2D[y][x] = getMTFValue(fr, gaussians) / dcValue; 
            }
        }
        return mtf2D;
    }
    
    public static ImageProcessor applyMTF(ImageProcessor originalImage, double[][] mtf2D) {
        int height = originalImage.getHeight();
        int width = originalImage.getWidth();

        ImageProcessor ipFloat = originalImage.convertToFloatProcessor();
        float[] pixels = (float[]) ipFloat.getPixels();
        double[][] fftData = new double[height][width * 2];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                fftData[y][2 * x] = pixels[y * width + x];
            }
        }
        
        DoubleFFT_2D fft = new DoubleFFT_2D(height, width);
        fft.complexForward(fftData);
        
        double[][] mtfShifted = NPScalc.fftshift(mtf2D);
        

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double factor = mtfShifted[y][x];
                fftData[y][2 * x] *= factor;
                fftData[y][2 * x + 1] *= factor;
            }
        }

        fft.complexInverse(fftData, true);

        float[] finalImage = new float[width * height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                finalImage[y * width + x] = (float) fftData[y][2 * x];
            }
        }
        return new FloatProcessor(width, height, finalImage);
    }
    
    /**
     * Retorna a matriz de parametros para a MTF obtidos empíricamente por Bujila et al. (2017)
     * * @param kernelName O nome do kernel
     * @return double[][] armazenando parâmetros de até 4 termos gaussianos para modelo de MTF
     */
    public static double[][] getMTFparameters(String kernelName) {
        
        switch (kernelName.toUpperCase()) {

            case "SOFT":
                return new double[][] {
                    {1.2442e+03, -2.8546e+01, 9.8841e+00},
                    {7.7973e-01,  1.0908e+00, 3.5111e+00}
                };

            case "STANDARD":
                return new double[][] {
                    {9.3885e-01, -4.2840e-01, 3.2308e+00},
                    {1.2519e-01,  2.1679e+00, 2.0414e+00},
                    {3.0915e-01,  3.8610e+00, 2.6403e+00}
                };

            case "DETAIL":
                return new double[][] {
                    {6.9485e-01, -1.4513e+00, 3.4939e+00},
                    {2.6816e-02,  2.6578e+00, 1.9685e+00}
                };

            case "CHEST":
                return new double[][] {
                    {3.8010e+03, -2.8614e+00, 7.8577e-01},
                    {7.3351e-02, -2.5815e-01, 1.8215e+00}
                };

            case "LUNG":
                return new double[][] {
                    {6.7589e-02,  2.1542e+00, 1.3908e-02},
                    {1.1687e+00,  3.7561e+00, 2.6115e+00},
                    {6.5667e-01,  2.8016e+00, 4.0885e+00},
                    {9.1152e-01, -1.1016e+00, 4.1051e+00}
                };

            case "BONE":
                return new double[][] {
                    { 1.0810e-01, 8.1258e+00, 2.5573e+00},
                    {-9.3344e-01, 1.3119e+00, 3.8808e+00},
                    { 9.2722e-01, 3.8673e-01, 4.9096e+00}
                };

            case "BONEPLUS":
                return new double[][] {
                    {-9.8745e-02, 1.2421e+00, 2.1032e+00},
                    {-7.2080e+00, 1.0646e+00, 5.9054e+00},
                    { 3.1412e-01, 6.2049e+00, 1.8031e+00},
                    { 8.3188e+00, 1.1231e+00, 6.1514e+00}
                };

            case "EDGE":
                return new double[][] {
                    { 9.8810e-01, -1.7257e-01, 6.1247e+00},
                    {-2.0130e-01,  7.4007e+00, 3.3801e+00},
                    { 1.9079e+00,  1.0022e+00, 4.9966e+00}
                };
        }
		return null;
    }
}
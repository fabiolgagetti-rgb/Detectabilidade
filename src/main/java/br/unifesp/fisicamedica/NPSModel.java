package br.unifesp.fisicamedica;

import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.jtransforms.fft.DoubleFFT_2D;
import java.util.Random;

public class NPSModel {

    /**
     * Calcula o valor da NPS para uma dada frequência.
     * * @param frequency Frequência radial em lp/mm.
     * @param params Array de parametros do modelo na ordem: [a, b, c, g, h, k].
     * @param variance A variância desejada para a imagem final.
     * @return O valor da NPS absoluta na frequência dada.
     */
    public static double getNPSValue(double frequency, double[] params, double variance) {
        
        // Converte a frequência de lp/mm para lp/cm para o cálculo com parametros empiricos de Bujila et al. (2017)
        double u = frequency * 10.0;
        
        double a = params[0];
        double b = params[1];
        double c = params[2];
        double g = params[3];
        double h = params[4];
        double k = params[5];
        
        // u=0 para evitar return NaN se b for fracionário/negativo
        if (u == 0.0) return 0.0;

        double termExp = Math.exp(-c * Math.pow(Math.abs(Math.pow(u, g) - h), k));
        double nnpsCm2 = a * Math.pow(u, b) * termExp;
        
        // O valor de nnpsCm2, com os parâmetros empiricos dados por Bujila et al. tem unidade de cm²
        // Para trazer para unidade de mm², multiplicamos por 100
        double nnpsMm2 = nnpsCm2 * 100.0;
        
        // Multiplica a NPS normalizada pela variância alvo para obter a NPS alvo
        return variance * nnpsMm2;
    }
    
    /**
     * Gera a matriz 2D da NPS radial.
     */
    public double[][] gerarNPS2D(int width, int height, Calibration cal, double[] params, double variance) {
        double[][] nps2D = new double[height][width];
        
        double pixelSizeMm = cal.pixelWidth; 
        double freqStepX = 1.0 / (width * pixelSizeMm);
        double freqStepY = 1.0 / (height * pixelSizeMm);
        
        int centerX = width / 2;
        int centerY = height / 2;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double fx = (x - centerX) * freqStepX;
                double fy = (y - centerY) * freqStepY;
                
                // fr calculado em lp/mm
                double fr = Math.sqrt(fx * fx + fy * fy);
                nps2D[y][x] = getNPSValue(fr, params, variance);
            }
        }
        return nps2D;
    }

    /**
     * Gera uma imagem de ruído com as características espectrais da NPS alvo.
     */
    public ImageProcessor texturedNoiseIP(double[][] nps2D, double pixelSize) {
        int height = nps2D.length;
        int width = nps2D[0].length;

        Random random = new Random();
        
        // Prepara dados para FFT Complexa
        double[][] fftData = new double[height][width * 2];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // Parte real com ruído branco/gaussiano
                fftData[y][2 * x] = random.nextGaussian(); 
                // Parte imaginária = 0
                fftData[y][2 * x + 1] = 0.0; 
            }
        }
        
        DoubleFFT_2D fft = new DoubleFFT_2D(height, width);
        
        //Transforma o ruído branco para o domínio da frequência
        fft.complexForward(fftData);

        //Prepara o filtro do ruido (Raiz quadrada da NPS) para modular o ruido aleatório com as magnitudes da NPS na frequencia
        double[][] npsSqrt = new double[height][width];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                npsSqrt[y][x] = nps2D[y][x] > 0 ? Math.sqrt(nps2D[y][x]) : 0.0;
            }
        }
        
        // Desloca freq=0 do centro da imagem para o canto, alinhando com a saída da FFT
        double[][] npsSqrtShifted = NPScalc.fftshift(npsSqrt);
        
        // Filtragem no domínio da frequência, criando a correlação/textura desejada do ruido
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double factor = npsSqrtShifted[y][x] / pixelSize;
                fftData[y][2 * x] *= factor;
                fftData[y][2 * x + 1] *= factor;
            }
        }
        
        fft.complexInverse(fftData, true);

        //Extrai a parte real para os pixels da imagem do ruído texturizado
        float[] pixelsFinais = new float[width * height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                pixelsFinais[y * width + x] = (float) fftData[y][2 * x];
            }
        }
        
        return new FloatProcessor(width, height, pixelsFinais);
    }
    
    /**
     * Retorna os parâmetros [a, b, c, g, h, k] para a NNPS 
     * com base no nome do kernel de reconstrução.
     * @param kernelName O nome do kernel
     * @return double[] com os parametros obtidos empiricamente por Bujila et al. (2017)
     */
    public static double[] getNPSCoefficients(String kernelName) {
        
        switch (kernelName.toUpperCase()) {
            
        case "SOFT":
                return new double[] {
                    8.9886e-03, 6.6698e-01, 1.4968e-02, 1.5603e+00, 3.8641e-01, 1.8771e+00
                };

            case "STANDARD":
                return new double[] {
                    6.0671e-03, 7.2201e-01, 6.6041e-24, 2.3026e-01, -2.6551e+00, 3.7845e+01
                };

            case "DETAIL":
                return new double[] {
                    3.8133e-03, 6.8606e-01, 1.6902e-02, 1.2087e+00, 7.2924e-01, 2.0170e+00
                };

            case "CHEST":
                return new double[] {
                    2.2764e-03, 7.3836e-01, 2.5600e-03, 2.2723e+00, 1.1036e+01, 1.4003e+00
                };

            case "LUNG":
                return new double[] {
                    5.7602e-04, 1.7366e+00, 1.2895e-03, 2.0725e-01, 1.4757e-01, 2.3349e+01
                };

            case "BONE":
                return new double[] {
                    8.7379e-04, 6.5581e-01, 2.6195e-06, 2.6514e+00, 1.8753e+02, 2.1996e+00
                };

            case "BONEPLUS":
                return new double[] {
                    2.1297e-04, 1.4932e+00, 2.9477e-03, 6.9186e-01, 2.2119e+00, 6.0085e+00
                };

            case "EDGE":
                return new double[] {
                    2.0572e-04, 1.1053e+00, 1.9844e-09, 4.4849e+00, 3.6214e+04, 1.8237e+00
                };
        }
		return null;
    }
}
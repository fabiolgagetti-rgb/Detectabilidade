package br.unifesp.fisicamedica;

import org.apache.commons.math3.special.Erf;

public class TheoreticalDetectabilityTest {

    /**
     * Método principal de teste. Pode ser chamado a partir de um botão na interface
     * ou rodado isoladamente no modo de desenvolvimento.
     */
    public static void runTest() {
        String kernel = "LUNG";          // Tente "LUNG", "STANDARD", etc.
        int numPoints = 33;              // Equivale aos 33 bins empíricos
        double df = 0.032;               // Exatamente 1 / (64 * tamanho_do_pixel)
        double variance = 100.0;         // Variância do ruído (SD = 10 HU -> Var = 100)

        // Parâmetros da Tarefa e do Observador
        double taskRadius = 1.5;         // mm (Raio do insert)
        double taskContrast = 30.0;      // HU
        double viewDistance = 500.0;     // mm
        double displaySize = 250.0;      // mm
        double FOV = 250.0;              // mm

        double[][] mtfParams = MTFModel.getMTFparameters(kernel);
        double[] npsParams = NPSModel.getNPSCoefficients(kernel);

        double[] freqArray = new double[numPoints];
        double[] mtfData = new double[numPoints];
        double[] npsData = new double[numPoints];

        for (int i = 0; i < numPoints; i++) {
            double f = i * df;
            freqArray[i] = f;
            
            mtfData[i] = MTFModel.getMTFValue(f, mtfParams);
            npsData[i] = NPSModel.getNPSValue(f, npsParams, variance);
        }

        DetectabilityIndexCalculator dPrimeCalc = new DetectabilityIndexCalculator(
                mtfData, freqArray, 
                npsData, freqArray, 
                taskRadius, taskContrast, viewDistance, displaySize, FOV
        );

        double dPrime = dPrimeCalc.compute();
        double auc = 0.5 + 0.5 * Erf.erf(dPrime / 2.0);


        System.out.println("Kernel        : " + kernel);
        System.out.println("Bins Freq.    : " + numPoints + " pontos");
        System.out.println("Resolução (df): " + df + " cy/mm");
        System.out.println("Nyquist Máx   : " + freqArray[numPoints - 1] + " cy/mm");
        System.out.println("Tarefa        : " + taskContrast + " HU | Raio: " + taskRadius + " mm");
        System.out.println("--------------------------------------------------");
        System.out.printf("d' Calculado  : %.6f\n", dPrime);
        System.out.printf("auc Calculado  : %.6f\n", auc);

        System.out.println("==================================================");
    }

    public static void main(String[] args) {
        runTest();
    }
}
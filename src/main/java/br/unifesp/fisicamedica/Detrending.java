package br.unifesp.fisicamedica;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;

/**
 * Classe utilitária para realizar o "detrending" de dados 2D (ROIs).
 * Remove tendências de baixa frequência ajustando uma superfície polinomial
 * e subtraindo-a dos dados originais.
 */
public class Detrending {

    /**
     * Realiza o detrending de uma ROI 2D ajustando um polinômio de 2ª ordem
     * e subtraindo a superfície ajustada.
     * A equação do polinômio é: f(x, y) = c0 + c1*x + c2*y + c3*x^2 + c4*y^2 + c5*xy
     *
     * @param roiData A ROI 2D original, como uma matriz de doubles.
     * @return Uma nova matriz 2D contendo os dados da ROI após o detrending.
     */
    public static double[][] detrendWith2ndOrderPolynomial(float[][] roiData) {
        int rows = roiData.length;
        if (rows == 0) {
            return new double[0][0];
        }
        int cols = roiData[0].length;
        int numPoints = rows * cols;
        int numCoefficients = 6;

        if (numPoints < numCoefficients) {
            throw new IllegalArgumentException("A ROI tem menos pontos de dados (" + numPoints + ") do que coeficientes a serem ajustados (" + numCoefficients + ").");
        }

        RealMatrix A = new Array2DRowRealMatrix(numPoints, numCoefficients);

        RealVector b = new ArrayRealVector(numPoints);

        int k = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double x = j;
                double y = i;

                A.setEntry(k, 0, 1.0);
                A.setEntry(k, 1, x);
                A.setEntry(k, 2, y);
                A.setEntry(k, 3, x * x);
                A.setEntry(k, 4, y * y);
                A.setEntry(k, 5, x * y);

                b.setEntry(k, roiData[i][j]);
                k++;
            }
        }

        DecompositionSolver solver = new QRDecomposition(A).getSolver();
        RealVector coefficients = solver.solve(b);

        double[][] detrendedData = new double[rows][cols];
        double c0 = coefficients.getEntry(0);
        double c1 = coefficients.getEntry(1);
        double c2 = coefficients.getEntry(2);
        double c3 = coefficients.getEntry(3);
        double c4 = coefficients.getEntry(4);
        double c5 = coefficients.getEntry(5);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double x = j;
                double y = i;

                double fittedValue = c0 + c1*x + c2*y + c3*x*x + c4*y*y + c5*x*y;
                
                detrendedData[i][j] = roiData[i][j] - fittedValue;
            }
        }

        return detrendedData;
    }
}
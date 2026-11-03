package br.unifesp.fisicamedica;

public class DetectabilityIndexCalculator {

    private double[] ttfData;
    private double[] ttfFreq;
    private double[] npsData;
    private double[] npsFreq;
    private double taskRadius;
    private double taskContrast;
    private double viewDistance1;
    private double displaySize1;
    private double FOV1;

        
    public DetectabilityIndexCalculator(double[] mtfData, double[] mtfFreq, double[] npsData, double[] npsFreq,
                                double taskRadius, double taskContrast, double viewDistance, double displaySize, double FOV) {
        this.ttfData = mtfData;
        this.ttfFreq = mtfFreq;
        this.npsData = npsData;
        this.npsFreq = npsFreq;
        this.taskRadius = taskRadius;
        this.taskContrast = taskContrast;
        this.viewDistance1 = viewDistance;
        this.displaySize1 = displaySize;
        this.FOV1=FOV;
    }

    /**
     * Interpola linearmente o valor da TTF para uma frequência alvo da NPS.
     * @param targetFreq A frequência da malha da NPS.
     * @return O valor da TTF interpolado (ou 0 se exceder o limite medido).
     */
    private double interpolateMTF(double targetFreq) {
        if (targetFreq >= this.ttfFreq[this.ttfFreq.length - 1]) {
            return 0.0;
        }
        if (targetFreq <= this.ttfFreq[0]) {
            return this.ttfData[0];
        }

        for (int j = 0; j < this.ttfFreq.length - 1; j++) {
            if (targetFreq >= this.ttfFreq[j] && targetFreq <= this.ttfFreq[j + 1]) {
                double f0 = this.ttfFreq[j];
                double f1 = this.ttfFreq[j + 1];
                double mtf0 = this.ttfData[j];
                double mtf1 = this.ttfData[j + 1];

                return mtf0 + (targetFreq - f0) * (mtf1 - mtf0) / (f1 - f0);
            }
        }
        return 0.0;
    }

    /**
     * Executa o cálculo completo do d'.
     * Utiliza integração numérica radial pela Regra do Trapézio.
     * @return O valor de d'.
     */
    /**
     * Executa o cálculo completo do d'.
     * Utiliza integração numérica radial (Regra do Trapézio para malhas não-uniformes).
     * @return O valor de d'.
     */
    public double compute() {
        double[] sSquared = calculateTaskFunctionSquared(this.npsFreq, this.taskRadius, this.taskContrast);
        double[] eSquared = calculateEyeModelSquared(this.npsFreq, viewDistance1, displaySize1, FOV1);

        double integralNumeratorSum = 0.0;
        double integralDenominatorSum = 0.0;

        double prev_r = this.npsFreq[0];
        double ttf_val_0 = interpolateMTF(prev_r);
        double ttf2_0 = ttf_val_0 * ttf_val_0;
        double e4_0 = eSquared[0] * eSquared[0];
        
        double prevNum = sSquared[0] * ttf2_0 * eSquared[0] * prev_r;
        double prevDen = sSquared[0] * ttf2_0 * this.npsData[0] * e4_0 * prev_r;

        for (int i = 1; i < this.npsFreq.length; i++) {
            double r = this.npsFreq[i];
            
            double ttf_val = interpolateMTF(r);
            double ttf2 = ttf_val * ttf_val;
            double e4 = eSquared[i] * eSquared[i];
            
            double currentNum = sSquared[i] * ttf2 * eSquared[i] * r;
            double currentDen = sSquared[i] * ttf2 * this.npsData[i] * e4 * r;
            
            double dr = r - prev_r;

            integralNumeratorSum += ((prevNum + currentNum) / 2.0) * dr;
            integralDenominatorSum += ((prevDen + currentDen) / 2.0) * dr;

            prev_r = r;
            prevNum = currentNum;
            prevDen = currentDen;
        }

        double numerator = integralNumeratorSum;
        double denominator = Math.sqrt(integralDenominatorSum);

        if (denominator == 0.0) {
            return Double.POSITIVE_INFINITY;
        }
        
        return Math.sqrt(2 * Math.PI) * numerator / denominator;
    }

    private static double[] calculateTaskFunctionSquared(double[] rho, double radius, double contrast) {
        double[] s2 = new double[rho.length];
        org.apache.commons.math3.special.BesselJ besselJ1 = new org.apache.commons.math3.special.BesselJ(1);

        for (int i = 0; i < rho.length; i++) {
            double s_val;
            if (rho[i] == 0.0) {
                s_val = contrast * (Math.PI * radius * radius);
            } else {
                double argument = 2 * Math.PI * radius * rho[i];
                s_val = contrast * (radius / rho[i]) * besselJ1.value(argument);
            }
            s2[i] = s_val * s_val;
        }
        return s2;
    }

    private static double[] calculateEyeModelSquared(double[] rho, double viewingDistance, double displaySize, double fov) {
        double[] e2 = new double[rho.length];
        
        double a1 = 1.5;
        double a2 = 0.98;
        double a3 = 0.68;
        
        double peak_p = Math.pow(a1 / (a2 * a3), 1.0 / a3); 
        double max_raw_val = Math.pow(peak_p, a1) * Math.exp(-a2 * Math.pow(peak_p, a3));
        
        double n = 1.0 / max_raw_val;

        for (int i = 0; i < rho.length; i++) {
            double f = rho[i];
            
            if (f == 0.0) {
                e2[i] = 0.0;
                continue;
            }
            double p = (f * fov * viewingDistance * Math.PI) / (displaySize * 180.0);
            
            double e_raw = Math.pow(p, a1) * Math.exp(-a2 * Math.pow(p, a3));
            
            double e_norm = n * e_raw;
            
            e2[i] = e_norm*e_norm; 
        }
        
        return e2;
    }
}
package br.unifesp.fisicamedica;

import ij.ImagePlus;
import ij.measure.Calibration;

public class Tester {

    public static void main(String[] args) {
                
        int nSlicesHomogeneo = 10;
        int nSlicesInserts = 10;
        int matriz = 512;
        double fov = 250.0;
        double bgHU = -1000.0;
        double objectMeanHU = 0.0;
        double diametroMM = 200.0;
        
        double[] insertHU = {1000, -100, 50, -30};
        double insertRaioMM = 12.5;
        double insertDistanciaMM = diametroMM / 4.0;
        
        double trendAmpHU = 0.0;
        String kernelName = "SOFT"; 
        double variance = 9.58;

        PhantomImageCreator creator = new PhantomImageCreator();
                
        ImagePlus impFinal = creator.criarImagemCTPhantom(
                nSlicesHomogeneo,
                nSlicesInserts,
                matriz,
                fov,
                bgHU,
                objectMeanHU,
                diametroMM,
                trendAmpHU,
                insertHU,
                insertRaioMM,
                insertDistanciaMM,
                kernelName,
                variance
        );impFinal.show();
        
        String dicomSavePath = "C:/Users/fabio/OneDrive/Desktop/Simulated phantom/teste";
        String dicomTemplatePath = "C:/Users/fabio/OneDrive/Desktop/backups/Simulated phantom/Z9600";
        
        Calibration cal = impFinal.getCalibration();
        
        try {
            DicomSaver.saveAsDicomFromTemplate(impFinal, dicomTemplatePath, dicomSavePath, cal);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
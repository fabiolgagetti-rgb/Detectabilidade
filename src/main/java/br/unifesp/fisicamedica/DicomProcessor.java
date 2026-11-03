package br.unifesp.fisicamedica;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JDialog;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

import ij.IJ;
import ij.io.Opener;
import ij.plugin.DICOM;

public class DicomProcessor {
    public static Map<String, Map<String, Map<String, Map<String, SeriesData>>>> patientData = new HashMap<>();
    int processed = 0;
    long totalFiles=0;

    public DicomProcessor() {
    }

    public void analyzeDirectory(File dir, JDialog loadingDialog) {
        File[] files = dir.listFiles();
        if (files != null) {
        	try {
               if(totalFiles==0) totalFiles = Files.walk(dir.toPath()).filter(Files::isRegularFile).count();
            } catch (IOException e) {
                e.printStackTrace();
                totalFiles=0;
            }
            JProgressBar progressBar = (JProgressBar) loadingDialog.getContentPane().getComponent(1);
            for (File file : files) {
                if (file.isDirectory()) {
                    analyzeDirectory(file, loadingDialog);
                } else {
                    try {
						addFileToPatientData(file);
						processed++;
						int progress = (int) ((double) 100*processed / totalFiles);
						 SwingUtilities.invokeLater(() -> progressBar.setValue(progress));
					} catch (IOException e) {
						IJ.log("Error: " + file.getName() + " - " + e.getMessage());
						e.printStackTrace();
					}
                }
            }
        }
    }

    private void addFileToPatientData(File file) throws IOException {
    	int type = new Opener().getFileType(file.getAbsolutePath());
		if (!(type==1||type==2||type==14)||file.getName().equals("DICOMDIR")) {
			return;
		}
    	DICOM dicom = new DICOM();
    	String	metadata = dicom.getInfo(file.getAbsolutePath());
    	if (metadata==null) {return;}
        String patientName = getTag(metadata, "0010,0010");
        String studyDescription = getTag(metadata, "0008,1030");
        String protocol = getTag(metadata, "0018,1030");
        String seriesDescription = getTag(metadata, "0008,103E");
        String uniqueSeriesKey = getTag(metadata, "0020,000E");
        SeriesData seriesData = patientData
            .computeIfAbsent(patientName, k -> new HashMap<>())
            .computeIfAbsent(studyDescription, k -> new HashMap<>())
            .computeIfAbsent(protocol, k -> new HashMap<>())
            .computeIfAbsent(uniqueSeriesKey, k -> new SeriesData(seriesDescription, uniqueSeriesKey));
        seriesData.addFile(file, metadata);
    }
    
    public static String getTag(String metadata, String tag) {
        if (metadata == null) {
            return null;
        }
        int index1 = metadata.indexOf(tag);
        if (index1 == -1) {
            return null;
        }
        if (index1 + 11 < metadata.length() && metadata.charAt(index1 + 11) == '>') {
            index1 = metadata.indexOf(tag, index1 + 10);
            if (index1 == -1) {
                return null;
            }
        }
   index1 = metadata.indexOf(":", index1);
        if (index1 == -1) {
            return null; 
        }
        int index2 = metadata.indexOf("\n", index1);
        if (index2 == -1) {
            return null;
        }
        if(metadata.substring(index1 + 1, index2)==null) {
        	return null;
        }
        return metadata.substring(index1 + 1, index2).trim();
    }
    public Map<String, Map<String, Map<String, Map<String, SeriesData>>>> getPatientData() {
    	if(patientData==null) {
        } return patientData;
    }
}
package br.unifesp.fisicamedica;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import org.dcm4che3.data.Attributes;
import org.dcm4che3.data.Tag;
import org.dcm4che3.data.UID;
import org.dcm4che3.data.VR;
import org.dcm4che3.io.DicomInputStream;
import org.dcm4che3.io.DicomOutputStream;
import org.dcm4che3.util.UIDUtils;

import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.ImageProcessor;

import java.util.Locale;

public class DicomSaver {

    /**
     * Salva um ImagePlus como uma série DICOM, usando um arquivo DICOM existente como modelo de cabeçalho.
     *
     * @param imp O ImagePlus com o stack de imagens a ser salvo.
     * @param templateDicomPath O caminho para o arquivo DICOM que servirá como modelo.
     * @param outputDirectoryPath O diretório onde a nova série DICOM será salva.
     */
    @SuppressWarnings("deprecation")
	public static void saveAsDicomFromTemplate(ImagePlus imp, String templateDicomPath, String outputDirectoryPath, Calibration cal) {
        File templateFile = new File(templateDicomPath);
        if (!templateFile.exists()) {
            System.err.println("Erro: Arquivo de modelo DICOM não encontrado em: " + templateDicomPath);
            return;
        }

        File outputDir = new File(outputDirectoryPath);
        if (!outputDir.exists()) {
            outputDir.mkdirs();
        }

        try {
            Attributes templateHeader;
            try (DicomInputStream dis = new DicomInputStream(templateFile)) {
                templateHeader = dis.readDataset(-1, -1);
            }

            ImageStack stack = imp.getStack();
            
            String seriesUID = UIDUtils.createUID();
            templateHeader.setString(Tag.SeriesInstanceUID, VR.UI, seriesUID);
            
            LocalDateTime now = LocalDateTime.now();
            String date = now.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
            String time = now.format(DateTimeFormatter.ofPattern("HHmmss.S"));

            for (int i = 1; i <= stack.getSize(); i++) {
                Attributes sliceHeader = new Attributes(templateHeader);
                ImageProcessor ip = stack.getProcessor(i);

                sliceHeader.setString(Tag.SOPInstanceUID, VR.UI, UIDUtils.createUID());
                sliceHeader.setString(Tag.InstanceNumber, VR.IS, Integer.toString(i));
                sliceHeader.setString(Tag.SeriesDate, VR.DA, date);
                sliceHeader.setString(Tag.SeriesTime, VR.TM, time);
                sliceHeader.setString(Tag.ContentDate, VR.DA, date);
                sliceHeader.setString(Tag.ContentTime, VR.TM, time);

                sliceHeader.setInt(Tag.Rows, VR.US, ip.getHeight());
                sliceHeader.setInt(Tag.Columns, VR.US, ip.getWidth());
                
                sliceHeader.setString(Tag.PixelSpacing, VR.DS, String.format(Locale.US, "%f\\%f", cal.pixelWidth, cal.pixelHeight));
                sliceHeader.setString(Tag.SliceThickness, VR.DS, String.format(Locale.US, "%f", cal.pixelDepth));
                double zPosition = (i - 1) * cal.pixelDepth;
                sliceHeader.setString(Tag.ImagePositionPatient, VR.DS, String.format(Locale.US, "0.0\\0.0\\%f", zPosition));
                sliceHeader.setDouble(Tag.RescaleIntercept, VR.DS, cal.getCoefficients()[0]);
                sliceHeader.setDouble(Tag.RescaleSlope, VR.DS, cal.getCoefficients()[1]);
                sliceHeader.setString(Tag.RescaleType, VR.LO, cal.getValueUnit());
                sliceHeader.setString(Tag.Modality, VR.CS, "CT");
                sliceHeader.setString(Tag.SOPClassUID, VR.UI, UID.CTImageStorage);
                sliceHeader.setString(Tag.RescaleType, VR.LO, "HU");
                sliceHeader.setInt(Tag.BitsAllocated, VR.US, 16);
                sliceHeader.setInt(Tag.BitsStored, VR.US, 16);
                sliceHeader.setInt(Tag.HighBit, VR.US, 15);
                sliceHeader.setInt(Tag.PixelRepresentation, VR.US, 0); // 1 = signed, 0 = unsigned

                Object pixelArray = ip.getPixels();
                byte[] pixelData;
                if (pixelArray instanceof short[]) {
                    short[] pixels = (short[]) pixelArray;
                    ByteBuffer byteBuffer = ByteBuffer.allocate(pixels.length * 2);
                    byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
                    for (short pixel : pixels) {
                        byteBuffer.putShort(pixel);
                    }
                    pixelData = byteBuffer.array();
                } else if (pixelArray instanceof float[]) {
                    float[] floatPixels = (float[]) pixelArray;
                    ByteBuffer byteBuffer = ByteBuffer.allocate(floatPixels.length * 2);
                    byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
                    for (float pixel : floatPixels) {
                        byteBuffer.putShort((short) Math.round(pixel));
                    }
                    pixelData = byteBuffer.array();
                } else {
                    throw new IOException("Tipo de pixel não suportado para salvamento DICOM: " + pixelArray.getClass().getName());
                }

                sliceHeader.setBytes(Tag.PixelData, VR.OW, pixelData);
                Attributes fileMetaInformation = new Attributes();
                fileMetaInformation.setString(Tag.MediaStorageSOPClassUID, VR.UI, UID.CTImageStorage);
                fileMetaInformation.setString(Tag.MediaStorageSOPInstanceUID, VR.UI, sliceHeader.getString(Tag.SOPInstanceUID));

             sliceHeader.setString(Tag.PatientName, VR.PN, "PHANTOM^VIRTUAL");
             sliceHeader.setString(Tag.PatientID, VR.LO, ""); 
             sliceHeader.setString(Tag.PatientBirthDate, VR.DA, "");
             sliceHeader.setString(Tag.PatientSex, VR.CS, "O");

             sliceHeader.setString(Tag.InstitutionName, VR.LO, ""); 
             sliceHeader.setString(Tag.InstitutionAddress, VR.ST, "");
             sliceHeader.setString(Tag.ReferringPhysicianName, VR.PN, "");
             sliceHeader.setString(Tag.PerformingPhysicianName, VR.PN, "");
             sliceHeader.setString(Tag.OperatorsName, VR.PN, "");

             sliceHeader.setString(Tag.StudyID, VR.SH, "001");
             sliceHeader.setString(Tag.AccessionNumber, VR.SH, "");
             sliceHeader.setString(Tag.StudyDescription, VR.LO, "");
             sliceHeader.setString(Tag.SeriesDescription, VR.LO, "");

             sliceHeader.setString(Tag.Manufacturer, VR.LO, "");
             sliceHeader.setString(Tag.ManufacturerModelName, VR.LO, "");
             sliceHeader.setString(Tag.StationName, VR.SH, "");
             sliceHeader.setString(Tag.DeviceSerialNumber, VR.LO, "000000");
             sliceHeader.setString(Tag.SoftwareVersions, VR.LO, "v1.0");

             sliceHeader.remove(Tag.OtherPatientIDs);
             sliceHeader.remove(Tag.OtherPatientNames);
                fileMetaInformation.setString(Tag.TransferSyntaxUID, VR.UI, UID.ExplicitVRLittleEndian);
                
                fileMetaInformation.setString(Tag.ImplementationClassUID, VR.UI, UIDUtils.createUID());
                fileMetaInformation.setString(Tag.ImplementationVersionName, VR.LO, "IQ_TOOLKIT_1.0");

                String fileName = String.format("slice-%04d.dcm", i);
                File file = new File(outputDir, fileName);
                
                try (DicomOutputStream dos = new DicomOutputStream(file)) {
                    dos.writeDataset(fileMetaInformation, sliceHeader);
                }
            }
            
            System.out.println("Série DICOM salva em: " + outputDir.getAbsolutePath());

        } catch (IOException e) {
            System.err.println("Erro ao salvar a série DICOM.");
            e.printStackTrace();
        }
    }
}
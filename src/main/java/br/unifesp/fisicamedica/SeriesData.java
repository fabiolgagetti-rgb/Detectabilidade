package br.unifesp.fisicamedica;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class SeriesData extends DicomProcessor{
    private String seriesDescription;
    private String uniqueSeriesKey;
    private Map<File, String> fileMetadata;

    public SeriesData(String seriesDescription, String uniqueSeriesKey) {
        this.seriesDescription = seriesDescription;
        this.uniqueSeriesKey = uniqueSeriesKey;
        this.fileMetadata = new LinkedHashMap<>();
    }

    public String getSeriesDescription() {
        return seriesDescription;
    }

    public String getUniqueSeriesKey() {
        return uniqueSeriesKey;
    }

    public List<File> getFiles() {
        return new ArrayList<>(fileMetadata.keySet());
    }

    public void addFile(File file, String metadata) {
        fileMetadata.put(file, metadata);
    }
    
    public String getMetadata(File file) {
        return fileMetadata.get(file);
    }
    
    public void sortFileMetadataByInstanceNumber() {
        fileMetadata = fileMetadata.entrySet().stream()
                .sorted((entry1, entry2) -> {
                    int instanceNumber1 = extractInstanceNumber(entry1.getValue());
                    int instanceNumber2 = extractInstanceNumber(entry2.getValue());

                    return Integer.compare(instanceNumber1, instanceNumber2);
                })
                .collect(LinkedHashMap::new,
                    (map, entry) -> map.put(entry.getKey(), entry.getValue()),LinkedHashMap::putAll);
    }
    private int extractInstanceNumber(String metadata) {
        if (metadata == null || metadata.isEmpty())return Integer.MAX_VALUE;
        String instanceNumberStr = getTag(metadata, "0020,0013");
        if (instanceNumberStr == null || instanceNumberStr.trim().isEmpty()) {
            return Integer.MAX_VALUE;
        }
        try {
            return Integer.parseInt(instanceNumberStr.trim());
        } catch (NumberFormatException e) {
            return Integer.MAX_VALUE;
        }
    }
    
    @Override
    public String toString() {
        return seriesDescription;
    }
}
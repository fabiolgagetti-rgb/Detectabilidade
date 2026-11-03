package br.unifesp.fisicamedica;

import java.awt.BasicStroke;
import java.util.List;
import java.util.Map;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class TtfGraphManager {

    private static ChartPanel chartPanel;

    /**
     * Cria um ChartPanel vazio para ser adicionado ao layout da UI.
     */
    public static ChartPanel createChartPanel() {
        JFreeChart emptyChart = ChartFactory.createXYLineChart(
            "", "Frequency (cycles/mm)", "TTF", new XYSeriesCollection());
        chartPanel = new ChartPanel(emptyChart);
        return chartPanel;
    }

    /**
     * Atualiza o gráfico de acordo com o tipo de curva selecionada (MTF, LSF, ESF).
     * Agora utiliza List<Integer> para referenciar os IDs dos inserts.
     */
    public static void updateChart(List<Integer> selectedIDs, String curveType, TtfCalc.mtfResult result) {
        if (chartPanel == null || selectedIDs == null || selectedIDs.isEmpty() || result == null) return;

        Map<Integer, double[]> dataMap;
        String xLabel;
        String yLabel;
        boolean isMTF = false;

        switch (curveType) {
            case "LSF":
                dataMap = result.lsfs;
                xLabel = "Position (pixels)";
                yLabel = "LSF";
                break;
            case "ESF":
                dataMap = result.esfs;
                xLabel = "Position (pixels)";
                yLabel = "ESF (HU)";
                break;
            case "TTF":
            default:
                dataMap = result.ttfs;
                xLabel = "Frequency (cycles/mm)";
                yLabel = "TTF";
                isMTF = true;
                break;
        }

        if (dataMap == null || dataMap.isEmpty()) {
            clearChart();
            return;
        }
        
        XYSeriesCollection dataset = new XYSeriesCollection();
        double maxFreq = 0.0;

        for (int i = 0; i < selectedIDs.size(); i++) {
            Integer currentID = selectedIDs.get(i);
            double[] yData = dataMap.get(currentID);
            double[] freqData = result.freqs.get(currentID); 

            if (yData == null || yData.length == 0) continue;

            String seriesName = "Task " + (currentID + 1);

            XYSeries series = new XYSeries(seriesName);
            
            if (isMTF) {
                if (freqData == null || freqData.length == 0) continue;

                int limit = Math.min(freqData.length, yData.length);
                for (int i1 = 0; i1 < limit; i1++) {
                    series.add(freqData[i1], yData[i1]);
                }
                if (limit > 0 && freqData[limit - 1] > maxFreq) {
                    maxFreq = freqData[limit - 1];
                }
            } else {
                int limit = yData.length;
                int center = limit / 2; 
                double maxAbsVal = 1.0;
                
                if (curveType.equals("LSF")) {
                    maxAbsVal = 0.0;
                    for (int i1 = 0; i1 < limit; i1++) {
                        double absVal = Math.abs(yData[i1]);
                        if (absVal > maxAbsVal) {
                            maxAbsVal = absVal;
                        }
                    }
                    if (maxAbsVal == 0) maxAbsVal = 1.0;
                }
                
                for (int i1 = 0; i1 < limit; i1++) {
                    double xPos = (i1 - center) * 0.5; 
                    double yVal = yData[i1] / maxAbsVal; 
                    series.add(xPos, yVal);
                }
            }
            
            dataset.addSeries(series);
        }

        if (dataset.getSeriesCount() == 0) {
            clearChart();
            return;
        }

        JFreeChart chart = ChartFactory.createXYLineChart(
            "",
            xLabel,
            yLabel,
            dataset,
            PlotOrientation.VERTICAL,
            true, true, false
        );

        XYPlot plot = chart.getXYPlot();
        
            plot.getDomainAxis().setAutoRange(true);
            plot.getRangeAxis().setAutoRange(true);
        
        
        plot.setDomainGridlineStroke(new BasicStroke(0.5f));
        plot.setRangeGridlineStroke(new BasicStroke(0.5f));

        chartPanel.setChart(chart);
    }
    
    public static void clearChart() {
         if (chartPanel == null) return;
         JFreeChart emptyChart = ChartFactory.createXYLineChart("",
                 "Frequency (cycles/mm)", "TTF", new XYSeriesCollection());
         chartPanel.setChart(emptyChart);
    }
}
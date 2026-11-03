package br.unifesp.fisicamedica;

import java.awt.BasicStroke;
import java.awt.Color;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class NpsGraphManager {
    
    private static ChartPanel chartPanel;

    /**
     * Cria um ChartPanel inicialmente vazio para estruturar a UI.
     */
    public static ChartPanel createChartPanel() {
        JFreeChart emptyChart = ChartFactory.createXYLineChart(
                "", 
                "Frequency (cycles/mm)", "NPS (HU²·mm²)", 
                new XYSeriesCollection());
        
        chartPanel = new ChartPanel(emptyChart);
        chartPanel.setPreferredSize(new java.awt.Dimension(450, 300)); 
        return chartPanel;
    }

    /**
     * Atualiza o gráfico com a curva do Perfil Radial do NPS.
     */
    public static void updateChart(double[] frequencies, double[] radialProfile) {
        if (chartPanel == null) return;

        if (frequencies == null || radialProfile == null || frequencies.length == 0) {
            clearChart();
            return;
        }
        
        XYSeriesCollection dataset = new XYSeriesCollection();
        addSeriesToDataset(dataset, frequencies, radialProfile, "Radial NPS");

        JFreeChart chart = ChartFactory.createXYLineChart(
                "",
                "Frequency (cycles/mm)", "NPS (HU²·mm²)",
                dataset,
                PlotOrientation.VERTICAL,
                false,
                true, false
        );

        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();

        renderer.setSeriesShapesVisible(0, false);
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesPaint(0, new Color(0, 102, 204)); 
        renderer.setSeriesStroke(0, new BasicStroke(1.5f));

        plot.setRenderer(renderer);
        
        plot.setDomainGridlineStroke(new BasicStroke(0.5f));
        plot.setRangeGridlineStroke(new BasicStroke(0.5f));
        
        plot.getDomainAxis().setRange(0, frequencies[frequencies.length - 1]);

        chartPanel.setChart(chart);
    }

    private static void addSeriesToDataset(XYSeriesCollection dataset, double[] xValues, double[] yValues, String seriesName) {
        XYSeries series = new XYSeries(seriesName);
        int limit = Math.min(xValues.length, yValues.length); 
        
        for (int i = 0; i < limit; i++) {
            series.add(xValues[i], yValues[i]);
        }
        dataset.addSeries(series);
    }
    
    /**
     * Limpa o gráfico em caso de falha no cálculo ou reset.
     */
    public static void clearChart() {
        if (chartPanel == null) return;
        JFreeChart emptyChart = ChartFactory.createXYLineChart(
                "", "Frequency (cycles/mm)", "NPS (HU²·mm²)", new XYSeriesCollection());
        chartPanel.setChart(emptyChart);
    }
}
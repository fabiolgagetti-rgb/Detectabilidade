package br.unifesp.fisicamedica;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class GraphManager {
    
    private static ChartPanel chartPanel;

    public static ChartPanel createChartPanel() {
        chartPanel = new ChartPanel(null);

        return chartPanel;
    }

    public static void updateChart(double[] frequencies, double[] radialProfile) {

      	if (frequencies==null || radialProfile == null ) return;
        
        XYSeriesCollection  dataset = new XYSeriesCollection ();
        addSeriesToDataset(dataset, frequencies, radialProfile, "NPS");

        JFreeChart chart = ChartFactory.createXYLineChart(
                "",
                "Frequency (mm-1)", "NPS (HU²mm²)",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false
        );
        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();

        for (int i = 0; i < plot.getSeriesCount(); i++) {
            renderer.setSeriesShapesVisible(i, false);
            renderer.setSeriesLinesVisible(i, true);
        }
        plot.setRenderer(renderer);

        chartPanel.setChart(chart);
        chartPanel.revalidate();
        chartPanel.repaint();
    }
    public static void addSeriesToDataset(XYSeriesCollection dataset, double[] xValues, double[] yValues, String seriesName) {
        XYSeries series = new XYSeries(seriesName);
        for (int i = 0; i < xValues.length; i++) {
            series.add(xValues[i], yValues[i]);
        }
        dataset.addSeries(series);
    }
}
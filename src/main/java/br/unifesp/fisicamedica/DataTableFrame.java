package br.unifesp.fisicamedica;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.util.List;

public class DataTableFrame extends JFrame {
    private static final long serialVersionUID = 1L;
    private DefaultTableModel model;
    private JTable table;

    /**
     * Construtor Universal para Tabelas de Dados.
     * @param title O título da janela.
     * @param columnNames Array com os nomes das colunas.
     * @param dataColumns Lista de arrays de dados. O índice 0 geralmente é o eixo X (Frequência).
     */
    public DataTableFrame(String title, String[] columnNames, List<double[]> dataColumns) {
        setTitle(title);
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setLayout(new BorderLayout());
        setPreferredSize(new Dimension(600, 400));

        model = new DefaultTableModel() {
            private static final long serialVersionUID = 1L;
            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };

        for (String colName : columnNames) {
            model.addColumn(colName);
        }

        fillTable(dataColumns);

        table = new JTable(model);
        configureTable();

        JScrollPane scrollPane = new JScrollPane(table);
        add(scrollPane, BorderLayout.CENTER);
        
        pack();
        setLocationRelativeTo(null);
    }

    private void fillTable(List<double[]> dataColumns) {
        model.setRowCount(0);
        if (dataColumns == null || dataColumns.isEmpty()) return;

        int numRows = dataColumns.get(0).length;

        for (int i = 0; i < numRows; i++) {
            Object[] rowData = new Object[dataColumns.size()];
            
            if (dataColumns.get(0) != null && i < dataColumns.get(0).length) {
                rowData[0] = String.format("%.3f", dataColumns.get(0)[i]);
            } else {
                rowData[0] = "-";
            }

            for (int col = 1; col < dataColumns.size(); col++) {
                double[] colData = dataColumns.get(col);
                if (colData != null && i < colData.length) {
                    rowData[col] = String.format("%.6g", colData[i]); 
                } else {
                    rowData[col] = "-";
                }
            }
            model.addRow(rowData);
        }
    }

    private void configureTable() {
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.setFillsViewportHeight(true);
        for (int i = 0; i < table.getColumnCount(); i++) {
            table.getColumnModel().getColumn(i).setPreferredWidth(120);
        }
    }
}
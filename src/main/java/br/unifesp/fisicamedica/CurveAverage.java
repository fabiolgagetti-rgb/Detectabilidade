package br.unifesp.fisicamedica;

import java.util.ArrayList;
import java.util.List;

public class CurveAverage {

	public static double[] averageEsfs(List<double[]> esfCurves) {
	    if (esfCurves == null || esfCurves.isEmpty()) return null;

	    List<double[]> validCurves = new ArrayList<>();
	    for (double[] esf : esfCurves) {
	        if (esf != null) validCurves.add(esf);
	    }
	    
	    if (validCurves.isEmpty()) return null;

	    int size = validCurves.get(0).length;
	    double[] sum = new double[size];
	    
	    for (double[] esf : validCurves) {   	
	        for (int i = 0; i < size; i++) {
	            sum[i] += esf[i];
	        }
	    }

	    double[] average = new double[size];
	    for (int i = 0; i < size; i++) {
	        average[i] = sum[i] / validCurves.size();
	    }
	    
	    return average;
	}
}
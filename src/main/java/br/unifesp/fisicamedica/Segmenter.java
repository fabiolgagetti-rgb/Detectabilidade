package br.unifesp.fisicamedica;

import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.measure.Calibration;
import ij.plugin.RoiEnlarger;
import ij.plugin.filter.EDM;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.process.FloodFiller;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class Segmenter {

	public static class DetectorResult {
	    List<Double> HU = new ArrayList<>();
	    List<Double> StdDev = new ArrayList<>();
	    List<Double> Threshold = new ArrayList<>();
	    List<Double> radius = new ArrayList<>();
	    List<List<Roi>> rois = new ArrayList<>();
	    List<Roi> roi = new ArrayList<>();
	    List<Point> Pos = new ArrayList<>();

	    List<List<Double>> modeContrast = new ArrayList<>();
	    List<List<Double>> meanI0 = new ArrayList<>();
	    List<List<Double>> meanI1 = new ArrayList<>();

	    // Per-structure contrast / noise / CNR, aligned with 'roi' (one entry per
	    // accepted structure). Registered at detection time for reporting in the UI.
	    //   contrast : full edge step height (HU-calibrated units of the image)
	    //   noise    : std-dev of the background band just outside the shape
	    //   cnr      : contrast / noise
	    public List<Double> contrast = new ArrayList<>();
	    public List<Double> noise = new ArrayList<>();
	    public List<Double> cnr = new ArrayList<>();

	    // Option A (PAE edge normals), unit normals [nPoints][2] aligned 1:1 with the
	    // final ROI polygon points, taken directly from the partial-area detector.
	    // Lets the ESF step use each point's own analytic normal instead of a
	    // neighbour-fit, so dropped points ("holes") no longer corrupt neighbouring
	    // normals. May be null/empty when unavailable.
	    //   edgeNormalsPerRoi : per-insert   (aligned with 'roi'),  set in processRecord
	    //   edgeNormals       : per-insert, per-slice (aligned with 'rois'), set in edgeFinder
	    List<double[][]> edgeNormalsPerRoi = new ArrayList<>();
	    List<List<double[][]>> edgeNormals = new ArrayList<>();
	}

    private static class OtsuResult {
        int threshold;
        double max_bcv;
    }

    // ------------------------------------------------------------------
    // Adaptive parameter selection for the subpixel edge detector.
    //
    // Replaces two previously hard-coded choices in processRecord:
    //   (1) a fixed Gaussian pre-blur of sigma = 1, and
    //   (2) gradThreshold = 0.5 * |mode - thresholdForEdge|.
    //
    // (1) The optimal pre-blur depends on the contrast-to-noise ratio (CNR):
    //     clean/soft-kernel data wants little blur (localization), noisy
    //     bone/lung-kernel data wants more (denoising). sigma is fit to the
    //     measured localization/denoising crossover:
    //         sigma = clamp(2.0 - 0.75*log10(CNR), 0.4, 1.6)
    //     (CNR 100 -> 0.5, 20 -> 1.0, 10 -> 1.25, 5 -> 1.5).
    //
    // (2) The old threshold mixed intensity units with a gradient-magnitude
    //     comparison and ignored noise. The detector thresholds the magnitude
    //     of a central difference of the 3x3-smoothed image, so the threshold
    //     is set from the wand-ROI noise std-dev: a central difference of two
    //     independent samples has std ~0.7*noise, and we sit ~3.5x above that
    //     floor (with a small contrast-relative minimum for the clean case).
    // ------------------------------------------------------------------

    private static double blurSigmaForCNR(double cnr) {
        if (!(cnr > 0) || Double.isInfinite(cnr)) return 0.4;
        double s = 2.0 - 0.75 * Math.log10(cnr);
        return Math.max(0.4, Math.min(1.6, s));
    }

    /**
     * Estimate image noise from a thin band JUST OUTSIDE the detected shape.
     *
     * Builds band = enlarge(shape, +bandPx) minus the shape interior, and returns
     * the std-dev of the ORIGINAL image there. This region is background (truly
     * homogeneous), so it reflects genuine image noise and is immune to internal
     * structures inside the object (which would otherwise inflate an interior-based
     * std and over-blur the edge). Falls back to the shrunk-interior std if the
     * outside band is empty or touches the image border badly.
     */
    private static double noiseOutsideShape(ImageProcessor origIp, Roi shape) {
        try {
            // Leave a GAP beyond the shape so the band clears the (blurred) edge
            // transition entirely, then take a band of width 'bandPx' past the gap.
            int gapPx = 6;      // clears the edge + its blur
            int bandPx = 8;     // thickness of the pure-background sampling band
            Roi inner = RoiEnlarger.enlarge(shape, gapPx);
            Roi outer = RoiEnlarger.enlarge(shape, gapPx + bandPx);
            ByteProcessor bandMask = new ByteProcessor(origIp.getWidth(), origIp.getHeight());
            bandMask.setColor(255);
            bandMask.fill(outer);
            bandMask.setColor(0);
            bandMask.fill(inner);                 // carve out shape+gap -> only the far band remains

            double sum = 0, sum2 = 0; long n = 0;
            for (int y = 0; y < origIp.getHeight(); y++) {
                for (int x = 0; x < origIp.getWidth(); x++) {
                    if (bandMask.get(x, y) != 0) {
                        double v = origIp.getf(x, y);
                        sum += v; sum2 += v * v; n++;
                    }
                }
            }
            if (n >= 20) {
                double mean = sum / n;
                double var = sum2 / n - mean * mean;
                if (var > 0) return Math.sqrt(var);
            }
        } catch (Exception ignore) {
            // fall through to interior fallback
        }
        // Fallback: previous behavior (shrunk interior std) if the outside band failed.
        Roi noiseRoi = RoiEnlarger.enlarge(shape, -5);
        origIp.setRoi(noiseRoi);
        double sd = origIp.getStatistics().stdDev;
        origIp.resetRoi();
        return sd;
    }

    public static DetectorResult edgeFinder(ImageStack stack, Calibration cal, Roi userRoi, boolean isAuto, boolean avgSlices) {
        ImageProcessor media = calcularMedia(stack);
        DetectorResult dR;

        if (isAuto) {
            double recordSizePx = media.getWidth() * 0.1;
            double overlap = 0.5;
            dR = autoDetect(media, cal, recordSizePx, overlap);
        } else {
            if (userRoi == null) {
                throw new IllegalArgumentException("Draw a ROI in the image or use the automatic detection method");
            }
            dR = roiDetect(media, media, cal, userRoi);
        }

        int nSlices = stack.size();
        dR.rois = new ArrayList<>();

        dR.edgeNormals = new ArrayList<>();

        for (int i = 0; i < dR.roi.size(); i++) {
            List<Roi> sliceRois = new ArrayList<>();
            List<double[][]> sliceNormals = new ArrayList<>();
            Roi avgRoi = dR.roi.get(i);
            double[][] avgNormals = (i < dR.edgeNormalsPerRoi.size()) ? dR.edgeNormalsPerRoi.get(i) : null;
            Rectangle bounds = dR.roi.get(i).getBounds();

            for (int z = 1; z <= nSlices; z++) {
                if (avgSlices) {
                    sliceRois.add(avgRoi);
                    sliceNormals.add(avgNormals);
                } else {
                    ImageProcessor sliceIp = stack.getProcessor(z);
                    ImageProcessor blurredSlice = sliceIp.duplicate();
                    blurredSlice.blurGaussian(1);

                    DetectorResult tempResult = new DetectorResult();

                    int padding = (int) Math.max(10, (bounds.getWidth()+bounds.getHeight())*0.1);
                    int rx = Math.max(0, bounds.x - padding);
                    int ry = Math.max(0, bounds.y - padding);
                    int rw = Math.min(sliceIp.getWidth() - rx, bounds.width + 2 * padding);
                    int rh = Math.min(sliceIp.getHeight() - ry, bounds.height + 2 * padding);
                    Roi sliceRoi = new Roi(rx, ry, rw, rh);

                    blurredSlice.setRoi(userRoi);
                    ImageStatistics roiStats = blurredSlice.getStatistics();
                    double roiMode = (roiStats.histMin + roiStats.mode * roiStats.binSize);

                    processRecord(blurredSlice, sliceIp, cal, sliceRoi, tempResult, null, roiMode);

                    if (!tempResult.roi.isEmpty()) {
                        sliceRois.add(tempResult.roi.get(0));
                        sliceNormals.add(tempResult.edgeNormalsPerRoi.isEmpty()
                                ? null : tempResult.edgeNormalsPerRoi.get(0));
                    }
                }
            }
            dR.rois.add(sliceRois);
            dR.edgeNormals.add(sliceNormals);
        }
        return dR;
    }

    private static DetectorResult autoDetect(ImageProcessor origIp, Calibration cal, double recordSizePx, double overlap) {
        DetectorResult dR = new DetectorResult();

        ImageProcessor blurredIp = origIp.duplicate();
        blurredIp.blurGaussian(1);

        PolygonRoi phantomRoi = segmenter(origIp, cal, -200, 5000, 0.2, false, 0, 20);

        blurredIp.setRoi(phantomRoi);
        ImageStatistics phantomStats = blurredIp.getStatistics();
        double mode = (phantomStats.histMin + phantomStats.mode * phantomStats.binSize);

        ByteProcessor mask = new ByteProcessor(blurredIp.getWidth(), blurredIp.getHeight());
        mask.setColor(0);
        mask.fill();

        Rectangle searchArea;


        if (phantomRoi != null) {
            searchArea = phantomRoi.getBounds();
        } else {
            searchArea = new Rectangle(0, 0, blurredIp.getWidth(), blurredIp.getHeight());;
        }

        int step = (int) Math.max(1, Math.round(recordSizePx * (1 - overlap)));
        int margin = (int) (recordSizePx * 0.8);
        int startX = Math.max(0, searchArea.x + margin);
        int startY = Math.max(0, searchArea.y + margin);
        int endX = (int) Math.min(blurredIp.getWidth() - recordSizePx, searchArea.x + searchArea.width - margin);
        int endY = (int) Math.min(blurredIp.getHeight() - recordSizePx, searchArea.y + searchArea.height - margin);

        for (int y = startY; y <= endY; y += step) {
            for (int x = startX; x <= endX; x += step) {
            	if (phantomRoi != null) {
                    boolean topLeft = phantomRoi.containsPoint(x, y);
                    boolean topRight = phantomRoi.containsPoint(x + recordSizePx, y);
                    boolean bottomLeft = phantomRoi.containsPoint(x, y + recordSizePx);
                    boolean bottomRight = phantomRoi.containsPoint(x + recordSizePx, y + recordSizePx);

                    if (!topLeft || !topRight || !bottomLeft || !bottomRight) {
                        continue;
                    }
                }
                if (containsObject(mask, x, y, recordSizePx)) continue;

                Roi recordBounds = new Roi(x, y, (int)recordSizePx, (int)recordSizePx);
                processRecord(blurredIp, origIp, cal, recordBounds, dR, mask, mode);
            }
        }
        return dR;
    }

    private static DetectorResult roiDetect(ImageProcessor ipOrig, ImageProcessor media, Calibration cal, Roi userRoi) {
        DetectorResult result = new DetectorResult();

        ImageProcessor blurredIp = media.duplicate();
        blurredIp.blurGaussian(1);

        blurredIp.setRoi(userRoi);
        ImageStatistics roiStats = blurredIp.getStatistics();
        double roiMode = (roiStats.histMin + roiStats.mode * roiStats.binSize);
        blurredIp.resetRoi();
        processRecord(blurredIp, ipOrig, cal, userRoi, result, null, roiMode);

        return result;
    }

    private static void processRecord(ImageProcessor blurredIp, ImageProcessor origIp, Calibration cal,
                                               Roi recordBounds, DetectorResult dR, ByteProcessor processedMask, double mode) {
        int imgWidth = blurredIp.getWidth();
        int imgHeight = blurredIp.getHeight();

        blurredIp.setRoi(recordBounds);
        ImageProcessor window = blurredIp.crop();
        blurredIp.resetRoi();

        ImageStatistics stats = window.getStatistics();
        OtsuResult or = Otsu(stats.histogram);
        double bcv = or.max_bcv;
        int temp = or.threshold;
        double thresholdForEdge = (stats.histMin + temp * stats.binSize + stats.binSize / 2.0);

        if (bcv < 3000) {
            EdgeDiagnostics.record(EdgeDiagnostics.Reason.BCV_GATE);
            return;
        }

        double minHU, maxHU;
        if (thresholdForEdge >= mode) {
            minHU = thresholdForEdge; maxHU = Double.MAX_VALUE;
        } else {
            minHU = -Double.MAX_VALUE; maxHU = thresholdForEdge;
        }

        Wand wand = new Wand(blurredIp);
        boolean found = false;

        int w = window.getWidth();
        int h = window.getHeight();
        int midX = w / 2;
        int midY = h / 2;

        for (int stepX = 0; stepX < w && !found; stepX++) {
            int i = (midX + stepX) % w;
            for (int stepY = 0; stepY < h && !found; stepY++) {
                int j = (midY + stepY) % h;
                double val = window.getValue(i, j);

                if (val > minHU && val < maxHU) {
                    int gx = (int) (recordBounds.getXBase() + i);
                    int gy = (int) (recordBounds.getYBase() + j);
                    wand.autoOutline(gx, gy, (int)Math.round(minHU), (int)Math.round(Math.min(maxHU, 65535)), 8 | 256);
                    if (wand.npoints > 0) {
                        found = true;
                    }
                }
            }
        }
        if (!found) {
            EdgeDiagnostics.record(EdgeDiagnostics.Reason.WAND_NOT_FOUND);
            return;
        }

        boolean touchesOriginalRoi = false;
        for (int p = 0; p < wand.npoints; p++) {
            if (recordBounds.contains((int)wand.xpoints[p], (int)wand.ypoints[p])) {
                touchesOriginalRoi = true;
                break;
            }
        }
        if (!touchesOriginalRoi) {
            EdgeDiagnostics.record(EdgeDiagnostics.Reason.WAND_OFF_ROI);
            return;
        }

        PolygonRoi wandRoi = new PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.FREEROI);

        Rectangle bounds = wandRoi.getBounds();
        if(bounds.x <= 0 || bounds.y <= 0 || bounds.x + bounds.width >= imgWidth || bounds.y + bounds.height >= imgHeight) {
             EdgeDiagnostics.record(EdgeDiagnostics.Reason.BOUNDS_TOUCH_IMAGE);
             return;
        }
        double sumX = 0, sumY = 0;
        int ptCount = wandRoi.getFloatPolygon().npoints;
        float[] xPts = wandRoi.getFloatPolygon().xpoints;
        float[] yPts = wandRoi.getFloatPolygon().ypoints;

        for (int p = 0; p < ptCount; p++) {
            sumX += xPts[p];
            sumY += yPts[p];
        }
        double centerX = sumX / ptCount;
        double centerY = sumY / ptCount;


        int padding = 5;
        int rx = Math.max(0, bounds.x - padding);
        int ry = Math.max(0, bounds.y - padding);
        int rw = Math.min(origIp.getWidth() - rx, bounds.width + 2 * padding);
        int rh = Math.min(origIp.getHeight() - ry, bounds.height + 2 * padding);

        // ------------------------------------------------------------------
        // Estimate the edge contrast and the local noise, then choose the
        // detector's pre-blur sigma and gradient threshold adaptively.
        //
        // contrast: full step height. thresholdForEdge sits mid-edge, so
        //           |mode - thresholdForEdge| is ~ half the contrast.
        // noise:    std-dev of a thin band JUST OUTSIDE the detected shape. The
        //           previous estimate used the whole wand-ROI INTERIOR, which for
        //           objects with internal structure (e.g. a bone-marrow disk inside
        //           a vertebra) counts that structure's step as "noise": the std-dev
        //           inflates, CNR collapses, and the adaptive pre-blur over-smooths
        //           the edge -> the measured MTF droops. The region just outside the
        //           shape is genuinely homogeneous (background), so it isolates true
        //           image noise regardless of what is inside the object.
        // ------------------------------------------------------------------
        double contrast = 2.0 * Math.abs(mode - thresholdForEdge);

        double noise = noiseOutsideShape(origIp, wandRoi);
        double cnr = noise > 0 ? contrast / noise : Double.POSITIVE_INFINITY;
        double blurSigma = blurSigmaForCNR(cnr);

        // Detector input: an ADAPTIVELY blurred copy of the ORIGINAL image,
        // cropped to the padded bounds. (The sigma=1 blurredIp is kept only for
        // the wand/Otsu/segmentation steps above.)
        ImageProcessor detectorFull = origIp.duplicate();
        detectorFull.blurGaussian(blurSigma);
        detectorFull.setRoi(rx, ry, rw, rh);
        ImageProcessor cropDetector = detectorFull.crop();
        detectorFull.resetRoi();

        FloatPolygon finalSubpixelPolygon = new FloatPolygon();

        // Gradient threshold, replacing the old 0.5*|mode - thresholdForEdge|.
        //
        // The detector thresholds the magnitude of a central difference of the
        // smoothed image. The threshold must stay BELOW the true edge gradient
        // at all CNR (otherwise it wipes out detection on low-contrast edges),
        // while rising enough with noise to reject texture.
        //
        // Previously the cap was 0.6*(contrast/3), modelling the edge gradient as
        // ~contrast/3. That model is wrong for CT: measured on MTF/NPS phantoms
        // (SOFT and Lung, Bujila 2017 coefficients) the true peak edge gradient
        // is ~contrast..3*contrast (SOFT ~304, Lung ~574 for contrast 100 due to
        // edge enhancement), i.e. ~10-17x larger. The contrast/3 cap therefore
        // forced the threshold BELOW the noise floor at high noise, admitting
        // stray noise edges (worst on SOFT at high SD). Capping at 0.5*contrast
        // tracks the real edge gradient: it stays well below the true edge while
        // still bounding the threshold, halves stray edges at high SD on SOFT,
        // and is neutral for Lung. Detection at low contrast/high noise is kept
        // alive by the single 0.6x fallback below.
        double gradThreshold = Math.min(
                Math.max(0.9 * noise, 0.12 * contrast),
                0.5 * contrast);

        // Restoration iterations reduced 10 -> 3 (convergence; avoids over-sharpening
        // blurred/soft-kernel edges). The 30x threshold-lowering loop is replaced by
        // a single, bounded fallback now that the threshold is measured rather than guessed.
        PartialAreaEdgeDetector.Result paeResult =
                PartialAreaEdgeDetector.detectAndRestoreIterative(cropDetector, gradThreshold, 2, 1);
        if (paeResult.edges.size() < ptCount) {
            paeResult = PartialAreaEdgeDetector.detectAndRestoreIterative(
                    cropDetector, gradThreshold * 0.6, 2, 1);
        }

        FloatProcessor restoredCrop = paeResult.restoredImage;
        Wand cleanWand = new Wand(restoredCrop);

        if (thresholdForEdge >= mode) { minHU = thresholdForEdge; maxHU = Double.MAX_VALUE; }
        else { minHU = -Double.MAX_VALUE; maxHU = thresholdForEdge; }

        cleanWand.autoOutline((int)(centerX - rx), (int)(centerY - ry), minHU, maxHU, 8 | 256);
        // PAE normals collected in lockstep with finalSubpixelPolygon points (Option A).
        List<double[]> finalNormals = new ArrayList<>();

        String edgeStrategy = edgeStrategy();

        if ("primary".equals(edgeStrategy)) {
            // ---- Strategy C: build the polygon from the PAE edge points directly ----
            // The PAE edge list IS the sub-pixel edge; it is gap-free wherever a true
            // edge exists. Rather than walk the wand and snap-or-drop (which creates
            // holes), we keep every PAE point that belongs to this insert's edge,
            // order them around the contour, and carry each point's analytic normal.
            // The wand is used only to (a) establish the interior reference (centroid)
            // for ordering and (b) reject stray edges far from the detected outline.
            // centroid in CROP coordinates (PAE edges and wand are crop-relative)
            buildPolygonFromPae(paeResult, cleanWand, centerX - rx, centerY - ry, rx, ry,
                    finalSubpixelPolygon, finalNormals);
        } else if ("fill".equals(edgeStrategy) && cleanWand.npoints > 0) {
            // ---- Strategy B (DEFAULT): hole-free wand anchoring + spike guard + smoothing ----
            // Keeps the wand-anchored geometry but eliminates holes at the source:
            // a wider match tolerance fills most gaps, a shape/deviation test rejects
            // spurious "spike" matches, genuine gaps are bridged by interpolation
            // (never by a bad snap), and the contour is lightly smoothed at the end.
            // Validated as the most robust across SOFT and LUNG at all noise levels.
            buildPolygonFilled(paeResult, cleanWand, rx, ry, finalSubpixelPolygon, finalNormals);
        } else if (("legacy".equals(edgeStrategy) || "normals".equals(edgeStrategy)) && cleanWand.npoints > 0) {
            // ---- Strategy A ("normals") / legacy: wand snap-or-drop polygon. ----
            // "legacy" uses neighbour-fit normals in the ESF step; "normals" (A) uses
            // the per-point PAE normals collected here. Both can leave holes.
            List<Point> denseWand = new ArrayList<>();
            for (int i = 0; i < cleanWand.npoints; i++) {
                int x1 = cleanWand.xpoints[i];
                int y1 = cleanWand.ypoints[i];
                int nextIdx = (i + 1) % cleanWand.npoints;
                int x2 = cleanWand.xpoints[nextIdx];
                int y2 = cleanWand.ypoints[nextIdx];

                denseWand.add(new Point(x1, y1));

                int steps = Math.max(Math.abs(x2 - x1), Math.abs(y2 - y1));
                if (steps > 1) {
                    for (int s = 1; s < steps; s++) {
                        int interpX = x1 + s * (x2 - x1) / steps;
                        int interpY = y1 + s * (y2 - y1) / steps;
                        denseWand.add(new Point(interpX, interpY));
                    }
                }
            }
            for (Point p : denseWand) {
                double px = p.x;
                double py = p.y;
                PartialAreaEdgeDetector.EdgePixel bestEdge = null;
                double minDist = Double.MAX_VALUE;

                for (PartialAreaEdgeDetector.EdgePixel ep : paeResult.edges) {
                    double dist = Math.hypot(ep.x - px, ep.y - py);
                    if (dist < minDist && dist < 3) {
                        minDist = dist;
                        bestEdge = ep;
                    }
                }
                EdgeDiagnostics.recordAnchor(bestEdge == null);
                if (bestEdge != null) {
                    finalSubpixelPolygon.addPoint((float)(bestEdge.x + rx), (float)(bestEdge.y + ry));
                    // store this point's analytic PAE normal (unit vector)
                    finalNormals.add(new double[]{ bestEdge.nx, bestEdge.ny });
                }
            }
        }
        double circularity = getCircularity(finalSubpixelPolygon);
        if (circularity < 0.3) {
            EdgeDiagnostics.record(EdgeDiagnostics.Reason.LOW_CIRCULARITY);
            return;
        }

        PolygonRoi finalRoi =  new PolygonRoi(finalSubpixelPolygon, Roi.POLYLINE);

        List<Double> edgeContrasts = new ArrayList<>();

        for (PartialAreaEdgeDetector.EdgePixel ep : paeResult.edges) {
             double[] flat = contrastAwayFromEdge(origIp, ep.x + rx, ep.y + ry, ep.nx, ep.ny, ep.i0, ep.i1);
             double valI0 = cal.getCValue(flat[0]);
             double valI1 = cal.getCValue(flat[1]);
             double c = Math.abs(valI1 - valI0);
             edgeContrasts.add(c);
        }

        List<Double> expectedContrasts = new ArrayList<>();
        List<Double> expectedI0 = new ArrayList<>();
        List<Double> expectedI1 = new ArrayList<>();

        if (!edgeContrasts.isEmpty()) {
             double minC = Collections.min(edgeContrasts);
             double maxC = Collections.max(edgeContrasts);

             double binSize = 15.0;
             int numBins = (int) Math.ceil((maxC - minC) / binSize) + 1;

             int[] histogram = new int[numBins];
             double[] sumInBin = new double[numBins];

             double[] sumI0InBin = new double[numBins];
             double[] sumI1InBin = new double[numBins];

             for (PartialAreaEdgeDetector.EdgePixel ep : paeResult.edges) {
                 double[] flat = contrastAwayFromEdge(origIp, ep.x + rx, ep.y + ry, ep.nx, ep.ny, ep.i0, ep.i1);
                 double valI0 = cal.getCValue(flat[0]);
                 double valI1 = cal.getCValue(flat[1]);
                 double c = Math.abs(valI1 - valI0);

                 int bin = (int) ((c - minC) / binSize);
                 histogram[bin]++;
                 sumInBin[bin] += c;
                 sumI0InBin[bin] += valI0;
                 sumI1InBin[bin] += valI1;
             }

             int maxBin1 = -1;
             int maxCount1 = 0;
             for (int i = 0; i < numBins; i++) {
                 if (histogram[i] > maxCount1) {
                     maxCount1 = histogram[i];
                     maxBin1 = i;
                 }
             }

             if (maxCount1 > 0) {
                 expectedContrasts.add(sumInBin[maxBin1] / maxCount1);
                 expectedI0.add(sumI0InBin[maxBin1] / maxCount1);
                 expectedI1.add(sumI1InBin[maxBin1] / maxCount1);
             }

             int maxBin2 = -1;
             int maxCount2 = 0;
             int minSeparationBins = (int) (20.0 / binSize);

             for (int i = 0; i < numBins; i++) {
                 if (Math.abs(i - maxBin1) > minSeparationBins) {
                     if (histogram[i] > maxCount2) {
                         maxCount2 = histogram[i];
                         maxBin2 = i;
                     }
                 }
             }

             if (maxBin2 != -1 && maxCount2 >= (maxCount1 * 0.2)) {
                 expectedContrasts.add(sumInBin[maxBin2] / maxCount2);
                 expectedI0.add(sumI0InBin[maxBin2] / maxCount2);
                 expectedI1.add(sumI1InBin[maxBin2] / maxCount2);
             }
        }

        if (expectedContrasts.isEmpty() || expectedContrasts.get(0) < 5) {
            EdgeDiagnostics.record(EdgeDiagnostics.Reason.CONTRAST_GATE);
            return;
        }
        EdgeDiagnostics.record(EdgeDiagnostics.Reason.ACCEPTED);

        dR.modeContrast.add(expectedContrasts);
        dR.meanI0.add(expectedI0);
        dR.meanI1.add(expectedI1);

        blurredIp.setRoi(wandRoi);
        ImageStatistics stats1 = blurredIp.getStatistics();
        double insertHU = Math.round(cal.getCValue(stats1.histMin + stats1.mode * stats1.binSize));

        dR.HU.add(insertHU);
        dR.Threshold.add(cal.getCValue(thresholdForEdge));

        dR.roi.add(finalRoi);
        dR.edgeNormalsPerRoi.add(finalNormals.toArray(new double[0][]));

        // Register per-structure contrast/noise/CNR for later UI display.
        dR.contrast.add(contrast);
        dR.noise.add(noise);
        dR.cnr.add(cnr);

        dR.Pos.add(new Point((int)Math.round(centerX), (int)Math.round(centerY)));

        if (processedMask != null) {
            processedMask.setColor(255);
            processedMask.fill(finalRoi);
        }
    }

    /**
     * Option C: build the edge polygon directly from the PAE edge points.
     *
     * The PAE edge list is the sub-pixel edge and is gap-free wherever a true
     * edge exists, so using it directly avoids the "holes" the wand snap-or-drop
     * anchoring produces. We:
     *   1. keep PAE points near the wand outline (reject strays far from it),
     *   2. order them around the contour by polar angle about the wand centroid
     *      (the inserts are star-shaped about their centre; the circularity gate
     *      downstream already assumes near-circular edges),
     *   3. emit polygon points + each point's own analytic PAE normal.
     *
     * Because each point carries its own normal, ordering only affects polygon
     * topology (used by the circularity gate), never the ESF normals, so a mild
     * mis-order cannot corrupt the ESF.
     *
     * @param paeResult detector output (edge pixels in crop coordinates)
     * @param wand      the clean wand outline (crop coords) -> centroid + gating
     * @param centerX,centerY wand-polygon centroid in CROP coordinates
     * @param rx,ry     crop origin offset to map back to image coordinates
     * @param outPoly   receives the ordered sub-pixel polygon (image coords)
     * @param outNormals receives the per-point unit normals, aligned 1:1 with outPoly
     */
    private static void buildPolygonFromPae(PartialAreaEdgeDetector.Result paeResult, Wand wand,
            double centerX, double centerY, int rx, int ry,
            FloatPolygon outPoly, List<double[]> outNormals) {
        if (paeResult == null || paeResult.edges == null || paeResult.edges.isEmpty()) return;

        // Interior reference: centroid of the wand outline (crop coords). centerX/Y
        // are already that centroid mapped from the wand polygon; use directly.
        double cx = centerX, cy = centerY;

        // Gating radius: median wand radius +/- a tolerance, to drop stray PAE edges
        // (e.g. internal-structure edges) that are not on this insert's boundary.
        double meanR = 0; int wn = Math.max(1, wand.npoints);
        for (int i = 0; i < wand.npoints; i++) {
            meanR += Math.hypot(wand.xpoints[i] - cx, wand.ypoints[i] - cy);
        }
        meanR /= wn;
        double rTol = Math.max(5.0, 0.5 * meanR);   // generous band around the wand radius

        // Collect in-band PAE points. The PAE edge is ~1-2 px thick and the detector
        // may emit several points at nearly the same polar angle (inner/outer side,
        // restoration duplicates). A raw polar-angle sort of those zig-zags radially
        // and inflates the perimeter (crashing the circularity gate). So bin by angle
        // and keep ONE representative per bin: the point whose radius is closest to the
        // local contour radius, giving a clean single-valued, ordered contour.
        int nBins = Math.max(180, (int) (2 * Math.PI * meanR));   // ~<=1 point/px of arc
        PartialAreaEdgeDetector.EdgePixel[] best = new PartialAreaEdgeDetector.EdgePixel[nBins];
        double[] bestErr = new double[nBins];
        java.util.Arrays.fill(bestErr, Double.MAX_VALUE);

        for (PartialAreaEdgeDetector.EdgePixel ep : paeResult.edges) {
            double r = Math.hypot(ep.x - cx, ep.y - cy);
            if (Math.abs(r - meanR) > rTol) continue;             // reject strays
            double ang = Math.atan2(ep.y - cy, ep.x - cx);        // [-pi, pi)
            int bin = (int) ((ang + Math.PI) / (2 * Math.PI) * nBins);
            if (bin < 0) bin = 0; else if (bin >= nBins) bin = nBins - 1;
            double err = Math.abs(r - meanR);                     // closeness to contour radius
            if (err < bestErr[bin]) { bestErr[bin] = err; best[bin] = ep; }
        }

        int used = 0;
        for (int b = 0; b < nBins; b++) {                         // already in angular order
            PartialAreaEdgeDetector.EdgePixel ep = best[b];
            if (ep == null) continue;                             // empty bin: no PAE edge there
            outPoly.addPoint((float) (ep.x + rx), (float) (ep.y + ry));
            outNormals.add(new double[]{ ep.nx, ep.ny });
            used++;
        }
        // Every emitted point is a real PAE edge point (no snap-or-drop holes). Empty
        // angular bins are genuine edge absences, not anchoring holes; record no holes.
        for (int i = 0; i < used; i++) EdgeDiagnostics.recordAnchor(false);
    }

    /**
     * Option B: hole-FREE wand anchoring with spike guard and smoothing.
     *
     * Keeps the wand-anchored contour (minimal change to the trusted geometry) but
     * removes the "holes" the strict 3 px snap-or-drop leaves. Pipeline:
     *
     *   1. Densify the wand outline (as before) to get an ordered point sequence.
     *   2. For each densified point, snap to the nearest PAE edge within a WIDER
     *      tolerance (ttf.fillTolPx, default 6). This alone fills most gaps but can
     *      introduce "spikes": a point near a PAE gap may grab a pixel from another
     *      part of the contour, jutting off the smooth edge.
     *   3. SPIKE GUARD: using the wand centroid as reference, reject any snapped
     *      point whose radius deviates from the LOCAL median radius by more than
     *      ttf.fillSpikePx (default 2.5). Rejected points become gaps (not spikes).
     *   4. FILL gaps by linear interpolation of BOTH position and normal between the
     *      nearest accepted points on each side (circular). A genuine gap is bridged
     *      smoothly instead of forcing a bad snap -> no holes, no spikes.
     *   5. SMOOTH: a light moving-average (ttf.fillSmooth half-window, default 1) on
     *      the positions to remove residual jitter. Normals are re-derived from the
     *      PAE matches (interpolated across fills) so they stay analytic where real.
     *
     * Emits points in image coordinates with aligned unit normals.
     */
    private static void buildPolygonFilled(PartialAreaEdgeDetector.Result paeResult, Wand wand,
            int rx, int ry, FloatPolygon outPoly, List<double[]> outNormals) {
        if (wand.npoints <= 0 || paeResult == null || paeResult.edges == null
                || paeResult.edges.isEmpty()) return;

        double fillTol   = getDoubleProp("ttf.fillTolPx", 6.0);
        double spikeTol  = getDoubleProp("ttf.fillSpikePx", 2.5);
        int smoothHalf   = (int) getDoubleProp("ttf.fillSmooth", 1);

        // 1) densify wand outline (crop coords)
        List<Point> dense = new ArrayList<>();
        for (int i = 0; i < wand.npoints; i++) {
            int x1 = wand.xpoints[i], y1 = wand.ypoints[i];
            int ni = (i + 1) % wand.npoints;
            int x2 = wand.xpoints[ni], y2 = wand.ypoints[ni];
            dense.add(new Point(x1, y1));
            int steps = Math.max(Math.abs(x2 - x1), Math.abs(y2 - y1));
            for (int s = 1; s < steps; s++) {
                dense.add(new Point(x1 + s * (x2 - x1) / steps, y1 + s * (y2 - y1) / steps));
            }
        }
        int n = dense.size();
        if (n < 3) return;

        // centroid of the densified wand (interior reference for the spike guard)
        double cx = 0, cy = 0;
        for (Point p : dense) { cx += p.x; cy += p.y; }
        cx /= n; cy /= n;

        // 2) snap each densified point to nearest PAE edge within the wider tolerance
        double[] px = new double[n], py = new double[n], pnx = new double[n], pny = new double[n];
        boolean[] ok = new boolean[n];
        double[] rad = new double[n];
        for (int i = 0; i < n; i++) {
            double qx = dense.get(i).x, qy = dense.get(i).y;
            PartialAreaEdgeDetector.EdgePixel best = null;
            double minDist = Double.MAX_VALUE;
            for (PartialAreaEdgeDetector.EdgePixel ep : paeResult.edges) {
                double d = Math.hypot(ep.x - qx, ep.y - qy);
                if (d < minDist && d < fillTol) { minDist = d; best = ep; }
            }
            if (best != null) {
                px[i] = best.x; py[i] = best.y; pnx[i] = best.nx; pny[i] = best.ny;
                rad[i] = Math.hypot(best.x - cx, best.y - cy);
                ok[i] = true;
            }
        }

        // 3) spike guard: reject points whose radius deviates from the LOCAL median
        //    (robust to the spike itself). Window is +/- a few points along the order.
        boolean[] spike = new boolean[n];
        int win = 4;
        for (int i = 0; i < n; i++) {
            if (!ok[i]) continue;
            List<Double> local = new ArrayList<>();
            for (int w = -win; w <= win; w++) {
                int j = ((i + w) % n + n) % n;
                if (ok[j]) local.add(rad[j]);
            }
            if (local.size() < 3) continue;
            Collections.sort(local);
            double med = local.get(local.size() / 2);
            if (Math.abs(rad[i] - med) > spikeTol) spike[i] = true;
        }
        for (int i = 0; i < n; i++) if (spike[i]) ok[i] = false;   // demote spikes to gaps

        long holes = 0;
        for (int i = 0; i < n; i++) if (!ok[i]) holes++;
        if (holes == n) return;   // nothing usable

        // 4) fill gaps by circular linear interpolation of position AND normal
        for (int i = 0; i < n; i++) {
            if (ok[i]) continue;
            int prev = i, fwd = i, backSteps = 0, fwdSteps = 0;
            do { prev = (prev - 1 + n) % n; backSteps++; } while (!ok[prev] && backSteps < n);
            do { fwd = (fwd + 1) % n; fwdSteps++; } while (!ok[fwd] && fwdSteps < n);
            double t = (double) backSteps / (backSteps + fwdSteps);
            px[i]  = (1 - t) * px[prev]  + t * px[fwd];
            py[i]  = (1 - t) * py[prev]  + t * py[fwd];
            double ix = (1 - t) * pnx[prev] + t * pnx[fwd];
            double iy = (1 - t) * pny[prev] + t * pny[fwd];
            double len = Math.hypot(ix, iy);
            if (len == 0) { pnx[i] = pnx[prev]; pny[i] = pny[prev]; }
            else { pnx[i] = ix / len; pny[i] = iy / len; }
        }

        // 5) light position smoothing (moving average) along the ordering
        double[] sx = px, sy = py;
        if (smoothHalf > 0) {
            sx = new double[n]; sy = new double[n];
            for (int i = 0; i < n; i++) {
                double ax = 0, ay = 0; int c = 0;
                for (int w = -smoothHalf; w <= smoothHalf; w++) {
                    int j = ((i + w) % n + n) % n;
                    ax += px[j]; ay += py[j]; c++;
                }
                sx[i] = ax / c; sy[i] = ay / c;
            }
        }

        for (int i = 0; i < n; i++) {
            outPoly.addPoint((float) (sx[i] + rx), (float) (sy[i] + ry));
            outNormals.add(new double[]{ pnx[i], pny[i] });
            EdgeDiagnostics.recordAnchor(false);   // hole-free by construction
        }
    }

    /**
     * Edge-polygon construction strategy. Default is "fill" (Strategy B): hole-free
     * wand anchoring with spike guard + smoothing, validated as the most robust
     * across SOFT and LUNG at all noise levels.
     *
     * Selection (first match wins):
     *   -Dttf.edgeStrategy=legacy|normals|fill|primary   explicit override
     *   -Dttf.paePrimary=true   -> "primary"  (Strategy C, back-compat toggle)
     *   -Dttf.paeFill=true      -> "fill"     (Strategy B, back-compat toggle)
     *   -Dttf.paeNormals=true   -> "normals"  (Strategy A, back-compat toggle)
     *   (nothing set)           -> "fill"     (new default)
     *
     *   legacy  : wand snap-or-drop polygon + neighbour-fit normals (pre-change behavior)
     *   normals : wand snap-or-drop polygon + per-point PAE normals (Strategy A)
     *   fill    : hole-free wand polygon (Strategy B) [DEFAULT]
     *   primary : polygon from PAE points directly (Strategy C)
     */
    static String edgeStrategy() {
        String s = System.getProperty("ttf.edgeStrategy");
        if (s != null) {
            s = s.trim().toLowerCase();
            if (s.equals("legacy") || s.equals("normals") || s.equals("fill") || s.equals("primary")) {
                return s;
            }
        }
        if (Boolean.getBoolean("ttf.paePrimary")) return "primary";
        if (Boolean.getBoolean("ttf.paeFill"))    return "fill";
        if (Boolean.getBoolean("ttf.paeNormals"))  return "normals";
        return "fill";   // Strategy B is the default
    }

    private static double getDoubleProp(String key, double def) {
        String v = System.getProperty(key);
        if (v == null) return def;
        try { return Double.parseDouble(v); } catch (NumberFormatException e) { return def; }
    }

    /**
     * Overshoot-robust edge contrast for a single sub-pixel edge pixel.
     *
     * The PAE detector's ep.i0/ep.i1 are taken from the smoothed intensities
     * INSIDE the detector window (+/-4 px around the edge). For edge-enhancing
     * reconstruction kernels (Lung, Bone, Edge, ...) that window straddles the
     * overshoot/undershoot ringing, so i1-i0 overestimates the true step height
     * (e.g. ~122 HU measured vs a 100 HU step on Lung). That inflated value then
     * fails to match the ESF endpoint contrast (measured in the flat regions),
     * so every ESF profile is rejected and the TTF cannot be computed at low noise.
     *
     * This samples the image a standoff distance beyond the ringing, on both
     * sides of the edge along its normal, averaging a few samples for stability.
     * Falls back to the PAE i0/i1 if sampling leaves the image.
     *
     * @param ipRaw  raw image to sample (same coordinate frame as x,y)
     * @param x,y    edge position (image coordinates)
     * @param nx,ny  unit edge normal
     * @return {valLow, valHigh} as RAW intensities (not yet calibrated)
     */
    private static double[] contrastAwayFromEdge(ImageProcessor ipRaw, double x, double y,
                                                 double nx, double ny,
                                                 double fallbackLow, double fallbackHigh) {
        // Standoff past the overshoot lobes; averaged over a short run to reduce noise.
        final double standoff = 5.0;   // px; beyond the typical ringing extent
        final int    nAvg     = 3;     // samples averaged per side
        final double sampleStep = 1.0; // px between averaged samples

        double sumA = 0, sumB = 0;
        int cntA = 0, cntB = 0;
        int w = ipRaw.getWidth(), h = ipRaw.getHeight();
        for (int s = 0; s < nAvg; s++) {
            double d = standoff + s * sampleStep;
            double ax = x + nx * d, ay = y + ny * d;   // one side (outer/inner)
            double bx = x - nx * d, by = y - ny * d;   // the other side
            if (ax >= 0 && ay >= 0 && ax < w && ay < h) { sumA += ipRaw.getInterpolatedValue(ax, ay); cntA++; }
            if (bx >= 0 && by >= 0 && bx < w && by < h) { sumB += ipRaw.getInterpolatedValue(bx, by); cntB++; }
        }
        if (cntA == 0 || cntB == 0) {
            return new double[]{ Math.min(fallbackLow, fallbackHigh), Math.max(fallbackLow, fallbackHigh) };
        }
        double va = sumA / cntA, vb = sumB / cntB;
        return new double[]{ Math.min(va, vb), Math.max(va, vb) };
    }

    private static double getCircularity(FloatPolygon fp) {
        int n = fp.npoints;
        if (n < 3) return 0;
        double area = 0;
        double perimeter = 0;

        for (int i = 0; i < n; i++) {
            int next = (i + 1) % n;
            float x1 = fp.xpoints[i], y1 = fp.ypoints[i];
            float x2 = fp.xpoints[next], y2 = fp.ypoints[next];

            area += (x1 * y2) - (x2 * y1);
            perimeter += Math.hypot(x2 - x1, y2 - y1);
        }
        area = Math.abs(area / 2.0);

        if (perimeter == 0) return 0;
        return 4.0 * Math.PI * (area / (perimeter * perimeter));
    }

    private static ImageProcessor calcularMedia(ImageStack stack) {
        int width = stack.getWidth();
        int height = stack.getHeight();
        FloatProcessor avgSlice = new FloatProcessor(width, height);
        int size = stack.size();

        if (stack.size() <= 0) return stack.getProcessor(1);

        for (int i = 1; i <= size; i++) {
            ImageProcessor ip = stack.getProcessor(i);
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    avgSlice.setf(x, y, avgSlice.getf(x, y) + ip.getf(x, y));
                }
            }
        }

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                avgSlice.setf(x, y, avgSlice.getf(x, y) / size);
            }
        }
        return avgSlice;
    }

    private static boolean containsObject(ByteProcessor mask, int x, int y, double windowSizePx) {
        int endX = Math.min((int)(x + windowSizePx), mask.getWidth());
        int endY = Math.min((int)(y + windowSizePx), mask.getHeight());
        for (int i = x; i < endX; i++) {
            for (int j = y; j < endY; j++) {
                if (mask.getPixel(i, j) == 255) return true;
            }
        }
        return false;
    }

    public static OtsuResult Otsu(int [] data ) {
		int ih;
		int threshold=-1;
		int num_pixels=0;
		double total_mean;
		double bcv, term;
		double max_bcv;
		double [] cnh = new  double [data.length];
		double [] mean = new  double [data.length];
		double [] histo = new  double [data.length];

		for ( ih = 0; ih < data.length; ih++ ) num_pixels += data[ih];
		term = 1.0 / ( double ) num_pixels;

		for ( ih = 0; ih < data.length; ih++ ) histo[ih] = term * data[ih];

		cnh[0] = histo[0];
		for ( ih = 1; ih < data.length; ih++ ) cnh[ih] = cnh[ih - 1] + histo[ih];

		mean[0] = 0.0;
		for ( ih = 1; ih <data.length; ih++ ) mean[ih] = mean[ih - 1] + ih * histo[ih];

		total_mean = mean[data.length-1];

		threshold = 0;
		max_bcv = 0.0;

		for ( ih = 0; ih < data.length; ih++ ) {
			double denom = cnh[ih] * ( 1.0 - cnh[ih] );
			if (denom <= 0.0) continue;

			bcv = total_mean * cnh[ih] - mean[ih];
			bcv *= bcv / denom;

			if ( max_bcv < bcv ) {
				max_bcv = bcv;
				threshold = ih;
			}
		}

	    OtsuResult or = new OtsuResult();
	    or.max_bcv = max_bcv;
	    or.threshold = threshold;
	    return or;
	}

    public static PolygonRoi segmenter(ImageProcessor ip, Calibration cal, double minHU, double minArea, double minCircularity, boolean useManager, int slicenumber, int radius) {
    	PolygonRoi largerRoi = null;
        double largerArea = 0;
        ByteProcessor Mask = new ByteProcessor(ip.getWidth(), ip.getHeight());
        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
            	double value = cal.getCValue(ip.getf(x, y));
                if (value >= minHU) {
                    Mask.putPixel(x, y, 255);
                } else {
                    Mask.putPixel(x, y, 0);
                }
            }
        }
        ByteProcessor filledMask = fillHoles(Mask);

        ByteProcessor erodedMask = filledMask;
        if(radius != 0) {
            erodedMask = erode(filledMask, radius);
        }

        Wand wand = new Wand(erodedMask);

        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                if (erodedMask.getPixel(x, y) == 255) {
                    wand.autoOutline(x, y, 255, 255, 8 | 256);
                    if (wand.npoints > 0) {
                        PolygonRoi erodedRoi = new PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.FREEROI);
                      Rectangle bounds = erodedRoi.getBounds();
                      boolean touchesEdge = (bounds.x <= 0 ||
                                bounds.y <= 0 ||
                                bounds.x + bounds.width >= ip.getWidth() ||
                                bounds.y + bounds.height >= ip.getHeight());

						if (touchesEdge) {
						continue;
						}
						PolygonRoi dilatedRoi = erodedRoi;
                        if(radius != 0) {
                            dilatedRoi = dilate(filledMask, erodedRoi, radius);
                        }

                        ip.setRoi(dilatedRoi);
                        ImageStatistics stats = ImageStatistics.getStatistics(ip, 0, cal);
                        double area = stats.area;
                        double perimeter = dilatedRoi.getLength();
                        double circularity = perimeter == 0.0 ? 0.0 : 4.0 * Math.PI * (stats.pixelCount / (perimeter * perimeter));
                        if (area >= minArea && circularity >= minCircularity) {
                            if (area > largerArea) {
                                largerArea = area;
                                largerRoi = dilatedRoi;
                            }
                        }
                        ip.resetRoi();
                        erodedMask.setValue(0);
                        erodedMask.fill(erodedRoi);
                    }
                }
            }
        }
        if(largerRoi != null) {
            PolygonRoi finalRoi = largerRoi;
            ip.setRoi(finalRoi);
            if (useManager) {
                RoiManager roiManager = RoiManager.getInstance();
                if (roiManager == null) {
                    roiManager = new RoiManager();
                }
                finalRoi.setPosition(slicenumber);
                roiManager.addRoi(finalRoi);
            }
            return finalRoi;
        } else {
            return null;
        }
    }

    public static ByteProcessor erode(ByteProcessor mask, int radius) {
        ByteProcessor erodedMask = new ByteProcessor(mask.getWidth(), mask.getHeight());
        FloatProcessor floatEdm = new EDM().makeFloatEDM(mask, 0, true);
        float[] fPixels = (float[]) floatEdm.getPixels();
        int width = erodedMask.getWidth();
        int height = erodedMask.getHeight();
        for (int y = 0; y < height; y++) {
             for (int x = 0, p = x + y * width; x < width; x++, p++) {
                if (fPixels[p] > radius) {
                    erodedMask.putPixel(x, y, 255);
                } else {
                    erodedMask.putPixel(x, y, 0);
                }
            }
        }
        return erodedMask;
    }

    public static PolygonRoi dilate(ByteProcessor mask, Roi roi, int radius) {
         ByteProcessor dilatedMask = new ByteProcessor(mask.getWidth(), mask.getHeight());
         dilatedMask.setColor(255); dilatedMask.fill();
         dilatedMask.setColor(0); dilatedMask.fill(roi);

         FloatProcessor edm = new EDM().makeFloatEDM(dilatedMask.duplicate(), 0, false);
         float[] edmPixels = (float[]) edm.getPixels();
         byte[] dilatedPixels = (byte[]) dilatedMask.getPixels();
         byte[] originalPixels = (byte[]) mask.getPixels();

         for (int y = 0; y < mask.getHeight(); y++) {
             for (int x = 0; x < mask.getWidth(); x++) {
                 int index = y * mask.getWidth() + x;
                 if (edmPixels[index] <= radius && originalPixels[index] == (byte) 255) {
                     dilatedPixels[index] = (byte) 0;
                 }
             }
         }

         Wand wand = new Wand(dilatedMask);
         for (int y = 0; y < dilatedMask.getHeight(); y++) {
             for (int x = 0; x < dilatedMask.getWidth(); x++) {
                 if (dilatedMask.getPixel(x, y) == 0) {
                     wand.autoOutline(x, y, 0, 0, 8 | 256);
                     if (wand.npoints > 0) {
                          roi = new PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.FREEROI);
                          return (PolygonRoi) roi;
                     }
                 }
             }
         }
         return null;
    }

    public static ByteProcessor fillHoles(ByteProcessor ip) {
        ByteProcessor filledIp = (ByteProcessor) ip.duplicate();

        int width = filledIp.getWidth();
        int height = filledIp.getHeight();
        int background = 0;
        int foreground = 255;
        int tempValue = 127;

        FloodFiller ff = new FloodFiller(filledIp);
        filledIp.setColor(tempValue);

        for (int x = 0; x < width; x++) {
            if (filledIp.getPixel(x, 0) == background) ff.fill(x, 0);
            if (filledIp.getPixel(x, height - 1) == background) ff.fill(x, height - 1);
        }
        for (int y = 0; y < height; y++) {
            if (filledIp.getPixel(0, y) == background) ff.fill(0, y);
            if (filledIp.getPixel(width - 1, y) == background) ff.fill(width - 1, y);
        }
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int pixel = filledIp.getPixel(x, y);
                if (pixel == tempValue) {
                    filledIp.putPixel(x, y, background);
                } else if (pixel == background) {
                    filledIp.putPixel(x, y, foreground);
                }
            }
        }
        return filledIp;
    }

    public static List<Roi> createCircularRecordRois(ImageProcessor ip, Calibration cal, PolygonRoi centerRoi, int recordNumber, double recordSizePixels, double recordRadiusPixels, boolean useManager, int slicenumber) {
        List<Roi> recordRois = new ArrayList<>();

        float[] xPoints = centerRoi.getFloatPolygon().xpoints;
        float[] yPoints = centerRoi.getFloatPolygon().ypoints;
        double centerX = 0, centerY = 0;
        for (int i = 0; i < xPoints.length; i++) {
            centerX += xPoints[i];
            centerY += yPoints[i];
        }
        centerX /= xPoints.length;
        centerY /= yPoints.length;

        for (int i = 0; i < recordNumber; i++) {
            double angle = 2 * Math.PI * i / recordNumber;
            double roiCenterX = centerX + recordRadiusPixels * Math.cos(angle);
            double roiCenterY = centerY + recordRadiusPixels * Math.sin(angle);

            Roi roi = new Roi(
                (int)(roiCenterX - recordSizePixels/2),
                (int)(roiCenterY - recordSizePixels/2),
                (int)recordSizePixels,
                (int)recordSizePixels);

            recordRois.add(roi);
            if (useManager) {
                RoiManager roiManager = RoiManager.getInstance();
                if (roiManager == null) roiManager = new RoiManager();
                roi.setPosition(slicenumber);
                roiManager.addRoi(roi);
            }
        }
        return recordRois;
    }

}

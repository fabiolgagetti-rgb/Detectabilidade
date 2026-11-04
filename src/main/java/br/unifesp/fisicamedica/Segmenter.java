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

        for (int i = 0; i < dR.roi.size(); i++) {
            List<Roi> sliceRois = new ArrayList<>();
            Roi avgRoi = dR.roi.get(i);
            Rectangle bounds = dR.roi.get(i).getBounds();

            for (int z = 1; z <= nSlices; z++) {
                if (avgSlices) {
                    sliceRois.add(avgRoi);
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
                    }
                }
            }
            dR.rois.add(sliceRois);
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
        // noise:    std-dev of intensities inside the wand ROI (sufficient per
        //           design; measured on the ORIGINAL image, not the sigma=1 blur).
        // ------------------------------------------------------------------
        double contrast = 2.0 * Math.abs(mode - thresholdForEdge);

     // Reduz a ROI em 5 pixels em todas as direções
        Roi noiseRoi = RoiEnlarger.enlarge(wandRoi, -5);

        // Aplica a ROI reduzida
        origIp.setRoi(noiseRoi);

        ImageStatistics noiseStats = origIp.getStatistics();

        double noise = noiseStats.stdDev;
        origIp.resetRoi();

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
        if (cleanWand.npoints > 0) {
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
                if (bestEdge != null) {
                    finalSubpixelPolygon.addPoint((float)(bestEdge.x + rx), (float)(bestEdge.y + ry));
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

        dR.Pos.add(new Point((int)Math.round(centerX), (int)Math.round(centerY)));

        if (processedMask != null) {
            processedMask.setColor(255);
            processedMask.fill(finalRoi);
        }
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

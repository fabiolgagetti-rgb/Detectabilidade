package br.unifesp.fisicamedica;

import java.util.EnumMap;
import java.util.Map;

/**
 * Lightweight, opt-in diagnostics channel for understanding WHY processRecord
 * discards a candidate record. Avoids changing the public signatures of
 * Segmenter/TtfCalc (which the ImageJ plugin relies on).
 *
 * Usage in a validation harness:
 *   EdgeDiagnostics.enable();      // start counting
 *   ... run Segmenter.edgeFinder ...
 *   // after the run, inspect EdgeDiagnostics.snapshot() / lastReason()
 *   EdgeDiagnostics.reset();       // clear counters between conditions
 *
 * When disabled (the default), recording is a no-op, so normal plugin use is
 * unaffected and there is no performance cost.
 */
public final class EdgeDiagnostics {

    /** Why a single candidate record was rejected (or accepted). */
    public enum Reason {
        ACCEPTED,              // record passed all gates and produced a ROI
        BCV_GATE,              // Otsu between-class variance below 3000 (weak/edgeless window)
        WAND_NOT_FOUND,        // no seed pixel in [minHU,maxHU] produced a wand outline
        WAND_OFF_ROI,          // wand outline did not touch the record ROI
        BOUNDS_TOUCH_IMAGE,    // wand bounds touched the image border
        LOW_CIRCULARITY,       // final subpixel polygon circularity < 0.3
        CONTRAST_GATE          // expected contrast empty or < 5 HU
    }

    private static boolean enabled = false;
    private static final Map<Reason, Integer> counts = new EnumMap<>(Reason.class);
    private static Reason lastReason = null;

    private EdgeDiagnostics() {}

    public static void enable()  { enabled = true;  }
    public static void disable() { enabled = false; }
    public static boolean isEnabled() { return enabled; }

    public static synchronized void reset() {
        counts.clear();
        lastReason = null;
    }

    /** Called at each decision point in processRecord. No-op unless enabled. */
    public static synchronized void record(Reason reason) {
        if (!enabled) return;
        lastReason = reason;
        counts.merge(reason, 1, Integer::sum);
    }

    public static synchronized Reason lastReason() { return lastReason; }

    public static synchronized Map<Reason, Integer> snapshot() {
        return new EnumMap<>(counts);
    }

    public static synchronized int count(Reason r) {
        return counts.getOrDefault(r, 0);
    }

    /** One-line summary of all non-zero reason counts, e.g. "BCV_GATE=18 LOW_CIRCULARITY=4". */
    public static synchronized String summary() {
        StringBuilder sb = new StringBuilder();
        for (Reason r : Reason.values()) {
            int c = counts.getOrDefault(r, 0);
            if (c > 0) {
                if (sb.length() > 0) sb.append(' ');
                sb.append(r).append('=').append(c);
            }
        }
        return sb.length() == 0 ? "(none)" : sb.toString();
    }
}

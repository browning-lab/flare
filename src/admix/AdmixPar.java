/*
 * Copyright 2021-2023 Brian L. Browning
 *
 * This file is part of the flare program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package admix;

import blbutil.Const;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code AdmixPar} represents the command line parameters for the
 * admix program.</p>
 *
 * <p>Instances of class {@code AdmixPar} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AdmixPar {

    private final String[] args;

    // Required parameters
    private final File ref;
    private final File ref_panel;
    private final File gt;
    private final File map;
    private final String out;

    // Optional parameters
    private final boolean array;
    private final float min_maf;
    private final int min_mac;
    private final boolean probs;
    private final float gen;
    private final File model;
    private final boolean em;
    private final boolean update_p;
    private final String gt_samples;
    private final File gtSamplesFile;
    private final File gt_ancestries;
    private final File excludemarkers;
    private final boolean panel_probs;
    private final float panel_cm;
    private final int nthreads;
    private final long seed;

    // Undocumented parameters
    private final float panel_weight;
    private final int states;
    private final float ibs_step;
    private final float ibs_buffer;
    private final int ibs_haps;
    private final float ibs_recycle;
    private final int em_its;
    private final int em_haps;
    private final float delta_mu;
    private final float delta_p;
    private final int panel_markers;
    private final int panel_haps;
    private final int panel_ne;
    private final boolean debug;

    private static final boolean DEF_ARRAY = false;
    private static final float DEF_MIN_MAF = 0.005f;
    private static final int DEF_MIN_MAC = 50;
    private static final boolean DEF_PROBS = false;
    static final float DEF_GEN = 10f;
    private static final boolean DEF_EM = true;
    static final boolean DEF_UPDATE_P = false;
    static final boolean DEF_PANEL_PROBS = false;
    static final float DEF_PANEL_CM = 0.5f;
    private static final int DEF_THREADS = Runtime.getRuntime().availableProcessors();
    private static final long DEF_SEED = -99999;

    static final float DEF_PANEL_WEIGHT = 0.80f;
    static final int DEF_STATES = 100;
    static final float DEF_IBS_STEP = 0.01f;
    static final float DEF_IBS_BUFFER = 2.0f;
    static final int DEF_IBS_HAPS = 4;
    static final float DEF_IBS_RECYCLE = 4.0f;
    static final int DEF_EM_ITS = 20;
    static final int DEF_EM_HAPS = 100;
    static final float DEF_DELTA_MU = 0.03f;
    static final float DEF_DELTA_P = 0.01f;
    static final int DEF_PANEL_MARKERS = 10;
    static final int DEF_PANEL_HAPS = 1000;
    static final int DEF_PANEL_NE = 100000;
    static final boolean DEF_DEBUG = false;

    /**
     * Constructs an {@code AdmixPar} instance that represents the
     * command line parameters.  See the {@code usage()} method for a
     * description of the command line parameters. The constructor exits with
     * an error message if a command line parameter name is not recognized.
     *
     * @param args the command line parameters
     * @throws IllegalArgumentException if a command line parameter is
     * incorrectly specified
     * @throws NullPointerException if {@code args == null}
    */
    public AdmixPar(String[] args) {
        int IMAX = Integer.MAX_VALUE;
        long LMIN = Long.MIN_VALUE;
        long LMAX = Long.MAX_VALUE;
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;
        this.args = args.clone();
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        // Required parameters
        ref = Validate.getFile(Validate.stringArg("ref", argsMap, true, null,
                null));
        ref_panel = Validate.getFile(Validate.stringArg("ref-panel", argsMap, true,
                null, null));
        gt = Validate.getFile(Validate.stringArg("gt", argsMap, true, null,
                null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, true, null, null));
        out = Validate.stringArg("out", argsMap, true, null, null);

        // Optional parameters
        array = Validate.booleanArg("array", argsMap, false, DEF_ARRAY);
        min_maf = Validate.floatArg("min-maf", argsMap, false, DEF_MIN_MAF, -FMAX, Math.nextDown(0.5f));
        min_mac = Validate.intArg("min-mac", argsMap, false, DEF_MIN_MAC, 0, IMAX);
        probs = Validate.booleanArg("probs", argsMap, false, DEF_PROBS);
        gen = Validate.floatArg("gen", argsMap, false, DEF_GEN, 1, IMAX);
        model = Validate.getFile(Validate.stringArg("model", argsMap, false,
                null, null));
        em = Validate.booleanArg("em", argsMap, false, DEF_EM);
        update_p = Validate.booleanArg("update-p", argsMap, false, DEF_UPDATE_P);
        gt_samples = Validate.stringArg("gt-samples", argsMap, false, null, null);
        gtSamplesFile = AdmixPar.gtSamplesFile(gt_samples);
        gt_ancestries = Validate.getFile(
                Validate.stringArg("gt-ancestries", argsMap, false, null, null));
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));
        panel_probs = Validate.booleanArg("panel-probs", argsMap, false, DEF_PANEL_PROBS);
        panel_cm = Validate.floatArg("panel-cm", argsMap, false, DEF_PANEL_CM, FMIN, FMAX);
        nthreads = Validate.intArg("nthreads", argsMap, false, DEF_THREADS, 1, IMAX);
        seed = Validate.longArg("seed", argsMap, false, DEF_SEED, LMIN, LMAX);

        // Undocumented parameters
        panel_weight = Validate.floatArg("panel-weight", argsMap, false, DEF_PANEL_WEIGHT, FMIN, 1.0f);
        states = Validate.intArg("states", argsMap, false, DEF_STATES, 1, IMAX);
        ibs_step = Validate.floatArg("ibs-step", argsMap, false, DEF_IBS_STEP, FMIN, FMAX);
        ibs_buffer = Validate.floatArg("ibs-buffer", argsMap, false, DEF_IBS_BUFFER, FMIN, FMAX);
        ibs_haps = Validate.intArg("ibs-haps", argsMap, false, DEF_IBS_HAPS, 1, IMAX);
        ibs_recycle = Validate.floatArg("ibs-recycle", argsMap, false, DEF_IBS_RECYCLE, FMIN, FMAX);
        em_its = Validate.intArg("em-its", argsMap, false, DEF_EM_ITS, 0, IMAX);
        em_haps = Validate.intArg("em-haps", argsMap, false, DEF_EM_HAPS, 1, IMAX);
        delta_mu = Validate.floatArg("delta-mu", argsMap, false, DEF_DELTA_MU, 0.0f, 1.0f);
        delta_p = Validate.floatArg("delta-p", argsMap, false, DEF_DELTA_P, 0.0f, 1.0f);
        panel_markers = Validate.intArg("panel-markers", argsMap, false, DEF_PANEL_MARKERS, 2, IMAX);
        panel_haps = Validate.intArg("panel-haps", argsMap, false, DEF_PANEL_HAPS, 1, IMAX);
        panel_ne = Validate.intArg("panel-ne", argsMap, false, DEF_PANEL_NE, 1, IMAX);
        debug = Validate.booleanArg("debug", argsMap, false, DEF_DEBUG);
        Validate.confirmEmptyMap(argsMap);
    }

    private static File gtSamplesFile(String gtSamples) {
        if (gtSamples==null) {
            return null;
        }
        else if (gtSamples.startsWith("^")) {
            return Validate.getFile(gtSamples.substring(1));
        } else {
            return Validate.getFile(gtSamples);
        }
    }

    /**
     * Returns the command line arguments.
     * @return the command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns a string representation of the flare command and version.  The
     * exact details of the representation are unspecified and subject to change.
     * @return a string representation of the flare command and version
     */
    public static String flareCommand() {
        String commandLine = ProcessHandle.current().info().commandLine().orElse("");
        StringBuilder sb = new StringBuilder(commandLine);
        sb.append("  # ");
        sb.append(AdmixMain.VERSION);
        return sb.toString();
    }

    /**
     * Returns a string describing the command line arguments.
     * The format of the returned string is unspecified and subject to change.
     * @return a string describing the command line arguments.
     */
    public static String usage() {
        String nl = Const.nl;
        return "Syntax: " + AdmixMain.COMMAND + " [arguments in format: parameter=value]" + nl
                + nl
                + "Required Parameters:" + nl
                + "  ref=<VCF file with phased reference genotypes>        (required)" + nl
                + "  ref-panel=<file with reference sample to panel map>   (required)" + nl
                + "  gt=<VCF file with phased genotypes to be analyzed>    (required)" + nl
                + "  map=<PLINK map file with cM units>                    (required)" + nl
                + "  out=<output file prefix>                              (required)" + nl + nl

                + "Optional Parameters:" + nl
                + "  array=<genotypes are from a SNP array: true/false>    (default: " + DEF_ARRAY + ")" + nl
                + "  min-maf=<minimum MAF in reference VCF file>           (default: " + DEF_MIN_MAF + ")" + nl
                + "  min-mac=<minimum MAC in reference VCF file>           (default: " + DEF_MIN_MAC + ")" + nl
                + "  probs=<report ancestry probs: true/false>             (default: " + DEF_PROBS + ")" + nl
                + "  gen=<number of generations since admixture>           (default: " + DEF_GEN + ")" + nl
                + "  model=<file with model parameters>                    (optional)" + nl
                + "  em=<estimate model parameters using EM: true/false>   (default: " + DEF_EM + ")" + nl
                + "  update-p=<estimate p and rho using EM: true/false>    (default: " + DEF_UPDATE_P + ")" + nl
                + "  gt-samples=<file with sample IDs to analyze>          (optional)" + nl
                + "  gt-ancestries=<file with sample ancestry proportions> (optional)" + nl
                + "  excludemarkers=<file with markers to exclude>         (optional)" + nl
                + "  panel-probs=<estimate panel probs: true/false>        (default: " + DEF_PANEL_PROBS + ")" + nl
                + "  panel-cm=<window size if panel-probs=true>            (default: " + DEF_PANEL_CM + ")" + nl
                + "  nthreads=<number of computational threads>            (default: all CPU cores)" + nl
                + "  seed=<seed for random number generations>             (default: " + DEF_SEED + ")" + nl;
    }

    // Required parameters

    /**
     * Returns the ref file.
     * @return the ref file
     */
    public File ref() {
        return ref;
    }

    /**
     * Returns the ref-panel file.
     * @return the ref-panel file
     */
    public File ref_panel() {
        return ref_panel;
    }

    /**
     * Returns the gt file.
     * @return the gt file
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the map file.
     * @return the map file
     */
    public File map() {
        return map;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    // Optional parameters

    /**
     * Returns the array parameter.
     * @return the array parameter
     */
    public boolean array() {
        return array;
    }

    /**
     * Returns the min-maf parameter.
     * @return the min-maf parameter
     */
    public float min_maf() {
        return min_maf;
    }

    /**
     * Returns the min-mac parameter.
     * @return the min-mac parameter
     */
    public int min_mac() {
        return min_mac;
    }

    /**
     * Returns the probs parameter.
     * @return the probs parameter
     */
    public boolean probs() {
        return probs;
    }

    /**
     * Returns the gen parameter.
     * @return the gen parameter
     */
    public float gen() {
        return gen;
    }

    /**
     * Returns the model file or {@code null} if no model parameter was
     * specified.
     * @return the model file
     */
    public File model() {
        return model;
    }

    /**
     * Returns the em parameter.
     * @return the em parameter
     */
    public boolean em() {
        return em;
    }

    /**
     * Returns the update-p parameter.
     * @return the update-p parameter
     */
    public boolean update_p() {
        return update_p;
    }

    /**
     * Returns the gt-samples parameter or {@code null}
     * if no gt-samples parameter was specified.
     *
     * @return the gt-samples parameter or {@code null}
     * if no gt-samples parameter was specified.
     */
    public String gt_samples() {
        return gt_samples;
    }

    /**
     * Returns the file specified with the gt-samples parameter or
     * {@code null} if no gt-samples parameter was specified.
     *
     * @return the file specified with the gt-samples parameter or
     * {@code null} if no gt-samples parameter was specified.
     */
    public File gtSamplesFile() {
        return gtSamplesFile;
    }

    /**
     * Returns the gt-ancestries file or {@code null}
     * if no gt-ancestries parameter was specified.
     *
     * @return the gt-ancestries file or {@code null}
     * if no gt-ancestries parameter was specified.
     */
    public File gt_ancestries() {
        return gt_ancestries;
    }

    /**
     * Returns the excludemarkers file or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the exclude-markers file or {@code null}
     * if no exclude-markers parameter was specified
     */
    public File excludemarkers() {
        return excludemarkers;
    }

    /**
     * Returns the panel-probs parameter
     * @return the panel-probs parameter
     */
    public boolean panel_probs() {
        return panel_probs;
    }

    /**
     * Returns the panel-cm parameter.
     * @return the panel-cm parameter
     */
    public float panel_cm() {
        return panel_cm;
    }

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter
     */
    public int nthreads() {
        return nthreads;
    }

    /**
     * Returns the seed parameter.
     * @return the seed parameter
     */
    public long seed() {
        return seed;
    }

    // Undocumented parameters

    /**
     * Returns the panel-weight parameter.
     * @return the panel-weight parameter
     */
    public float panel_weight() {
        return panel_weight;
    }

    /**
     * Returns the states parameter.
     * @return the states parameter
     */
    public int states() {
        return states;
    }

    /**
     * Returns the ibs-step parameter.
     * @return the ibs-step parameter
     */
    public float ibs_step() {
        return ibs_step;
    }

    /**
     * Returns the ibs-buffer parameter.
     * @return the ibs-buffer parameter
     */
    public float ibs_buffer() {
        return ibs_buffer;
    }

    /**
     * Returns the ibs-haps parameter.
     * @return the ibs-haps parameter
     */
    public int ibs_haps() {
        return ibs_haps;
    }

    /**
     * Returns the ibs-recycle parameter.
     * @return the ibs-recycle parameter
     */
    public float ibs_recycle() {
        return ibs_recycle;
    }

    /**
     * Returns the em-its parameter.
     * @return the em-its parameter
     */
    public int em_its() {
        return em_its;
    }

    /**
     * Returns the em-haps parameter.
     * @return the em-haps parameter
     */
    public int em_haps() {
        return em_haps;
    }

    /**
     * Return the delta-mu parameter
     * @return the delta-mu parameter
     */
    public float delta_mu() {
        return delta_mu;
    }

    /**
     * Return the delta-p parameter
     * @return the delta-p parameter
     */
    public float delta_p() {
        return delta_p;
    }

    /**
     * Returns the panel-markers parameter.
     * @return the panel-markers parameter
     */
    public int panel_markers() {
        return panel_markers;
    }

    /**
     * Return the panel-haps parameter
     * @return the panel-haps parameter
     */
    public int panel_haps() {
        return panel_haps;
    }

    /**
     * Return the panel-ne parameter
     * @return the panel-ne parameter
     */
    public int panel_ne() {
        return panel_ne;
    }

    /**
     * Returns the debug parameter.
     * @return the debug parameter
     */
    public boolean debug() {
        return debug;
    }
}

/*
 * Copyright 2021 Brian L. Browning
 * Copyright 2023 Genomics plc
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
import blbutil.Utilities;
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
    private final File anc_panel;
    private final boolean array;
    private final float min_maf;
    private final int min_mac;
    private final boolean probs;
    private final float gen;
    private final File model;
    private final boolean em;
    private final int nthreads;
    private final long seed;
    private final String gt_samples;
    private final File gtSamplesFile;
    private final File excludemarkers;

    // Undocumented parameters
    private final float panel_weight;
    private final int states;
    private final float ibs_step;
    private final float ibs_buffer;
    private final int ibs_haps;
    private final float ibs_recycle;
    private final int em_its;
    private final int em_haps;
    private final float em_anc_prob;
    private final float delta_mu;
    private final float min_mu;
    private final boolean debug;

    private static final boolean DEF_ARRAY = false;
    private static final float DEF_MIN_MAF = 0.005f;
    private static final int DEF_MIN_MAC = 50;
    private static final boolean DEF_PROBS = false;
    static final float DEF_GEN = 10f;
    private static final boolean DEF_EM = true;
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
    static final float DEF_EM_ANC_PROB = 0.9f;
    static final float DEF_DELTA_MU = 0.05f;
    static final float DEF_MIN_MU = 0.001f;
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
        anc_panel = null; // anc_panel is not currently a valid parameter
//      anc_panel = Validate.getFile(Validate.stringArg("anc-panel", argsMap, false,
//              null, null));
        array = Validate.booleanArg("array", argsMap, false, DEF_ARRAY);
        min_maf = Validate.floatArg("min-maf", argsMap, false, DEF_MIN_MAF, -FMAX, Math.nextDown(0.5f));
        min_mac = Validate.intArg("min-mac", argsMap, false, DEF_MIN_MAC, 0, IMAX);
        probs = Validate.booleanArg("probs", argsMap, false, DEF_PROBS);
        gen = Validate.floatArg("gen", argsMap, false, DEF_GEN, 1, IMAX);
        model = Validate.getFile(Validate.stringArg("model", argsMap, false,
                null, null));
        em = Validate.booleanArg("em", argsMap, false, DEF_EM);
        nthreads = Validate.intArg("nthreads", argsMap, false, DEF_THREADS, 1, IMAX);
        seed = Validate.longArg("seed", argsMap, false, DEF_SEED, LMIN, LMAX);
        gt_samples = Validate.stringArg("gt-samples", argsMap, false, null, null);
        gtSamplesFile = AdmixPar.gtSamplesFile(gt_samples);
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));

        // Undocumented parameters
        panel_weight = Validate.floatArg("panel-weight", argsMap, false, DEF_PANEL_WEIGHT, FMIN, 1.0f);
        states = Validate.intArg("states", argsMap, false, DEF_STATES, 1, IMAX);
        ibs_step = Validate.floatArg("ibs-step", argsMap, false, DEF_IBS_STEP, FMIN, FMAX);
        ibs_buffer = Validate.floatArg("ibs-buffer", argsMap, false, DEF_IBS_BUFFER, FMIN, FMAX);
        ibs_haps = Validate.intArg("ibs-haps", argsMap, false, DEF_IBS_HAPS, 1, IMAX);
        ibs_recycle = Validate.floatArg("ibs-recycle", argsMap, false, DEF_IBS_RECYCLE, FMIN, FMAX);
        em_its = Validate.intArg("em-its", argsMap, false, DEF_EM_ITS, 0, IMAX);
        em_haps = Validate.intArg("em-haps", argsMap, false, DEF_EM_HAPS, 1, IMAX);
        em_anc_prob = Validate.floatArg("em-anc-prob", argsMap, false, DEF_EM_ANC_PROB, 0.0f, 1.0f);
        delta_mu = Validate.floatArg("delta-mu", argsMap, false, DEF_DELTA_MU, 0.0f, 1.0f);
        min_mu = Validate.floatArg("min-mu", argsMap, false, DEF_MIN_MU, 0.0f, 1.0f);
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
     * Prints an error message and terminates the Java virtual machine if
     * the specified file is an input file.
     * @param filename a filename
     * @throws NullPointerException if {@code filename == null}
     */
    void verifyNotAnInputFile(String filename) {
        File file = new File(filename);
        if (file.equals(gt) || file.equals(ref) || file.equals(map)
                || file.equals(gtSamplesFile) || file.equals(excludemarkers)
                || file.equals(ref_panel) || file.equals(anc_panel)) {
            String err = "An output file has the same name as an input file";
            String info = Const.nl + "Error      :  " + err
                    + Const.nl     + "Filename   :  " + file;
            Utilities.exit(new Throwable(err), info);
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
     * Returns a string describing the command line arguments.
     * The format of the returned string is unspecified and subject to change.
     * @return a string describing the command line arguments.
     */
    public static String usage() {
        String nl = Const.nl;
        return "Syntax: " + AdmixMain.COMMAND + " [arguments in format: parameter=value]" + nl
                + nl
                + "Required Parameters:" + nl
                + "  ref=<VCF/BREF3 file with phased reference genotypes> (required)" + nl
                + "  ref-panel=<file with reference sample to panel map>  (required)" + nl
                + "  gt=<VCF file with phased genotypes to be analyzed>   (required)" + nl
                + "  map=<PLINK map file with cM units>                   (required)" + nl
                + "  out=<output file prefix>                             (required)" + nl + nl

                + "Optional Parameters:" + nl
//                + "  anc-panel=<file with ancestry to panels map>         (optional)" + nl
                + "  array=<genotypes are from a SNP array: true/false>   (default: " + DEF_ARRAY + ")" + nl
                + "  min-maf=<minimum MAF in reference file>              (default: " + DEF_MIN_MAF + ")" + nl
                + "  min-mac=<minimum MAC in reference file>              (default: " + DEF_MIN_MAC + ")" + nl
                + "  probs=<report ancestry probs: true/false>            (default: " + DEF_PROBS + nl
                + "  gen=<number of generations since admixture>          (default: " + DEF_GEN + ")" + nl
                + "  model=<file with model parameters>                   (optional)" + nl
                + "  em=<estimate model parameters using EM: true/false>  (default: " + DEF_EM + ")" + nl
                + "  nthreads=<number of computational threads>           (default: all CPU cores)" + nl
                + "  seed=<seed for random number generations>            (default: " + DEF_SEED + ")" + nl
                + "  gt-samples=<file with sample IDs to analyze>         (optional)" + nl
                + "  excludemarkers=<file with markers to exclude>        (optional)" + nl;
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
     * Returns the anc-panel file or {@code null} if no anc-panel parameter
     * was specified.
     * @return the anc-panel file
     */
    public File anc_panel() {
        return anc_panel;
    }

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
     * Returns the excludemarkers file or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the exclude-markers file or {@code null}
     * if no exclude-markers parameter was specified
     */
    public File excludemarkers() {
        return excludemarkers;
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
     * Returns the em-anc-prob parameter.
     * @return the em-anc-prob parameter
     */
    public float em_anc_prob() {
        return em_anc_prob;
    }

    /**
     * Return the delta-mu parameter
     * @return the delta-mu parameter
     */
    public float delta_mu() {
        return delta_mu;
    }

    /**
     * Return the min-mu parameter
     * @return the min-mu parameter
     */
    public float min_mu() {
        return min_mu;
    }

    /**
     * Returns the debug parameter.
     * @return the debug parameter
     */
    public boolean debug() {
        return debug;
    }
}

/*
 * Copyright 2021 Brian L. Browning
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
import blbutil.FileUtil;
import blbutil.Utilities;
import java.io.File;
import java.io.PrintWriter;
import java.util.Locale;
import java.util.Optional;

/**
 * <p>Class {@code AdmixMain} contains the main() method entry point for the
 * admix program.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixMain {

    static final String EXECUTABLE = "flare.jar";
    static final String REVISION = "flare.20Oct22.c2c.jar";
    static final String PROGRAM = EXECUTABLE + "  [ version 0.3.0, 20Oct22.c2c ]";
    static final String COPYRIGHT = "Copyright (C) 2022 Brian L. Browning";
    static final String COMMAND = "java -jar " + EXECUTABLE;

    private static final String HELP_MESSAGE = "Enter \"" + COMMAND
            + "\" to print a list of command line arguments";

    private AdmixMain() {
        // private constructor to prevent instantiation
    }

    /**
     * Entry point to the admix program.  See the {@code admix.AdmixPar}
     * class for details of the command line arguments.
     * @param args the command line arguments
     * @throws NullPointerException if {@code args == null}
     */
    public static void main(String[] args) {
	Locale.setDefault(Locale.US);
        AdmixPar par = getPar(args);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",
                String.valueOf(par.nthreads()));

        File logOut = AdmixUtils.clobberAndReturnFile(par.out() + ".log");
        try (PrintWriter log = FileUtil.printWriter(logOut)) {
            long t0 = System.nanoTime();
            Utilities.duoPrintln(log, startInfo(par));
            if (par.model()!=null) {
                printWarnings(par, log);
            }
            runAnalysis(par, log);
            Utilities.duoPrintln(log, endInfo(t0));
        }
        catch (Throwable t) {
            // Forces JVM exit if error if file-reading daemon thread exists
            Utilities.exit(t);
        }
    }

    private static void printWarnings(AdmixPar par, PrintWriter log) {
        boolean printWarnings = par.anc_panel()!=null
                || par.gen()!=AdmixPar.DEF_GEN
                || par.array()==true;
        if (printWarnings) {
             Utilities.duoPrintln(log, "");
        }
        if (par.anc_panel()!=null) {
            printAncPanelWarning(log);
        }
        if (par.gen()!=AdmixPar.DEF_GEN) {
            printGenWarning(log);
        }
        if (par.array()) {
            printMinMacNote(log);
        }
        if (printWarnings) {
            Utilities.duoPrintln(log, "");
        }
    }

    private static void printAncPanelWarning(PrintWriter log) {
        String msg = "Warning: Ignoring anc-panel file since a model file is specified";
        Utilities.duoPrintln(log, msg);
    }

    private static void printGenWarning(PrintWriter log) {
        String msg = "Warning: Ignoring 'gen' parameter since a model file is specified";
        Utilities.duoPrintln(log, msg);
    }

    private static void printMinMacNote(PrintWriter log) {
        String msg = "Note: the minor allele count filter is not applied when 'array=true'";
        Utilities.duoPrintln(log, msg);
    }

    private static void runAnalysis(AdmixPar par, PrintWriter log) {
        try (AdmixChromData.It chromIt = new AdmixChromData.It(par)) {
            FixedParams fixedParams = FixedParams.create(par,
                    chromIt.refSamples(), chromIt.targSamples());

            try (AdmixWriter admixWriter = new AdmixWriter(fixedParams)) {
                Optional<AdmixChromData> optChromData = chromIt.nextChrom();
                if (optChromData.isPresent()==false) {
                    missingDataError(par);
                }
                AdmixChromData chromData = optChromData.get();
                boolean includeRefHaps = oneAncPerPanel(fixedParams);
                SelectedHaps selectedHaps = new SelectedHaps(fixedParams, includeRefHaps);
                IbsHaps ibsHaps = new IbsHaps(chromData, selectedHaps);

                ParamsInterface params = ParamEstimator.getParams(fixedParams, ibsHaps, log);
                writeModelFile(params);

                AncEstimator.estAncestry(ibsHaps, params, admixWriter);
                optChromData = chromIt.nextChrom();
                while (optChromData.isPresent()) {
                    chromData = optChromData.get();
                    selectedHaps = selectedHaps.removeRefHaps();
                    ibsHaps = new IbsHaps(chromData, selectedHaps);
                    AncEstimator.estAncestry(ibsHaps, params, admixWriter);
                    optChromData = chromIt.nextChrom();
                }
            }
            Utilities.duoPrint(log, statistics(chromIt));
        }
    }

    private static boolean oneAncPerPanel(FixedParams fixedParams) {
        int nAnc = fixedParams.nAnc();
        if (nAnc!=fixedParams.nRefPanels()) {
            return false;
        }
        for (int i=0; i<nAnc; ++i) {
            if (fixedParams.ancPanels(i).size()!=1) {
                return false;
            }
        }
        return true;
    }

    private static void writeModelFile(ParamsInterface params) {
        AdmixPar par = params.fixedParams().par();
        String filename = par.out() + ".model";
        par.verifyNotAnInputFile(filename);
        File modelFile =  new File(filename);
        try (PrintWriter out = FileUtil.printWriter(modelFile)) {
            out.print(params.toString());
        }
    }

    private static AdmixPar getPar(String[] args) {
        if (args.length==0 || args[0].toLowerCase().startsWith("help")) {
            System.out.println(PROGRAM);
            System.out.println();
            System.out.println(AdmixPar.usage());
            System.exit(0);
        }
        AdmixPar par = new AdmixPar(args);
        checkOutputPrefix(par);
        par.verifyNotAnInputFile(par.out() + ".anc.vcf.gz");
        par.verifyNotAnInputFile(par.out() + ".log");
        return par;
    }

    private static void checkOutputPrefix(AdmixPar par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String err = "The out parameter cannot be a directory";
            String info = Const.nl + "Error      :  " + err
                    + Const.nl     + "Parameter  :  " + "out=" + par.out()
                    + Const.nl     + "Directory  :  " + outPrefix;
            Utilities.exit(new Throwable(err), info);
        }
    }

    private static void missingDataError(AdmixPar par) {
        String err = "Missing genotype data: ";
        String info = Const.nl + "Error          :  " + err
                + Const.nl     + "gt parameter   :  " + par.gt()
                + Const.nl     + "ref parameter  :  " + par.ref();
        Utilities.exit(new Throwable(err), info);
    }

    private static String startInfo(AdmixPar par) {
        StringBuilder sb = new StringBuilder(300);
        sb.append(COPYRIGHT);
        sb.append(Const.nl);
        sb.append(HELP_MESSAGE);
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("Program             :  ");
        sb.append(PROGRAM);
        sb.append(Const.nl);
        sb.append("Start Time          :  ");
        sb.append(Utilities.timeStamp());
        sb.append(Const.nl);
        sb.append("Max Memory          :  ");
        long maxMemory = Runtime.getRuntime().maxMemory();
        if (maxMemory != Long.MAX_VALUE) {
            long maxMb = maxMemory / (1024*1024);
            sb.append(maxMb);
            sb.append(" MB");
        }
        else {
            sb.append("[no limit])");
        }
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append(parameters(par));
        return sb.toString();
    }

    private static String parameters(AdmixPar par) {
        StringBuilder sb = new StringBuilder(150);
        sb.append("Parameters");
        addRequiredParams(par, sb);
        addOptionalParams(par, sb);
        addNonDefaultUndocumentedParams(par, sb);
        sb.append(Const.nl);
        return sb.toString();
    }

    private static void addRequiredParams(AdmixPar par, StringBuilder sb) {
        sb.append(Const.nl);
        sb.append("  ref               :  ");
        sb.append(par.ref());
        sb.append(Const.nl);
        sb.append("  ref-panel         :  ");
        sb.append(par.ref_panel());
        sb.append(Const.nl);
        sb.append("  gt                :  ");
        sb.append(par.gt());
        sb.append(Const.nl);
        sb.append("  map               :  ");
        sb.append(par.map());
        sb.append(Const.nl);
        sb.append("  out               :  ");
        sb.append(par.out());
    }

    private static void addOptionalParams(AdmixPar par, StringBuilder sb) {
        if (par.anc_panel()!=null) {
            sb.append(Const.nl);
            sb.append("  anc-panel         :  ");
            sb.append(par.anc_panel());
        }
        sb.append(Const.nl);
        sb.append("  array             :  ");
        sb.append(par.array());
        sb.append(Const.nl);
        sb.append("  min-maf           :  ");
        sb.append(par.min_maf());
        sb.append(Const.nl);
        if (par.array()==false) {
            sb.append("  min-mac           :  ");
            sb.append(par.min_mac());
            sb.append(Const.nl);
        }
        sb.append("  probs             :  ");
        sb.append(par.probs());
        sb.append(Const.nl);
        if (par.model()==null) {
            sb.append("  gen               :  ");
            sb.append(par.gen());
        }
        else {
            sb.append("  model             :  ");
            sb.append(par.model());
        }
        sb.append(Const.nl);
        sb.append("  em                :  ");
        sb.append(par.em());
        sb.append(Const.nl);
        sb.append("  nthreads          :  ");
        sb.append(par.nthreads());
        sb.append(Const.nl);
        sb.append("  seed              :  ");
        sb.append(par.seed());
        if (par.gt_samples()!=null) {
            sb.append(Const.nl);
            sb.append("  gt-samples        :  ");
            sb.append(par.gt_samples());
        }
        if (par.excludemarkers()!=null) {
            sb.append(Const.nl);
            sb.append("  excludemarkers    :  ");
            sb.append(par.excludemarkers());
        }
    }

    private static void addNonDefaultUndocumentedParams(AdmixPar par,
            StringBuilder sb) {
        if (par.panel_weight()!=AdmixPar.DEF_PANEL_WEIGHT) {
            sb.append(Const.nl);
            sb.append("  panel-weight      :  ");
            sb.append(par.panel_weight());
        }
        if (par.states()!=AdmixPar.DEF_STATES) {
            sb.append(Const.nl);
            sb.append("  states            :  ");
            sb.append(par.states());
        }
        if (par.ibs_step()!=AdmixPar.DEF_IBS_STEP) {
            sb.append(Const.nl);
            sb.append("  ibs-step          :  ");
            sb.append(par.ibs_step());
        }
        if (par.ibs_buffer()!=AdmixPar.DEF_IBS_BUFFER) {
            sb.append(Const.nl);
            sb.append("  ibs-buffer        :  ");
            sb.append(par.ibs_buffer());
        }
        if (par.ibs_haps()!=AdmixPar.DEF_IBS_HAPS) {
            sb.append(Const.nl);
            sb.append("  ibs-haps          :  ");
            sb.append(par.ibs_haps());
        }
        if (par.ibs_recycle()!=AdmixPar.DEF_IBS_RECYCLE) {
            sb.append(Const.nl);
            sb.append("  ibs-recycle       :  ");
            sb.append(par.ibs_recycle());
        }
        if (par.em_its()!=AdmixPar.DEF_EM_ITS) {
            sb.append(Const.nl);
            sb.append("  em-its            :  ");
            sb.append(par.em_its());
        }
        if (par.em_haps()!=AdmixPar.DEF_EM_HAPS) {
            sb.append(Const.nl);
            sb.append("  em-haps           :  ");
            sb.append(par.em_haps());
        }
        if (par.em_anc_prob()!=AdmixPar.DEF_EM_ANC_PROB) {
            sb.append(Const.nl);
            sb.append("  em-anc-prob       :  ");
            sb.append(par.em_its());
        }
        if (par.delta_mu()!=AdmixPar.DEF_DELTA_MU) {
            sb.append(Const.nl);
            sb.append("  delta-mu          :  ");
            sb.append(par.delta_mu());
        }
        if (par.min_mu()!=AdmixPar.DEF_MIN_MU) {
            sb.append(Const.nl);
            sb.append("  min-mu            :  ");
            sb.append(par.min_mu());
        }
        if (par.debug()!=AdmixPar.DEF_DEBUG) {
            sb.append(Const.nl);
            sb.append("  debug             :  ");
            sb.append(par.debug());
        }
    }

    private static String statistics(AdmixChromData.It dataIt) {
        StringBuilder sb = new StringBuilder(300);
        sb.append("Statistics");
        sb.append(Const.nl);
        sb.append("  reference samples :  ");
        sb.append(dataIt.refSamples().size());
        sb.append(Const.nl);
        sb.append("  target samples    :  ");
        sb.append(dataIt.targSamples().size());
        sb.append(Const.nl);
        sb.append("  markers           :  ");
        sb.append(dataIt.nMarkers());
        sb.append(Const.nl);
        return sb.toString();
    }

    private static String endInfo(long startNanoTime) {
        StringBuilder sb = new StringBuilder(300);
        long elapsedNanoTime = System.nanoTime() - startNanoTime;
        sb.append(Const.nl);
        sb.append("Wallclock Time      :  ");
        sb.append(Utilities.elapsedNanos(elapsedNanoTime));
        sb.append(Const.nl);
        sb.append("End Time            :  ");
        sb.append(Utilities.timeStamp());
        return sb.toString();
    }
}

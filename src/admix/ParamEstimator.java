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
import blbutil.MultiThreadUtils;
import blbutil.Utilities;
import ints.WrappedIntArray;
import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * <p>Class {@code ParamEstimator} has static methods for estimating
 * local ancestry model parameters.</p>
 *
 * <p>Instances of class {@code ParamEstimator} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ParamEstimator {

    private ParamEstimator() {
        // private class to prevent instantiation
    }

    /**
     * Returns estimated local ancestry model parameters.
     * @param fixedParams the fixed parameters for a local ancestry inference
     * analysis
     * @param ibsHaps IBS segments
     * @param log object to which log output will be written
     * @return estimated local ancestry model parameters
     * @throws NullPointerException if
     * {@code (fixedParams == null) || (ibsHaps == null) || (log == null)}
     */
    public static ParamsInterface getParams(FixedParams fixedParams,
            IbsHaps ibsHaps, PrintWriter log) {
        ParamsInterface initParams = initParams(ibsHaps);
        ParamsInterface params = initParams;
        if (fixedParams.par().em()) {
            AdmixPar par = fixedParams.par();
            double minMu = par.min_mu();
            double deltaMu = par.delta_mu();
            Random rand = new Random(par.seed());
            boolean muHasConverged = false;
            double[] prevMu = new double[initParams.fixedParams().nAnc()];
            for (int j=0, n=par.em_its(); j<=n && muHasConverged==false; ++j) {
                // (n+1) iterations because j=0 calculates init parameter score
                WrappedIntArray selectedTargHaps = selectedTargHaps(ibsHaps.chromData(), rand);
                ParamEstimateData ape = estimateMuT(ibsHaps, params, selectedTargHaps);
                if (fixedParams.par().debug()) {
                    writeDiagnosticOutput(j, params);
                }
                params = ape.estimatedParams(initParams);
                double[] nextMu = params.mu();
                muHasConverged = muHasConverged(prevMu, nextMu, minMu, deltaMu);
                prevMu = nextMu;
                if (fixedParams.par().debug() && muHasConverged) {
                    writeDiagnosticOutput((j+1), params);
                }
            }
        }
        return params;
    }

    private static boolean muHasConverged(double[] prevMu, double[] nextMu,
            double minMu, double deltaMu) {
        for (int j=0; j<prevMu.length; ++j) {
            if (prevMu[j]>=minMu || nextMu[j]>=minMu) {
                if (Math.abs(prevMu[j] - nextMu[j]) > prevMu[j]*deltaMu) {
                    return false;
                }
            }
        }
        return true;
    }

    private static ParamsInterface initParams(IbsHaps ibsHaps){
        SelectedHaps selectedHaps = ibsHaps.selectedHaps();
        FixedParams fixedParams = selectedHaps.fixedParams();
        if (fixedParams.par().model()!=null) {
            return new ModelFileParams(fixedParams);
        }
        else {
            if (selectedHaps.includeRefHaps()) {
                ParamsInterface defaultParams = new DefaultParams(fixedParams);
                int nAnc = fixedParams.nAnc();
                ParamEstimateData[] estimatedParams = new ParamEstimateData[nAnc];
                for (int i=0; i<nAnc; ++i) {
                    AncSpecificParams initParams = new AncSpecificParams(fixedParams, i);
                    estimatedParams[i] = estimateRhoP(ibsHaps, initParams,
                            selectedHaps.panelToSelectedHapsIndices(i));
                }
                return new PartiallyUpdatedParams(estimatedParams, defaultParams);
            }
            else {
                return new DefaultParams(fixedParams);
            }
        }
    }

    private static void writeDiagnosticOutput(int it, ParamsInterface params) {
        File diagnosticFile = new File(params.fixedParams().par().out() + ".diag.model");
        if (it==0) {
            diagnosticFile.delete();
        }
        try (PrintWriter out = FileUtil.printWriter(diagnosticFile, true)) {
            out.println(Const.nl + "Iteration: " + it);
            out.println("=============");
            out.println();
            out.println(params.toString());
        }
    }

    private static ParamEstimateData estimateRhoP(IbsHaps ibsHaps,
            ParamsInterface oldParams, WrappedIntArray selectedHapIndices) {
        ParamEstimateData ape = new ParamEstimateData(oldParams);
        AdmixData data = new AdmixData(ibsHaps.chromData(), oldParams);
        AtomicInteger index = new AtomicInteger(0);
        int nThreads = oldParams.fixedParams().par().nthreads();
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    int nHaps = selectedHapIndices.size();
                    AdmixHmm hmm = new AdmixHmm(data, ibsHaps);
                    int i = index.getAndIncrement();
                    while (i<nHaps) {
                        hmm.runFwdBwdRhoP(selectedHapIndices.get(i));
                        i = index.getAndIncrement();
                    }
                    hmm.updateRhoP(ape);
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        MultiThreadUtils.shutdownExecService(es);
        return ape;
    }

    private static ParamEstimateData estimateMuT(IbsHaps ibsHaps,
            ParamsInterface oldParams, WrappedIntArray selectedTargHaps) {
        ParamEstimateData ape = new ParamEstimateData(oldParams);
        AdmixData data = new AdmixData(ibsHaps.chromData(), oldParams);
        AtomicInteger index = new AtomicInteger(0);
        int nThreads = oldParams.fixedParams().par().nthreads();
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    int nHaps = selectedTargHaps.size();
                    AdmixHmm hmm = new AdmixHmm(data, ibsHaps);
                    int i = index.getAndIncrement();
                    while (i<nHaps) {
                        hmm.runFwdBwdMuT(selectedTargHaps.get(i));
                        i = index.getAndIncrement();
                    }
                    hmm.updateMuT(ape);
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        MultiThreadUtils.shutdownExecService(es);
        return ape;
    }

    private static WrappedIntArray selectedTargHaps(AdmixChromData chromData,
            Random rand) {
        int nTargHaps = chromData.nTargHaps();
        int maxHaps = chromData.par().em_haps();
        int[] ia = IntStream.range(0, nTargHaps)
                .parallel()
                .toArray();
        if (nTargHaps<=maxHaps) {
            return new WrappedIntArray(ia);
        }
        else {
            Utilities.shuffle(ia, maxHaps, rand);
            int[] hapsToAnalyze = Arrays.copyOf(ia, maxHaps);
            Arrays.sort(hapsToAnalyze);
            return new WrappedIntArray(hapsToAnalyze);
        }
    }
}
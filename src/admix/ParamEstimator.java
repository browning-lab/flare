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
     * @param ibsHaps IBS segments
     * @param log object to which log output will be written
     * @return estimated local ancestry model parameters
     * @throws NullPointerException if
     * {@code ((ibsHaps == null) || (log == null))}
     */
    public static ParamsInterface getParams(IbsHaps ibsHaps, PrintWriter log) {
        SampleData sampleData = ibsHaps.chromData().sampleData();
        ParamsInterface params = initParams(ibsHaps);
        AdmixPar par = sampleData.par();
        if (par.em()) {
            Random rand = new Random(par.seed());
            double[] prevMu = params.studyMu();
            double[][] prevP = params.p();
            boolean finished = false;
            for (int j=0, n=par.em_its(); j<=n && finished==false; ++j) {
                // (n+1) iterations because j=0 calculates init parameter score
                WrappedIntArray selectedTargHaps = selectedTargHaps(ibsHaps.chromData(), rand);
                ParamEstimateData paramData = estimateParams(ibsHaps, params, selectedTargHaps);
                double score = paramData.meanMaxAncProb();
                if (sampleData.par().debug()) {
                    writeDiagnosticOutput(j, params, score);
                }
                params = new EstimatedParams(paramData);
                double[] nextMu = params.studyMu();
                double[][] nextP = params.p();
                finished = parametersHaveConverged(par, prevMu, nextMu, prevP, nextP);
                prevMu = nextMu;
                prevP = nextP;
                if (sampleData.par().debug() && finished) {
                    writeDiagnosticOutput((j+1), params, Double.NaN);
                }
            }
        }
        return params;
    }

    private static boolean parametersHaveConverged(AdmixPar par,
            double[] prevMu, double[] nextMu,
            double[][] prevP, double[][] nextP) {
        double deltaMu = par.delta_mu();
        for (int i=0; i<prevMu.length; ++i) {
            if (Math.abs(prevMu[i] - nextMu[i]) > deltaMu) {
                return false;
            }
        }
        if (par.update_p()) {
            double deltaP = par.delta_p();
            for (int i=0; i<prevP.length; ++i) {
                double[] prev = prevP[i];
                double[] next = nextP[i];
                for (int j=0; j<prev.length; ++j) {
                    if (Math.abs(prev[j] - next[j]) > deltaP) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    private static ParamsInterface initParams(IbsHaps ibsHaps){
        ObservedHaps observedHaps = ibsHaps.observedHaps();
        SampleData sampleData = observedHaps.sampleData();
        if (sampleData.par().model()!=null) {
            return new ModelFileParams(sampleData);
        }
        else {
            assert observedHaps.includeRefHaps()==true;
            int nAnc = sampleData.nAnc();
            ParamEstimateData[] paramData = new ParamEstimateData[nAnc];
            for (int i=0; i<nAnc; ++i) {
                AncSpecificParams ancParams = new AncSpecificParams(sampleData, i);
                paramData[i] = estimateRho(ibsHaps, ancParams,
                        observedHaps.panelToObservedHapsIndices(i));
            }
            return new PartiallyUpdatedParams(sampleData, paramData);
        }
    }

    private static ParamEstimateData estimateRho(IbsHaps ibsHaps,
            ParamsInterface oldParams, WrappedIntArray selectedHapIndices) {
        ParamEstimateData paramData = new ParamEstimateData(oldParams);
        AdmixData data = new AdmixData(ibsHaps.chromData(), oldParams);
        AtomicInteger index = new AtomicInteger(0);
        int nThreads = oldParams.sampleData().par().nthreads();
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    int nHaps = selectedHapIndices.size();
                    AdmixHmm hmm = new AdmixHmm(data, ibsHaps);
                    int i = index.getAndIncrement();
                    while (i<nHaps) {
                        hmm.runFwdBwdRho(selectedHapIndices.get(i));
                        i = index.getAndIncrement();
                    }
                    hmm.updateRho(paramData);
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        MultiThreadUtils.shutdownExecService(es);
        return paramData;
    }

    private static ParamEstimateData estimateParams(IbsHaps ibsHaps,
            ParamsInterface oldParams, WrappedIntArray selectedTargHaps) {
        ParamEstimateData paramData = new ParamEstimateData(oldParams);
        AdmixData data = new AdmixData(ibsHaps.chromData(), oldParams);
        AtomicInteger index = new AtomicInteger(0);
        int nThreads = oldParams.sampleData().par().nthreads();
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    int nHaps = selectedTargHaps.size();
                    AdmixHmm hmm = new AdmixHmm(data, ibsHaps);
                    int i = index.getAndIncrement();
                    while (i<nHaps) {
                        hmm.runFwdBwdEstimateParams(selectedTargHaps.get(i));
                        i = index.getAndIncrement();
                    }
                    hmm.updateParams(paramData);
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        MultiThreadUtils.shutdownExecService(es);
        return paramData;
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

    private static void writeDiagnosticOutput(int it, ParamsInterface params, 
            double meanMaxAncProb) {
        File diagnosticFile = new File(params.sampleData().par().out() + ".diag.model");
        if (it==0) {
            diagnosticFile.delete();
        }
        try (PrintWriter out = FileUtil.printWriter(diagnosticFile, true)) {
            if (Double.isFinite(meanMaxAncProb)) {
                out.println(Const.nl + "Before iteration: " + it);
                out.println("==========================");
                out.println("meanMaxAncProb: " + (float) meanMaxAncProb);
            }
            else {
                out.println(Const.nl + "Iteration: " + it);
                out.println("=============");
            }
            out.println();
            out.println(params.toString());
        }
    }
}
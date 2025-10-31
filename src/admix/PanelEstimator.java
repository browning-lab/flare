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
import blbutil.MultiThreadUtils;
import blbutil.Utilities;
import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * <p>Class {@code PanelEstimator} estimates reference panel probabilities
 * and mean recombination intensities for each target haplotype in a set
 * of marker windows.</p>
 *
 * <p>Instances of class {@code PanelEstimator} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class PanelEstimator {

    private PanelEstimator() {
        // private constructor to prevent instantiation
    }

    /**
     * Print the panel probabilities header line to the specified
     * {@code PrintWriter}.
     * @param sampleData reference and target sample metadata
     * @param out an object to which output will be written
     * @throws NullPointerException if
     * {@code ((sampleData == null) || (out == null))}
     */
    public static void printPanelProbsHeader(SampleData sampleData, PrintWriter out) {
        out.print("##nRefHaps=");
        out.println(sampleData.nRefHaps());
        out.print("#WIND");
        out.print(Const.tab);
        out.print("COORD");
        out.print(Const.tab);
        out.print("GT_HAP");
        out.print(Const.tab);
        out.print("SWITCH");
        for (int j=0, n=sampleData.nRefPanels(); j<n; ++j) {
            out.print(Const.tab);
            out.print(sampleData.refPanelId(j));
        }
        out.println();
    }

    /**
     * Estimates reference panel probabilities and mean recombination
     * intensities for each target haplotype in a set of marker windows and
     * writes these estimates to the specified {@code PrintWriter}.
     * @param ibsHaps the IBS haplotype segments for constructing composite
     * reference haplotypes
     * @param out an object to which output will be written
     * @throws NullPointerException if
     * {@code ((ibsHaps == null) || (out == null))}
     */
    public static void writePanelProbs(IbsHaps ibsHaps, PrintWriter out) {
        AdmixChromData chromData = ibsHaps.chromData();
        EstimatedPanels estimatedPanels = new EstimatedPanels(chromData);
        int nSelectedTargHaps = estimatedPanels.nSelectedTargHaps();
        AtomicInteger index = new AtomicInteger(0);
        int nThreads = chromData.par().nthreads();
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    PanelHmm hmm = new PanelHmm(ibsHaps);
                    int i = index.getAndIncrement();
                    while (i<nSelectedTargHaps) {
                        hmm.runFwdBwd(i, estimatedPanels);
                        i = index.getAndIncrement();
                    }
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        MultiThreadUtils.shutdownExecService(es);

        estimatedPanels.writePanelProbabilities(out);
    }
}
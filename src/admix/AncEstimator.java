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

import blbutil.MultiThreadUtils;
import blbutil.Utilities;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * <p>Class {@code AncEstimator} has static methods for estimating
 * local ancestry model parameters and local ancestry.</p>
 *
 * <p>Instances of class {@code AncEstimator} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AncEstimator {

    private AncEstimator() {
        // private constructor to prevent instantiation
    }

   /**
     * Runs a local ancestry inference analysis and write the output to the
     * specified {@code AdmixWriter}.
     * @param ibsHaps IBS segments
     * @param params the local ancestry model parameters
     * @param admixWriter an object to which output will be written
     * @throws IllegalArgumentException if
     * {@code params.par() != admixWriter.par()}
     * @throws NullPointerException if
     * {@code (ibsHaps == null) || (params == null) || (admixWriter == null)}
     */
    public static void estAncestry(IbsHaps ibsHaps, ParamsInterface params,
            AdmixWriter admixWriter) {
        if (params.fixedParams().par() != admixWriter.par()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        EstimatedAncestry estAnc = estAncestry(ibsHaps, params);
        admixWriter.writeAncestry(estAnc);
    }

    private static EstimatedAncestry estAncestry(IbsHaps ibsHaps, ParamsInterface params) {
        AdmixData data = new AdmixData(ibsHaps.chromData(), params);
        int nTargHaps = ibsHaps.chromData().nTargHaps();
        EstimatedAncestry estAnc = new EstimatedAncestry(data);
        AtomicInteger index = new AtomicInteger(0);
        int nThreads = params.fixedParams().par().nthreads();
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    AdmixHmm hmm = new AdmixHmm(data, ibsHaps);
                    int i = index.getAndIncrement();
                    while (i<nTargHaps) {
                        hmm.runFwdBwdAnc(i, estAnc);
                        i = index.getAndIncrement();
                    }
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        MultiThreadUtils.shutdownExecService(es);
        return estAnc;
    }
}
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

import beagleutil.PbwtDivUpdater;
import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import vcf.RefGT;
import vcf.Steps;

/**
 * <p>Class {@code AdmixCodedSteps} divides phased genotype data
 * into non-overlapping chromosome intervals (steps), indexes the unique
 * allele sequences in each interval, and stores a map of haplotype
 * index to allele sequence index for each interval.</p>
 *
 * <p>Instances of class {@code CodedSteps} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixCodedSteps {

    private final RefGT targRefGT;
    private final Steps steps;
    private final IndexArray[] codedSteps;

    /**
     * Constructs a new {@code AdmixCodedSteps} instance from the specified
     * data.
     * @param chromData the input arguments and genotype data for a chromosome
     * @throws NullPointerException if {@code chromData == null}
     */
    public AdmixCodedSteps(AdmixChromData chromData) {
        this.targRefGT = chromData.targRefGT();
        this.steps = new Steps(chromData.map(), chromData.par().ibs_step());
        this.codedSteps = codedSteps(targRefGT, steps, chromData.par().nthreads());
    }

    private static IndexArray[] codedSteps(RefGT targRefGT, Steps steps, int nThreads) {
        int maxStepsPerBatch = 512;
        int nStepsPerBatch0 = (steps.size() + nThreads - 1)/nThreads;
        while (nStepsPerBatch0>maxStepsPerBatch) {
            nStepsPerBatch0 >>= 1;
        }
        int stepsPerBatch = nStepsPerBatch0;
        int nBatches = (steps.size() + (stepsPerBatch-1)) / stepsPerBatch;
        return IntStream.range(0, nBatches)
                .parallel()
                .boxed()
                .flatMap(batch -> codedSteps(targRefGT, steps, batch, stepsPerBatch))
                .toArray(IndexArray[]::new);
    }

    private static Stream<IndexArray> codedSteps(RefGT targRefGT, Steps steps,
            int batch, int batchSize) {
        int nHaps = targRefGT.nHaps();
        int startStep = batch*batchSize;
        int endStep = Math.min(startStep + batchSize, steps.size());
        int nSteps = endStep - startStep;
        IndexArray[] hapToSeq = new IndexArray[nSteps];
        int[] p = new int[nHaps];
        int[] d = new int[nHaps];
        int[] hap2Seq = new int[nHaps];
        PbwtDivUpdater pbwtUpdater = new PbwtDivUpdater(nHaps);

        int mStart = steps.start(startStep);
        for (int step=startStep; step<endStep; ++step) {
            int mEnd = steps.end(step);
            if ((mEnd-mStart)==1) {
                int nAlleles = targRefGT.marker(mStart).nAlleles();
                IntArray ia = rec(targRefGT, mStart);
                hapToSeq[step-startStep] = new IndexArray(ia, nAlleles);
            }
            else {
                setToIdentity(p);
                Arrays.fill(d, mStart);
                for (int m=mStart; m<mEnd; ++m) {
                    int nAlleles = targRefGT.marker(m).nAlleles();
                    pbwtUpdater.fwdUpdate(rec(targRefGT, m), nAlleles, m, p, d);
                }
                int index = -1;
                for (int j=0; j<nHaps; ++j) {
                    hap2Seq[p[j]] = (d[j]>mStart) ? ++index : index;
                }
                hapToSeq[step-startStep] = new IndexArray(hap2Seq, ++index);
            }
            mStart = mEnd;
        }
        return Arrays.stream(hapToSeq);
    }

    private static IntArray rec(final RefGT gt, final int marker) {
        return new IntArray() {

            @Override
            public int size() {
                return gt.nHaps();
            }

            @Override
            public int get(int index) {
                return gt.allele(marker, index);
            }
        };
    }

    private static void setToIdentity(int[] p) {
        for (int j=0; j<p.length; ++j) {
            p[j] = j;
        }
    }

    /**
     * Returns a map from haplotype index to allele sequence index
     * for the specified step
     * @param step a step index
     * @return a map from haplotype index to allele sequence index
     * for the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.steps().size()}
     */
    public IndexArray get(int step) {
        return codedSteps[step];
    }


    /**
     * Returns the phased reference and target genotypes with
     * target samples preceding reference samples.
     * @return the phased reference and target genotypes with
     * target samples preceding reference samples
     */
    public RefGT targRefGT() {
        return targRefGT;
    }

    /**
     * Returns the partition of the chromosome markers into non-overlapping
     * intervals.
     * @return the partition of the chromosome markers into non-overlapping
     * intervals
     */
    public Steps steps() {
        return steps;
    }
}

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

import beagleutil.PbwtDivUpdater;
import ints.IndexArray;
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.stream.IntStream;
import vcf.Steps;

/**
 * <p>Class {@code IbsHaps} partitions a chromosome into genomic intervals
 * and for each genomic interval and target haplotype stores a reference
 * haplotype that is identical by state with the target haplotype in a long
 * segment that contains the genomic interval if such a reference haplotype
 * exists.</p>
 *
 * <p>Instances of class {@code IbsHaps} are immutable.</p>
 *
 * <p>Reference: Durbin, R. 2014. Bioinformatics 30(9):1266–1272.
 * doi:10.1093/bioinformatics/btu014</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IbsHaps {

    private final AdmixChromData chromData;
    private final Steps steps;
    private final int nIbsHaps;
    private final SelectedHaps selectedHaps;
    private final WrappedIntArray[] ibsHaps; // [step][selected hap]

    /**
     * Constructs a new {@code AdmixIbsHaps} instance from the specified data.
     * @param chromData immutable input data for local ancestry inference
     * on a chromosome
     * @param selectedHaps the haplotypes for which IBS segments will be stored.
     * @throws NullPointerException if
     * {@code (chromData == null) || (selectedHaps == null)}
     */
    public IbsHaps(AdmixChromData chromData, SelectedHaps selectedHaps) {
        AdmixCodedSteps codedSteps = new AdmixCodedSteps(chromData);
        AdmixPar par = chromData.par();
        int nThreads = par.nthreads();
        int nSteps = codedSteps.steps().size();
        int batchSize = (nSteps + nThreads - 1) / nThreads;
        int nBatches = (nSteps + batchSize - 1) / batchSize;
        this.chromData = chromData;
        this.steps = codedSteps.steps();
        this.nIbsHaps = par.ibs_haps();
        this.selectedHaps = selectedHaps;
        final WrappedIntArray hapToSelectedHapsIndex
                = hapToSelectedHapsIndex(selectedHaps.selectedHaps(), chromData.nHaps());
        this.ibsHaps = IntStream.range(0, nBatches)
                .parallel()
                .mapToObj(batch -> ibsHaps(chromData, selectedHaps, codedSteps,
                        hapToSelectedHapsIndex, batch, batchSize))
                .flatMap(a -> Arrays.stream(a))
                .toArray(WrappedIntArray[]::new);

    }

    private static WrappedIntArray hapToSelectedHapsIndex(
            WrappedIntArray selectedHaps, int nHaps) {
        int[] hapToSelectedHapsIndex = IntStream.range(0, nHaps)
                .parallel()
                .map(j -> -1)
                .toArray();
        for (int j=0, n=selectedHaps.size(); j<n; ++j) {
            hapToSelectedHapsIndex[selectedHaps.get(j)] = j;
        }
        return new WrappedIntArray(hapToSelectedHapsIndex);
    }

    private static WrappedIntArray[] ibsHaps(AdmixChromData chromData,
            SelectedHaps selectedHaps, AdmixCodedSteps codedSteps,
            WrappedIntArray hapToSelectedHapsIndex, int batch, int batchSize) {
        AdmixPar par = chromData.par();
        int nSteps = codedSteps.steps().size();
        int start = batch*batchSize;
        int end = Math.min(start + batchSize, nSteps);
        int nOverlapSteps = (int) Math.rint(par.ibs_buffer() / par.ibs_step());
        int overlapStart = Math.max(0, start - nOverlapSteps);
        int overlapEnd = Math.min(end + nOverlapSteps, nSteps);

        PbwtDivUpdater pbwt = new PbwtDivUpdater(chromData.nHaps());
        int nSelectedHaps = selectedHaps.selectedHaps().size();
        int[][] ibsHaps0 = new int[end-start][nSelectedHaps*par.ibs_haps()];
        fwdIbsHaps(chromData, codedSteps, hapToSelectedHapsIndex, pbwt, ibsHaps0, overlapStart, start);
        bwdIbsHaps(chromData, codedSteps, hapToSelectedHapsIndex, pbwt, ibsHaps0, end, overlapEnd);
        return Arrays.stream(ibsHaps0)
                .map(ia -> new WrappedIntArray(ia))
                .toArray(WrappedIntArray[]::new);
    }

    private static void fwdIbsHaps(AdmixChromData chromData,
            AdmixCodedSteps codedSteps, WrappedIntArray hapToSelectedHapsIndex,
            PbwtDivUpdater pbwt, int[][] ibsHaps, int overlapStart, int start) {
        int nHaps = pbwt.nHaps();
        int[] a = IntStream.range(0, nHaps).toArray();
        int[] d = IntStream.range(0, nHaps+1).map(j -> overlapStart).toArray(); // last entry is sentinal
        for (int step=overlapStart; step<start; ++step) {
            IndexArray ia = codedSteps.get(step);
            pbwt.fwdUpdate(ia.intArray(), ia.valueSize(), step, a, d);
        }
        for (int j=0, step=start; j<ibsHaps.length; ++j, ++step) {
            IndexArray ia = codedSteps.get(step);
            pbwt.fwdUpdate(ia.intArray(), ia.valueSize(), step, a, d);
            setfwdIbsHaps(chromData, step, hapToSelectedHapsIndex, a, d, ibsHaps[j]);
        }
    }

    private static void setfwdIbsHaps(AdmixChromData chromData, int step,
            WrappedIntArray hapToSelectedHapsIndex,
            int[] a, int[] d, int[] ibsHaps) {
        assert d[0] == (step + 1);
        int nIbsHaps = chromData.par().ibs_haps();
        int nFwdIbsHaps = nIbsHaps >> 1;
        int nTargHaps = chromData.nTargHaps();
        d[a.length] = step + 1;  // set sentinal
        for (int i=0; i<a.length; ++i) {
            int selectedHapsIndex = hapToSelectedHapsIndex.get(a[i]);
            if (selectedHapsIndex>=0) {
                int index = 0;
                int start = selectedHapsIndex*nIbsHaps;
                int u = i;          // inclusive start
                int v = i + 1;      // exclusive end
                int uNextMatchStart = d[u];
                int vNextMatchStart = d[v];
                while (index<nFwdIbsHaps
                        && (uNextMatchStart<=step || vNextMatchStart<=step)) {
                    if (vNextMatchStart<=uNextMatchStart) {
                        if (v<a.length && a[v]>=nTargHaps) {
                            ibsHaps[start + (index++)] = a[v];
                        }
                        vNextMatchStart = Math.max(d[++v], vNextMatchStart);
                    }
                    else {
                        uNextMatchStart = Math.max(d[--u], uNextMatchStart);
                        if (a[u]>=nTargHaps) {
                            ibsHaps[start + (index++)] = a[u];
                        }
                    }
                }
                while (index<nFwdIbsHaps) {
                    ibsHaps[start + (index++)] = -1;
                }
            }
        }
    }

    private static void bwdIbsHaps(AdmixChromData chromData, AdmixCodedSteps codedSteps,
            WrappedIntArray hapToSelectedHapsIndex, PbwtDivUpdater pbwt,
            int[][] ibsHaps, int end, int overlapEnd) {
        int nHaps = pbwt.nHaps();
        int[] a = IntStream.range(0, nHaps).toArray();
        int[] d = IntStream.range(0, nHaps+1).map(j -> (overlapEnd-1)).toArray(); // last entry is sentinal
        for (int step=(overlapEnd-1); step>=end; --step) {
            IndexArray ia = codedSteps.get(step);
            pbwt.bwdUpdate(ia.intArray(), ia.valueSize(), step, a, d);
        }
        for (int j=(ibsHaps.length-1), step=(end-1); j>=0; --j, --step) {
            IndexArray ia = codedSteps.get(step);
            pbwt.bwdUpdate(ia.intArray(), ia.valueSize(), step, a, d);
            setBwdIbsHaps(chromData, step, hapToSelectedHapsIndex, a, d, ibsHaps[j]);
        }
    }

    private static void setBwdIbsHaps(AdmixChromData chromData,
            int step, WrappedIntArray hapToSelectedHapsIndex,
            int[] a, int[] d, int[] ibsHaps) {
        d[0] = d[a.length] = step - 1;  // set sentinals
        int nIbsHaps = chromData.par().ibs_haps();
        int nFwdIbsHaps = nIbsHaps >> 1;
        int nBwdIbsHaps = nIbsHaps - nFwdIbsHaps;
        int nTargHaps = chromData.nTargHaps();
        // no need to save and restore old d[0], d[a.length] values
        for (int i=0; i<a.length; ++i) {
            int selectedHapsIndex = hapToSelectedHapsIndex.get(a[i]);
            if (selectedHapsIndex>=0) {
                int index = 0;
                int start = selectedHapsIndex*nIbsHaps + nFwdIbsHaps;
                int u = i;          // inclusive start
                int v = i + 1;      // exclusive end
                int uNextMatchEnd = d[u];
                int vNextMatchEnd = d[v];
                while (index<nBwdIbsHaps
                        && (step<=uNextMatchEnd || step<=vNextMatchEnd)) {
                    if (uNextMatchEnd<=vNextMatchEnd) {
                        if (v<a.length && a[v]>=nTargHaps) {
                            ibsHaps[start + (index++)] = a[v];
                        }
                        vNextMatchEnd = Math.min(d[++v], vNextMatchEnd);
                    }
                    else {
                        uNextMatchEnd = Math.min(d[--u], uNextMatchEnd);
                        if (a[u]>=nTargHaps) {
                            ibsHaps[start + (index++)] = a[u];
                        }
                    }
                }
                while (index<nBwdIbsHaps) {
                    ibsHaps[start + (index++)] = -1;
                }
            }
        }
    }

    /**
     * Returns the immutable input data for local ancestry inference on a
     * chromosome.
     * @return the immutable input data for local ancestry inference on a
     * chromosome
     */
    public AdmixChromData chromData() {
        return chromData;
    }

    /**
     * Returns the haplotype indices with stored IBS segments.
     * The first reference haplotype has index
     * {@code this.chromData().nTargHaps()}.
     * @return the haplotype indices with stored IBS segments
     */
    public SelectedHaps selectedHaps() {
        return selectedHaps;
    }

    /**
     * Returns the intervals with IBS haplotype data.
     * @return the intervals with IBS haplotype data
     */
    public Steps steps() {
        return steps;
    }

    /**
     * Adds reference haplotypes that are identical by state with the
     * specified haplotype in the specified genomic interval to
     * the specified {@code AdmixStates} object.  Up to
     * {@code this.chromData().par().ibs_haps()} reference haplotypes will be
     * added via the {@code states.addIbsHap()} method.
     * @param hapListIndex an index in {@code this.selectedHaps().hapList()}
     * @param step an index of a genomic interval
     * @param states an object that constructs a HMM state space
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.selectedHaps().hapList().size()}
     * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.steps().size()}
     * @throws NullPointerException if {@code states == null}
     */
    public void addIbsHaps(int hapListIndex, int step, AdmixStates states) {
        int start = hapListIndex*nIbsHaps;
        int end = start + nIbsHaps;
        for (int j=start; j<end; ++j) {
            int ibsHap = ibsHaps[step].get(j);
            if (ibsHap>=0) {
                states.addIbsHap(ibsHap, step);
            }
        }
    }
}

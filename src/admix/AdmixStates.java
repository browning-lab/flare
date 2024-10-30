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

import ints.IntIntMap;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.stream.IntStream;
import beagleutil.CompHapSegment;
import blbutil.Utilities;
import ints.IntArray;
import ints.IntList;
import ints.WrappedIntArray;
import java.util.Arrays;
import vcf.Steps;
import vcf.RefGT;

/**
 * <p>Class {@code AdmixStates} has methods for constructing the state space
 * of a Li and the Stephens HMM for a haplotype.</p>
 *
 * <p>Instances of {@code AdmixStates} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AdmixStates {

    private static final int NIL = -103;
    private final AdmixChromData chromData;
    private final IbsHaps ibsHaps;
    private final WrappedIntArray hapList;
    private final IntArray refHapToPanel;
    private final Steps steps;
    private final int nSteps;
    private final RefGT targRefGT;
    private final int nTargHaps;
    private final int nMarkers;
    private final int maxStates;
    private final int minSteps;

    private final IntIntMap hapToLastIbsStep;
    private final PriorityQueue<CompHapSegment> q;
    private final IntList[] compHapHap;
    private final IntList[] compHapEnd;

    private final int[] compHapToListIndex;
    private final int[] compHapToHap;
    private final int[] compHapToPanel;
    private final int[] compHapToEnd;

    /**
     * Constructs a new {@code AdmixStates} instance from the specified data.
     * @param data immutable input data for a chromosome and the
     * analysis parameters for local ancestry inference
     * @param ibsHaps the IBS haplotype segments
     * @throws IllegalArgumentException if
     * {@code data.chromData() != ibsHaps.chromData()}
     * @throws NullPointerException if
     * {@code (chromData == null) || (ibsHaps == null)}
     */
    public AdmixStates(AdmixData data, IbsHaps ibsHaps) {
        if (data.chromData()!=ibsHaps.chromData()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        AdmixPar par = data.params().fixedParams().par();
        this.chromData = data.chromData();
        this.refHapToPanel = data.params().fixedParams().refHapToPanel();
        this.ibsHaps = ibsHaps;
        this.hapList = ibsHaps.selectedHaps().selectedHaps();
        this.steps = ibsHaps.steps();
        this.nSteps = steps.size();
        this.targRefGT = chromData.targRefGT();
        this.nTargHaps = chromData.nTargHaps();
        this.nMarkers = targRefGT.nMarkers();
        this.maxStates = Math.min(chromData.nRefHaps(), par.states());
        this.minSteps = (int) Math.ceil(par.ibs_recycle()/par.ibs_step());

        this.hapToLastIbsStep = new IntIntMap(maxStates);
        this.q = new PriorityQueue<>(maxStates);
        this.compHapHap = intLists(maxStates);
        this.compHapEnd = intLists(maxStates);

        this.compHapToListIndex = new int[maxStates];
        this.compHapToHap = new int[maxStates];
        this.compHapToPanel = new int[maxStates];
        this.compHapToEnd = new int[maxStates];
    }

    private static IntList[] intLists(int maxStates) {
        return IntStream.range(0, maxStates)
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
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
     * Returns the IBS haplotype segments.
     * @return the IBS haplotype segments
     */
    public IbsHaps ibsHaps() {
        return ibsHaps;
    }

    /**
     * Returns the maximum number of model states.
     * @return the maximum number of model states
     */
    public int maxStates() {
        return maxStates;
    }

    /**
     * Stores the Li and Stephens HMM for the specified target sample in
     * the specified arrays.  The contract for this method is undefined if
     * any row of either specified two-dimensional array has fewer than
     * {@code this.maxStates()} elements.
     *
     * @param hapListIndex an index in {@code this.selectedHaps().hapList()}
     * @param refPanel the two-dimensional array in which the reference panel
     * for each model state will be stored
     * @param nMismatches the two-dimensional array in which the number of
     * allele mismatches (0 or 1) for each model state will be stored
     * @return the number of state alleles at each marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code (hapListIndex < 0) ||
     * (hapListIndex >= this.ibsHaps().selectedHaps().hapList().size())}
     * @throws IndexOutOfBoundsException if
     * {@code refPanel.length < this.chromData().targRefGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code nMismatches.length < this.chromData().targRefGT().nMarkers()}
     * @throws NullPointerException if any array is {@code null}
     */
    public int ibsStates(int hapListIndex, short[][] refPanel, byte[][] nMismatches) {
        setCompRefHaps(hapListIndex);
        int hap = hapList.get(hapListIndex);
        int nStates = copyData(hap, refPanel, nMismatches);
        return nStates;
    }

    private void setCompRefHaps(int hapListIndex) {
        initializeFields();
        for (int step=0; step<nSteps; ++step) {
            ibsHaps.addIbsHaps(hapListIndex, step, this);
        }
        if (q.isEmpty()) {
            fillQWithRandomHaps(hapList.get(hapListIndex));
        }
    }

    private void initializeFields() {
        hapToLastIbsStep.clear();
        for (int j=0, n=q.size(); j<n; ++j) {
            compHapHap[j].clear();
            compHapEnd[j].clear();
        }
        q.clear();
    }

    void addIbsHap(int hap, int ibsStep) {
        if (hapToLastIbsStep.get(hap, NIL)==NIL) { // hap not currently in q
            updateHeadOfQ();
            if (q.size()==maxStates
                    || (q.isEmpty()==false && ibsStep - q.peek().lastIbsStep() >= minSteps)) {
                CompHapSegment head = q.poll();
                int startMarker = steps.start((head.lastIbsStep() + ibsStep) >>> 1);
                hapToLastIbsStep.remove(head.hap());
                compHapHap[head.compHapIndex()].add(hap);          // hap of new segment
                compHapEnd[head.compHapIndex()].add(startMarker);  // end of old segment
                head.updateSegment(hap, startMarker, ibsStep);
                q.offer(head);
            }
            else {
                int compHapIndex = q.size();
                int startMarker = 0;
                compHapHap[compHapIndex].add(hap); // hap of new segment
                q.offer(new CompHapSegment(hap, startMarker, ibsStep, compHapIndex));
            }
        }
        hapToLastIbsStep.put(hap, ibsStep);
    }

    private void updateHeadOfQ() {
        CompHapSegment head = q.peek();
        if (head!=null) {
            int lastIbsStep = hapToLastIbsStep.get(head.hap(), NIL);
            while (head.lastIbsStep()!=lastIbsStep) {
                head = q.poll();
                head.setLastIbsStep(lastIbsStep);
                q.offer(head);
                head = q.peek();
                lastIbsStep = hapToLastIbsStep.get(head.hap(), NIL);
            }
        }
    }

    private void flushQ() {
        updateHeadOfQ();
        CompHapSegment head = q.poll();
        while (head!=null) {
            hapToLastIbsStep.remove(head.hap());
            compHapHap[head.compHapIndex()].add(head.hap());       // hap of new segment
            if (head.lastIbsStep()>0) {
                int startMarker = steps.start((head.lastIbsStep() + nSteps - 1) >>> 1);
                compHapEnd[head.compHapIndex()].add(startMarker);  // end of old segment
            }
            updateHeadOfQ();
            head = q.poll();
        }
    }

    private int copyData(int hap, short[][] refPanel, byte[][] nMismatches) {
        int nCompHaps = q.size();
        // flushQ();    // xxx uncommenting decreases accuracy and changes the regression test
        initializeCopy(nCompHaps);
        for (int m=0; m<nMarkers; ++m) {
            short[] panel = refPanel[m];
            byte[] mismatches = nMismatches[m];
            if (panel.length < maxStates) {
                throw new IllegalArgumentException(String.valueOf(panel.length));
            }
            if (mismatches.length < maxStates) {
                throw new IllegalArgumentException(String.valueOf(mismatches.length));
            }
            int targAllele = targRefGT.allele(m, hap);
            for (int j=0; j<nCompHaps; ++j) {
                if (m==compHapToEnd[j]) {
                    int index = ++compHapToListIndex[j];
                    compHapToHap[j] = compHapHap[j].get(index);
                    compHapToPanel[j] = refHapToPanel.get(compHapToHap[j] - nTargHaps);
                    compHapToEnd[j] = compHapEnd[j].get(index);
                }
                refPanel[m][j] = (short) compHapToPanel[j];
                nMismatches[m][j] = targRefGT.allele(m, compHapToHap[j])==targAllele ? (byte) 0 : (byte) 1;
            }
        }
        return nCompHaps;
    }

    private void initializeCopy(int nCompHaps) {
        for (int j=0; j<nCompHaps; ++j) {
            compHapEnd[j].add(nMarkers); // add end of last segment
            compHapToListIndex[j] = 0;
            compHapToHap[j] = compHapHap[j].get(0);
            compHapToPanel[j] = refHapToPanel.get(compHapToHap[j]-nTargHaps);
            compHapToEnd[j] = compHapEnd[j].get(0);
        }
    }

    private void fillQWithRandomHaps(int hap) {
        assert q.isEmpty();
        int ibsStep = 0;
        int startMarker = 0;
        int compHapIndex = 0;
        int[] randomRefHaps = randomRefHaps(hap);
        for (int h : randomRefHaps) {
            assert (h!=hap && hapToLastIbsStep.get(h, NIL)==NIL);
            compHapHap[compHapIndex].add(h);
            q.add(new CompHapSegment(h, startMarker, ibsStep, compHapIndex++));
            hapToLastIbsStep.put(h, ibsStep);
        }
    }

    private int[] randomRefHaps(int hap) {
        Random rand = new Random(chromData.par().seed() + hap);
        int[] haps = IntStream.range(nTargHaps, targRefGT.nHaps())
                .filter(h -> (h!=hap))
                .toArray();
        Utilities.shuffle(haps, maxStates, rand);
        if (haps.length > maxStates) {
            haps = Arrays.copyOf(haps, maxStates);
        }
        return haps;
    }
}

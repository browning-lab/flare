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

import blbutil.FloatArray;
import ints.WrappedIntArray;
import java.util.Arrays;

/**
 * <p>Class {@code PanelHmm} implements the forward and backward
 * algorithms for a modified Li and Stephens hidden Markov model
 * that is used to estimate reference panel probabilities and
 * mean recombination intensities for a target haplotype in a set
 * of marker windows.</p>
 *
 * <p>Instances of class {@code PanelHmm} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PanelHmm {

    private final AdmixChromData chromData;
    private final int nMarkers;
    private final int nPanels;

    private final AdmixStates states;
    private final short[][] refPanel;
    private final byte[][] nMismatches;
    private final double[] pMismatch;
    private final FloatArray pRecomb;

    private int nStates = -1;
    private double invNStates = Double.NaN;
    private double hFactor = Double.NaN;
    private final double[] fwd;
    private final double[] bwd;
    private final double[] panelProbs;
    private final double[] recombIntensities;

    private final int windowSize;
    private final int nWindows;
    private final double[][] bwdCheckPts; // [checkpoint][state]
    private final double[][] bwdWindow;   // [offset][state]

    /**
     * Constructs a new {@code PanelHmm} instance for the specified data.
     *
     * @param ibsHaps the IBS haplotype segments for constructing composite
     * reference haplotypes
     * @throws NullPointerException if {@code (ibsHaps == null)}
     */
    public PanelHmm(IbsHaps ibsHaps) {
        this.chromData = ibsHaps.chromData();
        this.nMarkers = chromData.targRefGT().nMarkers();
        this.nPanels = chromData.sampleData().nRefPanels();
        double theta = ParamsInterface.liStephensPMismatch(chromData.nRefHaps());
        this.pMismatch = new double[] {1 - theta, theta};
        int ne = chromData.par().panel_ne();
        float recombIntensity = 0.04f*ne/chromData.nRefHaps(); // 0.01 converts from cM to Morgans
        this.pRecomb = chromData.map().pRecomb(recombIntensity);

        this.states = new AdmixStates(ibsHaps);
        this.refPanel = new short[nMarkers][states.maxStates()];
        this.nMismatches = new byte[nMarkers][states.maxStates()];

        this.fwd = new double[states.maxStates()];
        this.bwd = new double[states.maxStates()];
        this.panelProbs = new double[nMarkers * nPanels];
        this.recombIntensities = new double[nMarkers];
        this.windowSize = (int) Math.ceil(Math.sqrt(nMarkers));
        this.nWindows = (nMarkers + windowSize - 1)/windowSize;
        this.bwdCheckPts = new double[nWindows][states.maxStates()];
        this.bwdWindow = new double[windowSize][states.maxStates()];
    }

    /**
     * Returns the input genotype data for local ancestry inference on a
     * chromosome.
     * @return the genotype input data for local ancestry inference on a
     * chromosome
     */
    public AdmixChromData chromData() {
        return chromData;
    }

    /**
     * Returns the list of haplotypes with stored IBS segments
     * @return the list of haplotypes with stored IBS segments
     */
    public WrappedIntArray hapList() {
        return states.ibsHaps().observedHaps().hapList();
    }

    /**
     * Estimates reference panel probabilities and recombination intensities
     * for the specified target haplotype and stores these estimates
     * in the specified {@code EstimatedPanels} object.
     * @param index index of a selected target haplotype
     * @param estPanels the object in which estimated panel probabilities
     * will be stored
     * @throws IndexOutOfBoundsException if
     * {@code ((selectedIndex < 0) || (selectedIndex >= estPanels.nSelectedTargHaps())}
     * @throws NullPointerException if {@code estPanel == null}
     */
    public void runFwdBwd(int index, EstimatedPanels estPanels) {
        if (index<0 || index>=estPanels.nSelectedTargHaps()) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int targHap = estPanels.selectedTargHap(index);
        nStates = states.ibsStates(targHap, refPanel, nMismatches);
        invNStates = 1.0/nStates;
        this.hFactor = nStates / (nStates - 1.0);
        Arrays.fill(panelProbs, 0.0);
        Arrays.fill(recombIntensities, 0.0);
        Arrays.fill(fwd, 0, nStates, invNStates);
        setBwdCheckPoints();
        int nextWindowStart = 0;
        int window = -1;
        int offset = Integer.MIN_VALUE;
        double fwdSum = 1.0f;
        for (int m=0; m<nMarkers; ++m)  {
            if (m==nextWindowStart) {
                offset = -1;
                fillBwdWindow(++window);
                nextWindowStart += windowSize;
            }
            fwdSum = fwdUpdate(m, bwdWindow[++offset], fwdSum, hFactor);
        }
        estPanels.set(index, panelProbs, recombIntensities);
    }

    private void setBwdCheckPoints() {
        int windowIndex = nWindows;
        Arrays.fill(bwd, 0, nStates, invNStates);
        System.arraycopy(bwd, 0, bwdCheckPts[--windowIndex], 0, nStates);
        int nextCheckPt = ((nWindows-1)*windowSize) - 1; // checkpoint is last marker in window
        for (int m=nMarkers-2; m>=0; --m) {
            bwdUpdate(m);
            if (m==nextCheckPt) {
                System.arraycopy(bwd, 0, bwdCheckPts[--windowIndex], 0, nStates);
                nextCheckPt -= windowSize;
            }
        }
        assert windowIndex==0;
    }

    private void bwdUpdate(int m) {
        int mP1 = m+1;
        assert mP1 < nMarkers;
        double bwdSum = 0.0;
        byte[] mismatchCnt = nMismatches[mP1];
        for (int h=0; h<nStates; ++h) {
            bwd[h] *= pMismatch[mismatchCnt[h]];
            bwdSum += bwd[h];
        }
        float pRec = pRecomb.get(mP1);
        double scale = (1.0 - pRec)/bwdSum;
        double shift = pRec/nStates;
        for (int h=0; h<nStates; ++h) {
            bwd[h] = scale*bwd[h] + shift;
        }
    }

    private double fwdUpdate(int m, double[] bwdM, double fwdSum, double hFactor) {
        int offsetStart = m*nPanels;
        float pRec = pRecomb.get(m);
        double shift = pRec/nStates;
        double scale = (1.0 - pRec)/fwdSum;
        double noSwitchScale = ((1.0f - pRec) + shift)/fwdSum;
        float jointStateSum = 0.0f;
        double stateSum = 0.0;
        byte[] mismatchCnt = nMismatches[m];
        short[] hap2Panel = refPanel[m];
        double nextFwdSum = 0.0;
        for (int h=0; h<nStates; ++h) {
            double em = pMismatch[mismatchCnt[h]];
            jointStateSum += bwdM[h]*em*noSwitchScale*fwd[h];
            fwd[h] = m==0 ? em : em*(scale*fwd[h] + shift);
            nextFwdSum += fwd[h];
            stateSum += fwd[h]*bwdM[h];
            panelProbs[offsetStart + hap2Panel[h]] += fwd[h]*bwdM[h];
        }
        recombIntensities[m] = hFactor*(1.0f - jointStateSum/stateSum);
        scaleToSumTo1(panelProbs, offsetStart, offsetStart + nPanels);
        return nextFwdSum;
    }

    private void fillBwdWindow(int window) {
        // checkpoint is last marker in window
        int m = Math.min((window+1)*windowSize, nMarkers) - 1;
        int lastOffset = m - window*windowSize;
        System.arraycopy(bwdCheckPts[window], 0, bwd, 0, nStates);
        System.arraycopy(bwdCheckPts[window], 0, bwdWindow[lastOffset], 0, nStates);
        for (int j=(lastOffset-1); j>=0; --j) {
            bwdUpdate(--m);
            System.arraycopy(bwd, 0, bwdWindow[j], 0, nStates);
        }
    }

    private static void scaleToSumTo1(double[] da, int start, int end) {
        double sum = 0.0;
        for (int j=start; j<end; ++j) {
            sum += da[j];
        }
        double factor = 1.0/sum;
        for (int j=start; j<end; ++j) {
            da[j] *= factor;
        }
    }
}

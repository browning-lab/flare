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

import ints.WrappedIntArray;
import java.util.Arrays;

/**
 * <p>Class {@code AdmixHmm} implements the forward and backward
 * algorithms for a modified Li and Stephens hidden Markov model.</p>
 *
 * <p>Instances of class {@code AdmixHmm} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixHmm {

    private final AdmixChromData chromData;
    private final ParamsInterface params;
    private final int nMarkers;
    private final int nTargHaps;
    private final int nAnc;
    private final int nPanels;
    private final float minAncProb;
    private final double[][] qMu;

    private final AdmixStates states;
    private final short[][] refPanel;
    private final byte[][] nMismatches;

    private final double[][] fwd;
    private final double[][] bwd;
    private final double[] fwdSums;
    private final double[] bwdSums;
    private int windowSize;
    private int nWindows;
    private final double[][][] bwdCheckPts; // [checkpoint][anc][states]
    private final double[][][] bwdWindow;   // [offset][anc][states]

    private final AdmixHmmUpdater hmmUpdater;
    private final AdmixHmmUpdaterRho hmmUpdaterRho;
    private final AdmixHmmUpdaterT hmmUpdaterT;
    private final double[] ancProbs;

    private final double[][] scaledStateProbs;
    private final double[][] sumStateProbs;

    /**
     * Constructs a new {@code AdmixHmm} instance for the specified data.
     *
     * @param data the input data and analysis parameters
     * @param ibsHaps the object for constructing composite reference
     * haplotypes
     * @throws IllegalArgumentException if
     * {@code data.chromData() != ibsHaps.chromData()}
     * @throws NullPointerException if
     * {@code data == null || ibsHaps == null}
     */
    public AdmixHmm(AdmixData data, IbsHaps ibsHaps) {
        if (data.chromData()!=ibsHaps.chromData()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        this.chromData = data.chromData();
        this.params = data.params();
        this.nMarkers = chromData.targRefGT().nMarkers();
        this.nTargHaps = chromData.nTargHaps();
        this.nAnc = params.fixedParams().nAnc();
        this.nPanels = params.fixedParams().nRefPanels();
        this.minAncProb = params.fixedParams().par().em_anc_prob();
        this.qMu = ParamUtils.qMu(params);

        this.states = new AdmixStates(data, ibsHaps);
        this.refPanel = new short[nMarkers][states.maxStates()];
        this.nMismatches = new byte[nMarkers][states.maxStates()];

        this.fwd = new double[nAnc][states.maxStates()];
        this.bwd = new double[nAnc][states.maxStates()];
        this.fwdSums = new double[nMarkers];
        this.bwdSums = new double[nMarkers];
        this.windowSize = (int) Math.ceil(Math.sqrt(nMarkers));
        this.nWindows = (nMarkers + windowSize - 1)/windowSize;
        this.bwdCheckPts = new double[nWindows][nAnc][states.maxStates()];
        this.bwdWindow = new double[windowSize][nAnc][states.maxStates()];

        this.hmmUpdater = new AdmixHmmUpdater(data);
        this.hmmUpdaterRho = new AdmixHmmUpdaterRho(data);
        this.hmmUpdaterT = new AdmixHmmUpdaterT(data);
        this.ancProbs = new double[nMarkers*nAnc];

        this.scaledStateProbs = new double[nAnc][nPanels];
        this.sumStateProbs = new double[nAnc][nPanels];
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
        return states.ibsHaps().selectedHaps().hapList();
    }

    /**
     * Calculates and stores state data for estimating {@code mu} and {@code T},
     * the ancestry proportions and the number of generations since admixture
     * respectively.
     * @param hapListIndex an index in {@code this.hapList()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.hapList().size()}
     */
    public void runFwdBwdMuT(int hapListIndex) {
        int nHapStates = states.ibsStates(hapListIndex, refPanel, nMismatches);
        setBwdCheckPoints(nHapStates);
        fwdAlgMuT(nHapStates);
    }

    /**
     * Calculates and stores state data for estimating {@code rho} and
     * {@code p}, the pre-admixture, ancestry-specific haplotype switch rates
     * and the probability of reference panel conditional on ancestry
     * respectively.
     * @param hapListIndex an index in {@code this.hapList()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.hapList().size()}
     */
    public void runFwdBwdRhoP(int hapListIndex) {
        int nHapStates = states.ibsStates(hapListIndex, refPanel, nMismatches);
        setBwdCheckPoints(nHapStates);
        fwdAlgRhoP(nHapStates);
    }

    /**
     * Estimates local ancestry probabilities and stores ancestry probabilities
     * in the specified {@code EstimatedAncestry} object.
     * @param targHap a target haplotype index
     * @param estAnc an object for storing the estimated local ancestry
     * @throws IndexOutOfBoundsException if
     * {@code targHap < 0 || targHap >= this.chromData().nTargHaps()}
     * @throws IndexOutOfBoundsException if {@code estAnc.nTargHaps() < targHap}
     * @throws IndexOutOfBoundsException if
     * {@code estAnc.nMarkers()!= this.chromData().targRefGT().nMarkers()}
     * @throws NullPointerException if {@code estAnc == null}
     */
    public void runFwdBwdAnc(int targHap, EstimatedAncestry estAnc) {
        if (targHap<0 || targHap>=nTargHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(targHap));
        }
        int nHapStates = states.ibsStates(targHap, refPanel, nMismatches);
        setBwdCheckPoints(nHapStates);
        fwdAlgAnc(nHapStates);
        estAnc.set(targHap, ancProbs);
    }

    private void fwdAlgRhoP(int nHapStates) {
        double[] ancSums = params.mu();
        initFwd(fwd, nHapStates);
        int window = -1;
        int nextWindowStart = 0;
        int offset = -1;
        for (int m=0; m<nMarkers; ++m)  {
            if (m==nextWindowStart) {
                offset = -1;
                fillBwdWindow(++window, nHapStates);
                nextWindowStart = Math.min(nextWindowStart+windowSize, nMarkers);
            }
            ++offset;
            if (m>0) {
                double[][] bwdM1 = offset==0 ? bwdCheckPts[window-1] : bwdWindow[offset-1];
                hmmUpdaterRho.storeRhoData(m, bwdM1, bwdSums[m-1], bwdWindow[offset],
                        fwd, refPanel[m], nMismatches[m], nHapStates);
            }
            fwdSums[m] = hmmUpdater.fwdUpdate(m, fwd, ancSums, refPanel[m],
                    nMismatches[m], nHapStates);
            setAncAndStateProbs(m, fwd, bwdWindow[offset], nHapStates);
        }
    }

    private void fwdAlgMuT(int nHapStates) {
        double[] ancSums = params.mu();
        initFwd(fwd, nHapStates);
        int window = -1;
        int nextWindowStart = 0;
        int offset = -1;
        for (int m=0; m<nMarkers; ++m)  {
            if (m==nextWindowStart) {
                offset = -1;
                fillBwdWindow(++window, nHapStates);
                nextWindowStart = Math.min(nextWindowStart+windowSize, nMarkers);
            }
            ++offset;
            if (m>0) {
                double[][] bwdM1 = offset==0 ? bwdCheckPts[window-1] : bwdWindow[offset-1];
                hmmUpdaterT.storeTData(m, bwdM1, bwdSums[m-1], bwdWindow[offset],
                        fwd, refPanel[m], nMismatches[m], nHapStates);
            }
            fwdSums[m] = hmmUpdater.fwdUpdate(m, fwd, ancSums, refPanel[m],
                    nMismatches[m], nHapStates);
            setAncAndStateProbs(m, fwd, bwdWindow[offset], nHapStates);
        }
    }

    private void fwdAlgAnc(int nHapStates) {
        double[] ancSums = params.mu();
        initFwd(fwd, nHapStates);
        int window = -1;
        int nextWindowStart = 0;
        int offset = -1;
        for (int m=0; m<nMarkers; ++m)  {
            if (m==nextWindowStart) {
                offset = -1;
                fillBwdWindow(++window, nHapStates);
                nextWindowStart = Math.min(nextWindowStart+windowSize, nMarkers);
            }
            ++offset;
            fwdSums[m] = hmmUpdater.fwdUpdate(m, fwd, ancSums, refPanel[m],
                    nMismatches[m], nHapStates);
            setAncProbs(m, fwd, bwdWindow[offset], nHapStates);
        }
    }

    private void initFwd(double[][] fwd, int nHapStates) {
        short[] panel = refPanel[0];
        for (int i=0; i<nAnc; ++i) {
            for (int h=0; h<nHapStates; ++h) {
                fwd[i][h] = qMu[i][panel[h]];
            }
        }
    }

    private void setBwdCheckPoints(int nHapStates) {
        int windowIndex = nWindows;
        AdmixUtils.fill(bwd, 1.0f/(nAnc*nHapStates));
        bwdSums[nMarkers-1] = 1.0;
        AdmixUtils.copy(bwd, bwdCheckPts[--windowIndex], nHapStates);
        int nextCheckPt = ((nWindows-1)*windowSize) - 1; // checkpoint is last marker in window
        for (int m=nMarkers-2; m>=0; --m) {
            int mP1 = m+1;
            bwdSums[m] = hmmUpdater.bwdUpdate(m, bwd, refPanel[mP1],
                    nMismatches[mP1], nHapStates);
            if (m==nextCheckPt) {
                AdmixUtils.copy(bwd, bwdCheckPts[--windowIndex], nHapStates);
                nextCheckPt -= windowSize;
            }
        }
        assert windowIndex==0;
    }

    private void fillBwdWindow(int window, int nHapStates) {
        // checkpoint is last marker in window
        int m = Math.min((window+1)*windowSize, nMarkers) - 1;
        int lastOffset = m - window*windowSize;
        AdmixUtils.copy(bwdCheckPts[window], bwd, nHapStates);
        AdmixUtils.copy(bwdCheckPts[window], bwdWindow[lastOffset], nHapStates);
        for (int j=(lastOffset-1); j>=0; --j) {
            int mP1 = m;
            --m;
            hmmUpdater.bwdUpdate(m, bwd, refPanel[mP1], nMismatches[mP1],
                    nHapStates);
            AdmixUtils.copy(bwd, bwdWindow[j], nHapStates);
        }
    }

    private void setAncProbs(int m, double[][] fwd, double[][] bwd0,
            int nHapStates) {
        int offset = m*nAnc;
        double cumSum = 0f;
        for (int i=0; i<nAnc; ++i) {
            double[] f = fwd[i];
            double[] b = bwd0[i];
            double sum = 0f;
            for (int h=0; h<nHapStates; ++h) {
                sum += f[h]*b[h];
            }
            ancProbs[offset + i] = sum;
            cumSum += sum;
        }
        AdmixUtils.scale(ancProbs, offset, (offset + nAnc), (1.0/cumSum));
    }

    private void setAncAndStateProbs(int m, double[][] fwd, double[][] bwd0,
            int nHapStates) {
        int offset = m*nAnc;
        double cumSum = 0f;
        short[] panel = refPanel[m];
        double[] ancSums = new double[nAnc];
        for (int i=0; i<nAnc; ++i) {
            double[] f = fwd[i];
            double[] b = bwd0[i];
            double ancSum = 0f;
            for (int h=0; h<nHapStates; ++h) {
                int j = panel[h];
                double scaledStateProb = f[h]*b[h];
                ancSum += scaledStateProb;
                scaledStateProbs[i][j] += scaledStateProb;
            }
            ancProbs[offset + i] = ancSum;
            ancSums[i] = ancSum;
            cumSum += ancSum;
        }
        AdmixUtils.scale(ancProbs, offset, (offset + nAnc), (1.0/cumSum));
        updateSumStateProbs(cumSum, ancSums);
    }

    private void updateSumStateProbs(double cumSum, double[] ancSums) {
        double unscaleFactor = 1.0/cumSum;
        for (int i=0; i<nAnc; ++i) {
            ancSums[i] *= unscaleFactor;
	    if (ancSums[i]>=minAncProb) {
		for (int j=0; j<nPanels; ++j) {
		    sumStateProbs[i][j] += unscaleFactor*scaledStateProbs[i][j];
		}
	    }
            Arrays.fill(scaledStateProbs[i], 0.0);
        }
    }

    /**
     * Add the data for estimating {@code mu} and {@code T} to the specified
     * {@code ParamEstimateData} object, where {@code mu} and {@code T}
     * are the ancestry proportions and the number of generations
     * since admixture respectively.
     * @param paramData the object which estimates analysis parameters
     * @throws NullPointerException if {@code paramData == null}
     */
    public void updateMuT(ParamEstimateData paramData) {
        paramData.addTSwitchData(hmmUpdaterT.sumTSwitchProbs(),
                hmmUpdaterT.sumTGenDist());
        hmmUpdaterT.clear();
        updateStateProbs(paramData);
    }

    /**
     * Add the data for estimating {@code rho} and {@code p} to the specified
     * {@code ParamEstimateData} object, where {@code rho} and {@code p}
     * are the pre-admixture, ancestry-specific haplotype switch rates
     * and the probability of reference panel conditional on ancestry
     * respectively.
     * @param paramData the object which estimates analysis parameters
     * @throws NullPointerException if {@code paramData == null}
     */
    public void updateRhoP(ParamEstimateData paramData) {
        for (int i=0; i<nAnc; ++i) {
            paramData.addRhoSwitchData(i, hmmUpdaterRho.sumRhoSwitchProbs(i),
                    hmmUpdaterRho.sumRhoGenDist(i));
        }
        hmmUpdaterRho.clear();
        updateStateProbs(paramData);
    }

    private void updateStateProbs(ParamEstimateData paramData) {
        for (int i=0; i<nAnc; ++i) {
            for (int j=0; j<nPanels; ++j) {
                paramData.addStateData(i, j, sumStateProbs[i][j]);
            }
        }
        for (int i=0; i<nAnc; ++i) {
            for (int j=0; j<nPanels; ++j) {
                sumStateProbs[i][j] = 0.0;
            }
        }
    }
}

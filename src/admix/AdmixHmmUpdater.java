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

/**
 * <p>Class {@code AdmixHmmUpdater} has methods for one-step
 * updates of forward and backward HMM probabilities.</p>
 *
 * <p>Instances of {@code AdmixHmmUpdater} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixHmmUpdater {

    private final int nMarkers;
    private final int nAnc;
    private final ParamsInterface params;
    private final AdmixHmmProbs hmm;
    private final double[] iShift;
    private final double[] jShift;

    private final double[] mu;
    private final double[][] q;
    private final double[][][] pObserved;

    /**
     * Constructs a new {@code AdmixHmmUpdater} instance from the specified
     * data.
     * @param data the input data for local ancestry inference on a chromosome
     * @throws NullPointerException if {@code data == null}
     */
    public AdmixHmmUpdater(AdmixData data) {
        AdmixChromData chromData = data.chromData();
        FixedParams fixedParams = data.params().fixedParams();
        this.params = data.params();
        this.hmm = data.hmmProbs();
        this.nMarkers = chromData.targRefGT().nMarkers();
        this.nAnc = fixedParams.nAnc();

        this.iShift = new double[nAnc];
        this.jShift = new double[fixedParams.nRefPanels()];

        this.mu = params.mu();
        this.q = ParamUtils.q(params);
        this.pObserved = ParamUtils.pObserved(params);
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return nMarkers;
    }

    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return nAnc;
    }

    /**
     * Updates the forward state and ancestry probabilities, and scales both of
     * these probabilities to sum to 1.0f.  The unscaled sum of the updated
     * forward probabilities will be returned.  The contract for this method
     * is unspecified if the specified state probabilities do not sum to 1.0f
     * or if the specified ancestry probabilities do not sum to 1.0f.
     * @param m the marker index
     * @param fwd the scaled forward state probabilities at marker
     * {@code (m - 1)} which will be updated
     * @param ancProbs the scaled ancestry probabilities
     * @param hap2Panel a map from HMM reference haplotype to reference panel
     * @param mismatch the number of allele mismatches (0 or 1) for each HMM
     * state at marker {@code m}
     * @param nStates the number of HMM states
     * @return the sum of the calculated forward probabilities before the
     * probabilities are scaled to sum to 1.0
     * @throws IndexOutOfBoundsException if
     * {@code m < 0 || m >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code fwd.length < this.nAnc() || ancSums.length < this.nAnc()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2Panel.length < nStates || mismatch.length < nStates}
     * @throws IndexOutOfBoundsException if there exists {@code i} such that
     * {@code (0 <= i) && (i < this.nAnc()) && (fwd[i].length < nStates)}
     * @throws NullPointerException if any array is null
     */
    public double fwdUpdate(int m, double[][] fwd, double[] ancProbs,
            short[] hap2Panel, byte[] mismatch, int nStates) {
        double sumFwdProb = 0.0;
        for (int i=0; i<nAnc; ++i) {
            double scale = hmm.pNoRecTNoRecRho(m,i);
            setJShifts(m, i, ancProbs[i]);
            ancProbs[i] = 0f;
            double[] fwdi = fwd[i];
            for (int h=0; h<nStates; ++h) {
                int j = hap2Panel[h];
                double em = pObserved[i][j][mismatch[h]];
                fwdi[h] = em*(scale*fwdi[h] + jShift[j]);
                ancProbs[i] += fwdi[h];
            }
            sumFwdProb += ancProbs[i];
        }
        double invSumFwdProb = 1.0/sumFwdProb;
        AdmixUtils.scale(fwd, invSumFwdProb);
        AdmixUtils.scale(ancProbs, invSumFwdProb);
        return sumFwdProb;
    }

    private void setJShifts(int m, int i, double ancSum) {
        // computation of shift[j] assumes that sum(anc)==1.0
        // m = marker, i = ancestry, j = reference panel
        for (int j=0; j<jShift.length; ++j) {
            jShift[j] = hmm.pRecTqMu(m, i, j) + hmm.pNoRecTRecRho(m, i)*q[i][j]*ancSum;
        }
    }

    /**
     * Updates the backward state probabilities and scales them to sum to 1.0.
     * Returns the sum of the updated backward state probabilities from they
     * are scaled to sum to 1.0.
     * @param m the marker index
     * @param bwd the backward state probabilities for marker {@code (m + 1)}
     * @param hap2Panel a map from reference haplotype to reference panel at
     * marker {@code (m + 1)}
     * @param mismatch the number of allele mismatches (0 or 1) for each HMM
     * state at marker {@code (m + 1)}
     * @param nStates the number of HMM states
     * @return the sum of the calculated backward probabilities before
     * scaling them to sum to 1.0
     * @throws IndexOutOfBoundsException if
     * {@code m < 0 || m >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code bwd.length < this.nAnc()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2Panel.length < nStates || mismatch.length < nStates}
     * @throws IndexOutOfBoundsException if there exists {@code i} such that
     * {@code (0 <= i) && (i < this.nAnc()) && (bwd[i].length < nStates)}
     * @throws NullPointerException if any array is null
     */
    public double bwdUpdate(int m, double[][] bwd, short[] hap2Panel,
            byte[] mismatch, int nStates) {
        double bwdShiftSum = 0.0;
        for (int i=0; i<bwd.length; ++i) {
            iShift[i] = 0.0;
            for (int h=0; h<nStates; ++h) {
                int j = hap2Panel[h];
                double em =  pObserved[i][j][mismatch[h]];
                bwd[i][h] *= em;
                iShift[i] += bwd[i][h]*q[i][j];
            }
            bwdShiftSum += iShift[i]*mu[i];
        }
        int mP1 = m+1;
        double bwdProbSum = 0.0;
        for (int i=0; i<bwd.length; ++i) {
            double shift = hmm.pRecT(mP1)*bwdShiftSum + hmm.pNoRecTRecRho(mP1, i)*iShift[i];
            double scale = hmm.pNoRecTNoRecRho(mP1, i);
            for (int h=0; h<nStates; ++h) {
                bwd[i][h] = scale*bwd[i][h] + shift;
                bwdProbSum += bwd[i][h];
            }
        }
        AdmixUtils.scale(bwd, 1.0/bwdProbSum);
        return bwdProbSum;
    }
}

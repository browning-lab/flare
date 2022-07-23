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

import vcf.MarkerMap;

/**
 * <p>Class {@code AdmixHmmUpdaterT} computes and stores data for estimating
 * {@code T}, the number of generations since admixture.</p>
 *
 * <p>Instances of {@code AdmixHmmUpdateT} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixHmmUpdaterT {

    private final MarkerMap map;
    private final int nMarkers;
    private final int nAnc;
    private final ParamsInterface params;
    private final AdmixHmmProbs hmm;

    private double sumTSwitchProbs;
    private double sumTGenDist;

    private final double[] inv1Mmu;
    private final double[][][] pObserved;

    /**
     * Constructs a new {@code AdmixHmmUpdaterT} instance from the specified
     * data.
     * @param data the input data for local ancestry inference on a chromosome
     * @throws NullPointerException if {@code data == null}
     */
    public AdmixHmmUpdaterT(AdmixData data) {
        AdmixChromData chromData = data.chromData();
        FixedParams fixedParams = data.params().fixedParams();
        this.map = chromData.map();
        this.params = data.params();
        this.hmm = data.hmmProbs();
        this.nMarkers = chromData.targRefGT().nMarkers();
        this.nAnc = fixedParams.nAnc();

        this.sumTSwitchProbs = 0.0;
        this.sumTGenDist = 0.0;

        this.inv1Mmu = ParamUtils.inv1Mmu(params);
        this.pObserved = ParamUtils.pObserved(params);
    }

    /**
     * Clears the stored data so that the {@code this.sumTSwitchProbs()} and
     * {@code this.sumTGenDist()} methods return 0.0;
     */
    public void clear() {
        sumTSwitchProbs = 0.0;
        sumTGenDist = 0.0;
    }

    /**
     * Returns the stored sum of post-admixture switch probabilities.
     * @return the stored sum of post-admixture switch probabilities
     */
    public double sumTSwitchProbs() {
        return sumTSwitchProbs;
    }

    /**
     * Returns the stored sum of genetic distances.
     * @return the stored sum of genetic distances
     */
    public double sumTGenDist() {
        return sumTGenDist;
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
     * Computes and stores data for estimating {@code T}, the number of
     * generations since admixture.
     * The contract for this method is unspecified if any of the specified
     * forward and backward state probabilities do not sum to 1.0f.
     * @param m the marker index
     * @param bwdM1 scaled backward state probabilities for marker {@code (m - 1)}
     * @param bwdM1Sum the sum of the backward state probabilities for marker
     * {@code (m - 1)}
     * @param bwd scaled backward state probabilities for marker {@code m}
     * @param fwdM1 the scaled forward state probabilities at marker
     * {@code (m - 1)}
     * @param hap2Panel the panel containing each haplotype
     * @param mismatch the number of allele mismatches (0 or 1) for each HMM
     * state at marker {@code m}
     * @param nStates the number of HMM states
     * @throws IndexOutOfBoundsException if
     * {@code (m < 0) || (m >= this.nMarkers())}
     * @throws IndexOutOfBoundsException if
     * {@code (bwd.length < this.nAnc()) || (bwdP1.length < this.nAnc())
     * || (fwd.length < this.nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (hap2Panel.length < nStates) || (mismatch.length < nStates)}
     * @throws IndexOutOfBoundsException if there exists {@code i} such that
     * {@code (0 <= i) && (i < this.nAnc()
     * && ((bwd[i].length < nStates) || (bwdP1[i].length < nStates)
     * || (fwd[i].length < nStates))}
     * @throws NullPointerException if any array is null
     */
    public void storeTData(int m, double[][] bwdM1, double bwdM1Sum, double[][] bwd,
            double[][] fwdM1, short[] hap2Panel, byte[] mismatch, int nStates) {
        assert m>0;
        double morgans = 0.01*map.genDist().get(m);
        double sumStateProbs = 0.0;
        double tEst = 0.0;
        for (int i=0; i<nAnc; ++i) {
            double sumAncStateProbs = 0.0;
            double sumAncFwdM1 = 0.0;
            double noRecSum = 0.0;
            double recSum = 0.0;
            for (int h=0; h<nStates; ++h) {
                int j = hap2Panel[h];
                sumAncFwdM1 += fwdM1[i][h];
                // Multiply bwdM1 by bwdM1Sum so that bwd and bwdM1 are on same scale
                double stateProb = bwdM1Sum*bwdM1[i][h]*fwdM1[i][h];
                sumAncStateProbs += stateProb;
                double bwdEm = bwd[i][h]*pObserved[i][j][mismatch[h]];
                recSum += bwdEm*hmm.pHapChange(m, i, j);
                noRecSum += bwdEm*hmm.pNoRecTNoRecRho(m, i)*fwdM1[i][h];
            }
            double tNum = sumAncStateProbs - recSum*sumAncFwdM1 - noRecSum;
            tEst += Math.max(0.0, tNum)*inv1Mmu[i];
            sumStateProbs += sumAncStateProbs;
        }
        double invSumStateProbs = 1.0/sumStateProbs;
        sumTSwitchProbs += invSumStateProbs*tEst;
        sumTGenDist += morgans;
    }
}

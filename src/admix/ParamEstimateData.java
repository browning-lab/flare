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

import blbutil.Const;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * <p>Class {@code ParamEstimateData} estimates the analysis parameters for a
 * local ancestry inference analysis.</p>
 *
 * <p>Instances of class {@code ParamEstimateData} are thread-safe if
 * the requirements in each method's documentation is satisfied.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ParamEstimateData {

    private final ParamsInterface prevParams;
    private final int nAnc;
    private final int nPanels;
    private final DoubleAdder[][] stateProbs;
    private final DoubleAdder[] sumRhoSwitchProbs;
    private final DoubleAdder[] sumRhoGenDist;
    private final DoubleAdder sumTSwitchProbs;
    private final DoubleAdder sumTGenDist;

    /**
     * Constructs a new {@code ParamsEstimateData} instance for the specified
     * data.
     * @param params the current parameter values
     * @throws NullPointerException if {@code params == null}
     */
    public ParamEstimateData(ParamsInterface params) {
        this.prevParams = params;
        this.nAnc = params.fixedParams().nAnc();
        this.nPanels = params.fixedParams().nRefPanels();
        this.stateProbs = createDoubleAdderArray(nAnc, nPanels);
        this.sumRhoSwitchProbs = createDoubleAdderArray(nAnc);
        this.sumRhoGenDist = createDoubleAdderArray(nAnc);
        this.sumTSwitchProbs = new DoubleAdder();
        this.sumTGenDist = new DoubleAdder();
    }

    private static DoubleAdder[][] createDoubleAdderArray(int nRows, int nCols) {
        DoubleAdder[][] daa = new DoubleAdder[nRows][];
        for (int r=0; r<nRows; ++r) {
            daa[r] = createDoubleAdderArray(nCols);
        }
        return daa;
    }

    private static DoubleAdder[] createDoubleAdderArray(int length) {
        DoubleAdder[] daa = new DoubleAdder[length];
        for (int j=0; j<daa.length; ++j) {
            daa[j] = new DoubleAdder();
        }
        return daa;
    }

    private static double sum(DoubleAdder[] daa) {
        int sum = 0;
        for (DoubleAdder da : daa) {
            sum += da.sum();
        }
        return sum;
    }

    /**
     * Returns the analysis parameters used for parameter estimation.
     * @return the analysis parameters used for parameter estimation
     */
    public ParamsInterface prevParams() {
        return prevParams;
    }

    /**
     * Returns the estimated analysis parameters. Invocation in the
     * absence of concurrent calls to {@code this.addMaxAncData()},
     * {@code this.addStateData()}, {@code this.addRhoSwitchData()}, and
     * {@code this.TSwitchaddData()}, returns an accurate result, but an
     * accurate result is not guaranteed in the presence of concurrent updates.
     * @param initParams the initial parameters
     * @return the estimated analysis parameters
     */
    public ParamsInterface estimatedParams(ParamsInterface initParams) {
        return new EstimatedParams(this, initParams);
    }

    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return nAnc;
    }

    /**
     * Returns the number of reference panels.
     * @return the number of reference panels
     */
    public int nRefPanels() {
        return nPanels;
    }

    /**
     * Adds the specified sum of probabilities for a set of states to
     * {@code this}.
     * @param anc the ancestry
     * @param panel the reference panel
     * @param stateProb the sum of probabilities for a set of states
     * @throws IndexOutOfBoundsException if
     * {@code anc < 0 || anc >= this.nAnc()}
     * @throws IndexOutOfBoundsException if
     * {@code panel < 0 || panel >= this.nPanels()}
     */
    public void addStateData(int anc, int panel, double stateProb) {
        stateProbs[anc][panel].add(stateProb);
    }

    /**
     * Add the specified data to {@code this}.
     * @param anc the ancestry
     * @param rhoSwitchProb a sum of pre-admixture switch probabilities
     * @param rhoGenDist a sum of weighted genetic distances
     * @throws IndexOutOfBoundsException if
     * {@code anc < 0 || anc >= this.nAnc()}
     */
    public void addRhoSwitchData(int anc, double rhoSwitchProb, double rhoGenDist) {
        sumRhoSwitchProbs[anc].add(rhoSwitchProb);
        sumRhoGenDist[anc].add(rhoGenDist);
    }

    /**
     * Add the specified data to {@code this}.
     * @param tSwitchProb a sum of post-admixture switch probabilities
     * @param tGenDist a sum of genetic distances
     */
    public void addTSwitchData(double tSwitchProb, double tGenDist) {
        sumTSwitchProbs.add(tSwitchProb);
        sumTGenDist.add(tGenDist);
    }

    /**
     * Returns a list with length {@code this.nAnc()} whose {@code i}-th entry
     * is the estimated proportion of target samples having ancestry {@code i}.
     * The returned value is not an atomic snapshot. Invocation in the
     * absence of concurrent calls to {@code this.addStateData()} returns an
     * accurate result, but an accurate result is not guaranteed in
     * the presence of concurrent updates.
     * @return a list of the estimated ancestry proportions
     */
    public double[] mu() {
        double[] mu = new double[nAnc];
        double sum = 0.0;
        for (int i=0; i<nAnc; ++i) {
            mu[i] = sum(stateProbs[i]);
            sum += mu[i];
        }
        if (sum==0.0) {
            return prevParams.mu();
        }
        else {
            AdmixUtils.scale(mu, 1.0/sum);
            return mu;
        }
    }

    /**
     * Returns a list with length {@code this.nAnc()} whose {@code i}-th element
     * is the estimated pre-admixture switch rate for ancestry {@code i}.
     * The returned value is not an atomic snapshot. Invocation in the
     * absence of concurrent calls to {@code this.addRhoSwitchData()} returns
     * an accurate result, but an accurate result is not guaranteed in the
     * presence of concurrent udpates.
     * @return a list of the pre-admixture switch rates
     */
    public double[] rho() {
        double[] rho = new double[nAnc];
        for (int i=0; i<nAnc; ++i) {
            double val = sumRhoSwitchProbs[i].sum()/sumRhoGenDist[i].sum();
            rho[i] = (Double.isFinite(val) && val>0.0) ? val : prevParams.rho(i);
        }
        return rho;
    }

    /**
     * Returns an array with {@code this.nAnc()} rows and
     * {@code this.nRefPanels()} columns whose {@code [i][j]}-th entry is the
     * estimated probability that the HMM state with ancestry {@code i} has
     * reference haplotype from reference panel {@code j}.
     * The returned value is not an atomic snapshot. Invocation in the
     * absence of concurrent calls to {@code this.addStateData()} returns an
     * accurate result, but concurrent updated may not be included in the
     * returned value.
     * @return an array whose {@code [i][j]}-th entry is the estimated
     * probability that the HMM state with ancestry {@code i} has reference
     * haplotype from reference panel {@code j}
     */
    public double[][] p() {
        double[][] da = new double[nAnc][nPanels];
        for (int i=0; i<nAnc; ++i) {
            double[] row = da[i];
            double rowSum = 0.0;
            for (int j=0; j<nPanels; ++j) {
                row[j] = stateProbs[i][j].sum();
                rowSum += row[j];
            }
            if (rowSum==0.0) {
                for (int j=0; j<nPanels; ++j) {
                    row[j] = prevParams.p(i, j);
                }
            }
            else {
                AdmixUtils.scale(row, 1.0/rowSum);
            }
        }
        return da;
    }

    /**
     * Returns the estimated number of generations since admixture.
     * The returned value is not an atomic snapshot. Invocation in the
     * absence of concurrent calls to {@code this.addTSwitchData()} returns
     * an accurate result, but an accurate result is not guaranteed in
     * the presence of concurrent updates.
     * @return the estimated number of generations since admixture
     */
    public double t() {
        double value = sumTSwitchProbs.sum() / sumTGenDist.sum();
        return (Double.isFinite(value) && value>0.0) ? value : prevParams.T();
    }

    /**
     * Prints a string representation of {@code this}.  The exact details of
     * the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(1<<9);
        sb.append("# mu[ancestry]");
        sb.append(Const.nl);
        sb.append(Const.nl);
        AdmixUtils.appendLineWithVector(sb, mu());
        sb.append(Const.nl);
        sb.append("# p[ancestry][panel]");
        sb.append(Const.nl);
        sb.append(Const.nl);
        AdmixUtils.appendLinesWithMatrix(sb, p());
        sb.append(Const.nl);
        sb.append("# theta[ancestry][panel]");
        sb.append(Const.nl);
        sb.append(Const.nl);
        AdmixUtils.appendLinesWithMatrix(sb, prevParams.theta());
        sb.append(Const.nl);
        sb.append("# rho[ancestry]");
        sb.append(Const.nl);
        sb.append(Const.nl);
        AdmixUtils.appendLineWithVector(sb, rho());
        return sb.toString();
    }
}

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

import ints.IntArray;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code DefaultParams} represents default analysis parameters for a
 * local ancestry inference analysis.</p>
 *
 * <p>Instances of class {@code DefaultParams} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class DefaultParams implements ParamsInterface {

    private final FixedParams fixedParams;
    private final double t;                    // gen since admixture
    private final double[] mu;                 // ancestry proportions
    private final double[][] theta;            // miscopy probability matrix
    private final double[][] p;                // copying probability matrix
    private final double[] rho;                // population switch probabilities

    /**
     * Constructs a new {@code DefaultParams} instance for the specified
     * data.
     * @param fixedParams the fixed parameters for a local ancestry inference
     * analysis
     * @throws NullPointerException if {@code fixedParams == null}
     */
    public DefaultParams(FixedParams fixedParams) {
        this.fixedParams = fixedParams;
        this.t = fixedParams.par().gen();
        this.mu = initMu(fixedParams.nAnc());
        this.p = defaultP(fixedParams);
        this.theta = defaultTheta(fixedParams);
        double initRho = ((double) 200_000)/fixedParams.nRefHaps();
        this.rho = IntStream.range(0, fixedParams.nAnc())
                .mapToDouble(j -> initRho)
                .toArray();
    }

    private static double[] initMu(int nAnc) {
        double[] mu = new double[nAnc];
        Arrays.fill(mu, 1.0f/nAnc);
        return mu;
    }

    /**
     * Returns the default {@code theta} values. The returned array will have
     * one row per ancestry and one column per reference panel.
     * @param fixedParams the fixed parameters for a local ancestry inference
     * analysis
     * @return the default {@code theta} values
     * @throws NullPointerException if {@code fixedParams == null}
     */
    public static double[][] defaultTheta(FixedParams fixedParams) {
        int nAnc = fixedParams.nAnc();
        int nRefHaps = fixedParams.nRefHaps();
        int nRefPanels = fixedParams.nRefPanels();
        double theta0 = liStephensPMismatch(nRefHaps);
        double[] row = IntStream.range(0, nRefPanels)
                .mapToDouble(j -> theta0)
                .toArray();
        return IntStream.range(0, nAnc)
                .mapToObj(j -> row)
                .toArray(double[][]::new);
    }

    /**
     * <p>Return an approximation to the allele mismatch probability suggested
     * by Li and Stephens.  The approximation uses a Riemann sum approximation
     * of the natural log function.</p>
     *
     * <p>Refs:
     * Li N, Stephens M. Genetics 2003 Dec;165(4):2213-33 and
     * Marchini J, Howie B. Myers S, McVean G, Donnelly P. 2007;39(7):906-13.</p>
     *
     * @param nHaps the number of haplotypes
     * @return an approximation to the Li and Stephens allele mismatch
     * probability
     * @throws IllegalArgumentException if {@code nHaps < 1}
     */
    private static double liStephensPMismatch(int nHaps) {
        if (nHaps<1) {
            throw new IllegalArgumentException(String.valueOf(nHaps));
        }
        double theta = 1.0/((Math.log(nHaps) + 0.5));
        return theta/(2.0*(theta + nHaps));
    }

    /**
     * Returns the default {@code p} values.  The returned array will have
     * one row per ancestry and one column per reference panel.
     * @param fixedParams the fixed parameters for a local ancestry inference
     * analysis
     * @return the default {@code p} values
     * @throws NullPointerException if {@code fixedParams == null}
     */
    public static double[][] defaultP(FixedParams fixedParams) {
        float panelWeight = fixedParams.par().panel_weight();
        int nAnc = fixedParams.nAnc();
        int nRefPanels = fixedParams.nRefPanels();
        double[][] p = new double[nAnc][nRefPanels];
        for (int i=0; i<nAnc; ++i) {
            IntArray ancIPanels = fixedParams.ancPanels(i);
            int nRelevant = ancIPanels.size();
            double ancPanelsProb = nRelevant==nRefPanels ? 1.0 : panelWeight;
            double ancPanelValue = ancPanelsProb/nRelevant;
            if (nRelevant<nRefPanels) {
                double otherValue = (1.0 - ancPanelsProb)/(nRefPanels - nRelevant);
                Arrays.fill(p[i], otherValue);
            }
            for (int k=0; k<nRelevant; ++k) {
                p[i][ancIPanels.get(k)] = ancPanelValue;
            }
        }
        return p;
    }

    @Override
    public FixedParams fixedParams() {
        return fixedParams;
    }

    @Override
    public double T() {
        return t;
    }

    @Override
    public double[] mu() {
        return mu.clone();
    }

    @Override
    public double rho(int i) {
        return rho[i];
    }

    @Override
    public double[] rho() {
        return rho.clone();
    }

    @Override
    public double p(int i, int j) {
        return p[i][j];
    }

    @Override
    public double[][] p() {
        return AdmixUtils.cloneArray(p);
    }

    @Override
    public double[][] theta() {
        return AdmixUtils.cloneArray(theta);
    }

    @Override
    public String toString() {
        return ParamUtils.toString(this);
    }
}

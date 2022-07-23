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
 * <p>Class {@code EstimatedParams} stores estimated parameters for a
 * local ancestry inference analysis.</p>
 *
 * <p>Instances of class {@code EstimatedtParams} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class EstimatedParams implements ParamsInterface {

    private final FixedParams fixedParams;
    private final double t;                    // gen since admixture
    private final double[] mu;                 // ancestry proportions
    private final double[][] theta;            // miscopy probability matrix
    private final double[][] p;                // copying probability matrix
    private final double[] rho;                // population switch probabilities

    /**
     * Constructs a new {@code EstimatedParams} instance with the specified
     * data. The {@code theta} and {@code minP} parameters are set equal
     * to {@code initParams.theta()} and {@code initParams.minP()} respectively.
     * Other parameters are obtained from the specified {@code estimatedParams}
     * object.  The specified {@code ParamEstimateData} object is not modified
     * by this method.
     * @param estimatedParams the estimated parameter values
     * @param initParams the initial parameters
     * @throws NullPointerException if
     * {@code (estimatedParams == null) || (initParams == null)}
     */
    public EstimatedParams(ParamEstimateData estimatedParams,
            ParamsInterface initParams) {
        this.fixedParams = initParams.fixedParams();
        this.t = estimatedParams.t();
        this.mu = estimatedParams.mu();
        this.p = initParams.p();
        this.theta = initParams.theta();
        this.rho = initParams.rho();
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

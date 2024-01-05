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

/**
 * <p>Class {@code PartiallyUpdatedParams} represents partially updated
 * analysis parameters for a local ancestry inference analysis.
 *
 * <p>Instances of class {@code PartiallyUpdatedParams} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PartiallyUpdatedParams implements ParamsInterface {

    private final FixedParams fixedParams;
    private final double t;                    // gen since admixture
    private final double[] mu;                 // ancestry proportions
    private final double[][] theta;            // miscopy probability matrix
    private final double[][] p;                // copying probability matrix
    private final double[] rho;                // population switch probabilities

    /**
     * Constructs a new {@code PartiallyUpdatedParams} instance for the
     * specified data. The constructed {@code PartiallyUpdatedParams} object
     * will have the same parameter values as the specified {@code baseParams}
     * object, except for its {@code p()} and {@code rho()} parameters.
     * The {@code i}-th elements the constructed object's {@code p()} and
     * {@code rho()} parameters are obtained from {@code paramEstimates[i]}.
     * @param paramEstimates a list with estimated parameter values
     * for each ancestry
     * @param baseParams the base parameter values
     * @throws IllegalArgumentException if
     * {@code paramEstimates.length != baseParams.fixedParams().nAnc()}
     * @throws IllegalArgumentException if there exists {@code i} such that
     * {@code (0 <= i) && (i < baseParams.fixedParams().nAnc())
     * && (paramEstimates[i].prevParams().fixedParams() != baseParams.fixedParams())}
     * @throws NullPointerException if
     * {@code (paramEstimates == null) || (baseParams == null)}
     * @throws NullPointerException if there exists {@code i} such that
     * {@code (0 <= i) && (i < baseParams.fixedParams().nAnc())
     * && (paramEstimates[i] == null)}
     */
    public PartiallyUpdatedParams(ParamEstimateData[] paramEstimates,
            ParamsInterface baseParams) {
        checkParameters(paramEstimates, baseParams);
        ParamsInterface params = baseParams;
        this.fixedParams = params.fixedParams();
        this.t = baseParams.T();
        this.mu = baseParams.studyMu();
        this.p = updatedP(fixedParams, paramEstimates);
        this.theta = baseParams.theta();
        this.rho = updatedRho(fixedParams, paramEstimates);
    }

    private static  void checkParameters(ParamEstimateData[] paramEstimates,
            ParamsInterface baseParams) {
        FixedParams fixedParams = baseParams.fixedParams();
        int nAnc = fixedParams.nAnc();
        if (paramEstimates.length != nAnc) {
            throw new IllegalArgumentException(String.valueOf(paramEstimates.length));
        }
        for (int i=0; i<paramEstimates.length; ++i) {
            if (paramEstimates[i].prevParams().fixedParams() != fixedParams) {
                String msg = "inconsistent parameters for ancestry " + i;
                throw new IllegalArgumentException(msg);
            }
        }
    }

    private static double[][] updatedP(FixedParams fixedParams,
            ParamEstimateData[] estimatedParams) {
        double[][] updatedP = new double[fixedParams.nAnc()][fixedParams.nRefPanels()];
        for (int i=0; i<updatedP.length; ++i) {
            updatedP[i] = estimatedParams[i].p()[i];
        }
        return updatedP;
    }


    private static double[] updatedRho(FixedParams fixedParams,
            ParamEstimateData[] estimatedParams) {
        double[] updatedRho = new double[fixedParams.nAnc()];
        for (int i=0; i<updatedRho.length; ++i) {
            updatedRho[i] = estimatedParams[i].rho()[i];
        }
        return updatedRho;
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
    public double[] studyMu() {
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

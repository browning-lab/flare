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

import java.util.stream.IntStream;

/**
 * <p>Class {@code PartiallyUpdatedParams} represents partially updated
 * analysis parameters for a local ancestry inference analysis.
 *
 * <p>Instances of class {@code PartiallyUpdatedParams} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PartiallyUpdatedParams implements ParamsInterface {

    private final SampleData sampleData;
    private final double t;                    // gen since admixture
    private final double[] mu;                 // ancestry proportions
    private final double[][] theta;            // miscopy probability matrix
    private final double[][] p;                // copying probability matrix
    private final double[] rho;                // population switch probabilities

    /**
     * Constructs a new {@code PartiallyUpdatedParams} instance for the
     * specified data. The {@code i}-th element of the constructed
     * object's {@code p()} and {@code rho()} parameters are obtained
     * from {@code paramData[i]}.  The constructed object's {@code T()},
     * {@code mu()}, and {@code theta()} parameters are obtained from
     * {@code sampleData.par().gen()},
     * {@code ParamsInterface.defaultMu(sampleData.nAnc())},
     * and {@code ParamsInterface.defaultTheta(sampleData)} respectively.
     * @param sampleData reference and target sample metadata
     * @param paramData a list with estimated parameter values
     * for each ancestry
     * @throws IllegalArgumentException if
     * {@code (paramData.length != sampleData.nAnc())}
     * @throws IllegalArgumentException if there exists {@code i} such that
     * {@code ((0 <= i) && (i < paramData.length)
     * && (paramData[i].prevParams().sampleData() != sampleData))}
     * @throws NullPointerException if
     * {@code ((sampleData == null) || (paramData == null)) }
     * @throws NullPointerException if there exists {@code i} such that
     * {@code ((0 <= i) && (i < paramData.length) && (paramData[i] == null))}
     */
    public PartiallyUpdatedParams(SampleData sampleData,
            ParamEstimateData[] paramData) {
        checkParameters(sampleData, paramData);
        int nAnc = sampleData.nAnc();
        this.sampleData = sampleData;
        this.t = sampleData.par().gen();
        this.mu = ParamsInterface.defaultMu(nAnc);
        this.theta = ParamsInterface.defaultTheta(sampleData);
        this.p = IntStream.range(0, nAnc)
                .mapToObj(i -> paramData[i].p()[i])
                .toArray(double[][]::new);
        this.rho = IntStream.range(0, nAnc)
                .mapToDouble(i -> paramData[i].rho()[i])
                .toArray();

    }

    private static  void checkParameters(SampleData sampleData,
            ParamEstimateData[] paramData) {
        int nAnc = sampleData.nAnc();
        if (paramData.length != nAnc) {
            throw new IllegalArgumentException(String.valueOf(paramData.length));
        }
        for (int i=0; i<paramData.length; ++i) {
            if (paramData[i].prevParams().sampleData() != sampleData) {
                String msg = "inconsistent parameters for ancestry " + i;
                throw new IllegalArgumentException(msg);
            }
        }
    }

    @Override
    public SampleData sampleData() {
        return sampleData;
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

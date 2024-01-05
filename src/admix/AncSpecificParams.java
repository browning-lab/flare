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

import ints.IntArray;
import java.util.stream.IntStream;

/**
 * <p>Class {@code AncSpecificParams} represents the analysis parameters for a
 * local ancestry inference analyses that is used to estimate the copying
 * probabilities {@code p} and the population switch probabilities {@code rho}
 * for a fixed ancestry.</p>
 *
 * <p>Instances of class {@code AncSpecificParams} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AncSpecificParams implements ParamsInterface {

    private final FixedParams fixedParams;
    private final double t;                    // gen since admixture
    private final double[] mu;                 // ancestry proportions
    private final double[][] theta;            // miscopy probability matrix
    private final double[][] p;                // copying probability matrix
    private final double[] rho;                // population switch probabilities

    /**
     * Constructs a new {@code AncSpecificParams} instance from the specified
     * data.
     * @param fixedParams the fixed parameters for a local ancestry inference
     * analysis
     * @param anc the ancestry index
     * @throws IllegalArgumentException if
     * {@code anc < 0 || anc >= fixedParams.nAnc()}
     * @throws NullPointerException if {@code fixedParams == null}
     */
    public AncSpecificParams(FixedParams fixedParams, int anc) {
        if (anc<0 || anc>=fixedParams.nAnc()) {
            throw new IllegalArgumentException(String.valueOf(anc));
        }
        this.fixedParams = fixedParams;
        this.t = Double.MIN_VALUE;
        this.mu = IntStream.range(0, fixedParams.nAnc())
                .mapToDouble(j -> (j==anc ? 1.0 :  0.0))
                .toArray();
        double[] panelToProp = panelToProportion(fixedParams);
        this.p = IntStream.range(0, fixedParams.nAnc())
                .mapToObj(j -> panelToProp)
                .toArray(double[][]::new);
        this.theta = DefaultParams.defaultTheta(fixedParams);
        double initRho = ((double) 200_000)/fixedParams.nRefHaps();
        this.rho = IntStream.range(0, fixedParams.nAnc())
                .mapToDouble(j -> initRho)
                .toArray();
    }

    private static double[] panelToProportion(FixedParams fixedParams) {
        IntArray nPanelHaps = fixedParams.nPanelHaps();
        int nRefHaps = fixedParams.nRefHaps();
        return IntStream.range(0, nPanelHaps.size())
                .mapToDouble(j -> (double) nPanelHaps.get(j) / nRefHaps)
                .toArray();
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

/*
 * Copyright 2021 Brian L. Browning
 * 
 * Copyright 2023 Genomics plc
 * 
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
import java.util.stream.IntStream;
import vcf.MarkerMap;

/**
 * <p>Instances of class {@code AdmixHmmProbs} stores emission probabilities
 * and components of transition probabilities.</p>
 *
 * <p>Instances of {@code AdmixHmmProbs} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixHmmProbs {

    private final ParamsInterface params;
    private final MarkerMap markerMap;

    private final double[] pRecT;             // [marker]
    private final double[] invPNoRecT;        // [marker]
    private final double[][] pNoRecTRecRho;   // [ancestry][marker]
    private final double[][] pNoRecTNoRecRho; // [ancestry][marker]
    private final double[][] qMu;             // [ancestry][ref-panel]
    private final double[][] q;               // [ancestry][ref-panel]

    /**
     * Constructs a new {@code AdmixHmmUpdater} instance from the specified
     * data.
     * @param params the analysis parameters
     * @param map the genetic position of each marker
     * @throws NullPointerException if {@code params == null || map == null}
     */
    public AdmixHmmProbs(ParamsInterface params, MarkerMap map) {
        this.params = params;
        this.markerMap = map;

        this.pRecT = AdmixUtils.pRec(params.T(), map.genDist());
        this.invPNoRecT = invPNoRecT(pRecT);
        double[][] pRecRho = pRecRho(params, map.genDist());
        this.pNoRecTRecRho = pNoRecTRecRho(pRecT, pRecRho);
        this.pNoRecTNoRecRho = pNoRecTNoRecRho(pRecT, pRecRho);

        this.qMu = ParamUtils.qMu(params);
        this.q = ParamUtils.q(params);
    }

    private static double[][] pNoRecTRecRho(double[] pRecT, double[][] pRecRho) {
        int nAnc = pRecRho.length;
        int nMarkers = pRecT.length;
        double[][] pNoRecTRecRho = new double[nAnc][];
        for (int i=0; i<nAnc; ++i) {
            double[] pRecRhoI = pRecRho[i];
            pNoRecTRecRho[i] = IntStream.range(0, nMarkers)
                    .parallel()
                    .mapToDouble(m -> (1.0 - pRecT[m])*pRecRhoI[m])
                    .toArray();
        }
        return pNoRecTRecRho;
    }

    private static double[][] pNoRecTNoRecRho(double[] pRecT, double[][] pRecRho) {
        int nAnc = pRecRho.length;
        double[][] pNoRecTNoRecRho = new double[nAnc][];
        for (int i=0; i<nAnc; ++i) {
            double[] pRecRhoI = pRecRho[i];
            pNoRecTNoRecRho[i] = IntStream.range(0, pRecT.length)
                    .parallel()
                    .mapToDouble(m -> (1.0 - pRecT[m])*(1.0 - pRecRhoI[m]))
                    .toArray();
        }
        return pNoRecTNoRecRho;
    }

    private static double[] invPNoRecT(double[] pRecT) {
        return IntStream.range(0, pRecT.length)
                .parallel()
                .mapToDouble(m -> 1.0/(1.0 - pRecT[m]))
                .toArray();
    }

    private static double[][] pRecRho(ParamsInterface params, FloatArray genDist) {
        int nAnc = params.fixedParams().nAnc();
        double[][] pRecRho = new double[nAnc][];
        for (int i=0; i<nAnc; ++i) {
            pRecRho[i] = AdmixUtils.pRec(params.rho(i), genDist);
        }
        return pRecRho;
    }


    /**
     * Returns the analysis parameters for a local ancestry inference analysis.
     * @return the analysis parameters for a local ancestry inference analysis
     */
    public ParamsInterface params() {
        return params;
    }

    /**
     * Returns the number of markers
     * @return the number of markers
     */
    public int nMarkers() {
        return markerMap.genPos().size();
    }

    /**
     * Returns the probability of a post-admixture switch to a random
     * reference haplotype in the genomic interval between the marker preceding
     * the specified marker and the specified marker.
     * @param marker marker index
     * @return the probability of a post-admixture switch to a random reference
     * haplotype
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public double pRecT(int marker) {
        return pRecT[marker];
    }

    /**
     * Returns {@code 1.0/(1.0 - this.pRecT(marker))}
     * @param marker a marker index
     * @return {@code 1.0/(1.0 - this.pRecT(marker))}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public double invPNoRecT(int marker) {
        return invPNoRecT[marker];
    }

    /**
     * Returns {@code (1.0 - this.pRecT(marker)) * this.pRecRho(marker, ancestry)}
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @return {@code (1.0 - this.pRecT(marker)) * this.pRecRho(marker, ancestry)}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().fixedParams().nAnc())}
     */
    public double pNoRecTRecRho(int marker, int ancestry) {
        return pNoRecTRecRho[ancestry][marker];
    }

    /**
     * Returns {@code (1.0 - this.pRecT(marker)) * (1.0 - this.pRecRho(marker, ancestry))}
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @return {@code (1.0 - this.pRecT(marker)) * (1.0 - this.pRecRho(marker, ancestry))}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().fixedParams().nAnc())}
     */
    public double pNoRecTNoRecRho(int marker, int ancestry) {
        return pNoRecTNoRecRho[ancestry][marker];
    }

    /**
     * Returns {@code this.pRecT(marker) * this.params().qMu(ancestry, refPanel))}
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @param refPanel a reference panel index
     * @return {@code this.pRecT(marker) * this.params().qMu(ancestry, refPanel))}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().fixedParams().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || 
     * (refPanel >= this.params().fixedParams().nRefPanels())}
     */
    public double pRecTqMu(int marker, int ancestry, int refPanel) {
        return pRecT[marker] * qMu[ancestry][refPanel];
    }

    /**
     * Returns
     * {@code this.pRecTqMu(marker, ancestry, refPanel) + 1.0 - this.pRecT(marker)}
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @param refPanel a reference panel index
     * @return {@code this.pRecTqMu(marker, ancestry, refPanel) + 1.0 - this.pRecT(marker)}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().fixedParams().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || 
     * (refPanel >= this.params().fixedParams().nRefPanels())}
     */
    public double rhoFactor(int marker, int ancestry, int refPanel) {
        return pRecTqMu(marker, ancestry, refPanel) + 1.0 - pRecT[marker];
    }

    /**
     * Returns the probability that a state transition from the marker
     * preceding the specified marker to the specified marker results
     * in a change in reference haplotype but not ancestry when the
     * destination reference haplotype is in the specified reference panel.
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @param refPanel the reference panel index
     * @return the probability that a state transition from the marker
     * preceding the specified marker to the specified marker results
     * in a change in reference haplotype but not ancestry
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().fixedParams().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || 
     * (refPanel >= this.params().fixedParams().nRefPanels())}
     */
    public double pHapChange(int marker, int ancestry, int refPanel) {
        return pRecTqMu(marker, ancestry, refPanel) + pNoRecTRecRho[ancestry][marker] * q[ancestry][refPanel];
    }

    /**
     * Returns the probability that a state transition from the marker
     * preceding the specified marker to the specified marker does not result
     * in a change in HMM state for the specified ancestry and reference
     * panel.
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @param refPanel a reference panel index
     * @return the probability that a state transition does not result
     * in a change in state
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().fixedParams().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || 
     * (refPanel >= this.params().fixedParams().nRefPanels())}
     */
    public double pNoChange(int marker, int ancestry, int refPanel) {
        return pRecTqMu(marker, ancestry, refPanel) + pNoRecTRecRho[ancestry][marker] * q[ancestry][refPanel]
                + pNoRecTNoRecRho[ancestry][marker];
    }
}

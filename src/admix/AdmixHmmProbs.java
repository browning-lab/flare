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
    private final double[][] q;
    private final double[] studyMu;
    private final double[][] sampleMu;        // [targSample][ancestry]

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

        this.q = ParamUtils.q(params);
        this.studyMu = params.studyMu();
        this.sampleMu = params.sampleData().globalAncestries(studyMu);
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
        int nAnc = params.sampleData().nAnc();
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
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargSamples() {
        return params.sampleData().targSamples().size();
    }

    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return params.sampleData().nAnc();
    }

    /**
     * Returns the study ancestry proportion for the specified ancestry
     * @param ancestry an ancestry index
     * @return the ancestry proportion for the specified sample and ancestry
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0 || ancestry >= this.nAnc())}
     */
    public double studyMu(int ancestry) {
        return studyMu[ancestry];
    }

    /**
     * Returns the ancestry proportion for the specified sample and ancestry
     * @param sample a sample index
     * @param ancestry an ancestry index
     * @return the ancestry proportion for the specified sample and ancestry
     * @throws IndexOutOfBoundsException if
     * {@code (sample < 0 || sample >= this.nTargSamples())}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0 || ancestry >= this.nAnc())}
     */
    public double sampleMu(int sample, int ancestry) {
        return sampleMu[sample][ancestry];
    }

    /**
     * Initialize the HMM forward probabilities using per-sample ancestry
     * proportions. The forward probabilities are a two-dimensional array
     * whose rows correspond to ancestry and whose columns correspond to
     * reference haplotypes.
     * @param targHap the target haplotype
     * @param fwd the arrays in which the forward probabilities will be stored
     * @param refHapToPanel a map from reference haplotype to haplotype panel
     * @param nRefHaps the number of reference haplotypes
     * @throws IndexOutOfBoundsException
     * {@code (targHap < 0 || targHap >= 2*this.nTargSamples())}
     * @throws IndexOutOfBoundsException if {@code fwd.length < this.nAnc()}
     * @throws IndexOutOfBoundsException if there is an {@code i}
     * satisfying {@code (0 <= i && i < fwd.length)} such that
     * {@code (fwd[i].length < nHapStates)}
     * @throws IndexOutOfBoundsException if
     * {@code (refHapToPanel.length < nHapStates)}
     * @throws NullPointerException if {@code (fwd == null)} or if there
     * is an {@code i} satisfying {@code (0 <= i && i < fwd.length)} such
     * that {@code (fwd[i] == null)}
     * @throws NullPointerException if {@code (refHapToPanel == null)}
     */
    public void initFwd(int targHap, double[][] fwd, short[] refHapToPanel,
            int nRefHaps) {
        int nAnc = params.sampleData().nAnc();
        int sample = targHap >> 1;
        for (int i=0; i<nAnc; ++i) {
            for (int h=0; h<nRefHaps; ++h) {
                fwd[i][h] = sampleMu[sample][i]*q[i][refHapToPanel[h]];
            }
        }
    }

   /**
     * Initialize the HMM forward probabilities using study ancestry
     * proportions. The forward probabilities are a two-dimensional array
     * whose rows correspond to ancestry and whose columns correspond to
     * reference haplotypes.
     * @param fwd the arrays in which the forward probabilities will be stored
     * @param refHapToPanel a map from reference haplotype to haplotype panel
     * @param nRefHaps the number of reference haplotypes
     * @throws IndexOutOfBoundsException if {@code fwd.length < this.nAnc()}
     * @throws IndexOutOfBoundsException if there is an {@code i}
     * satisfying {@code (0 <= i && i < fwd.length)} such that
     * {@code (fwd[i].length < nHapStates)}
     * @throws IndexOutOfBoundsException if
     * {@code (refHapToPanel.length < nHapStates)}
     * @throws NullPointerException if {@code (fwd == null)} or if there
     * is an {@code i} satisfying {@code (0 <= i && i < fwd.length)} such
     * that {@code (fwd[i] == null)}
     * @throws NullPointerException if {@code (refHapToPanel == null)}
     */
    public void initFwd(double[][] fwd, short[] refHapToPanel, int nRefHaps) {
        int nAnc = params.sampleData().nAnc();
        for (int i=0; i<nAnc; ++i) {
            for (int h=0; h<nRefHaps; ++h) {
                fwd[i][h] = studyMu[i]*q[i][refHapToPanel[h]];
            }
        }
    }

    /**
     * Set shift parameter for each reference panel for the HMM forward update
     * using per-sample ancestry proportions.
     * @param targHap the target haplotype
     * @param marker the marker
     * @param ancestry the ancestry
     * @param pAnc the state probability sum for each ancestry at the previous
     * marker
     * @param shift the shift parameter for each reference panel in the
     * HMM forward update
     * @throws IndexOutOfBoundsException if
     * {@code (targHap < 0) || (targHap >=  2*this.params().sampleData().targSamples().size())}
     * @throws IndexOutOfBoundsException if
     * {@code (marker < 0) || (marker >= this.nMarkers())}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().sampleData().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || (refPanel >= this.params().sampleData().nRefPanels())}
     * @throws IndexOutOfBoundsException if {@code refPanel >= shift.length}
     * @throws NullPointerException if {@code shift == null}
     */
    public void setFwdShift(int targHap, int marker, int ancestry,
            double pAnc, double[] shift) {
        int sample = targHap >> 1;
        for (int j=0; j<shift.length; ++j) {
            shift[j] = pRecT[marker]*sampleMu[sample][ancestry]*q[ancestry][j]
                    + pNoRecTRecRho[ancestry][marker]*q[ancestry][j]*pAnc;
        }
    }

    /**
     * Set shift parameter for each reference panel for the HMM forward update
     * using study ancestry proportions.
     * @param marker the marker
     * @param ancestry the ancestry
     * @param pAnc the normalized state probability sum for each ancestry
     * at the previous marker
     * @param shift the shift parameter when updating HMM state probabilities
     * @throws IndexOutOfBoundsException if
     * {@code (marker < 0) || (marker >= this.nMarkers())}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().sampleData().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || (refPanel >= this.params().sampleData().nRefPanels())}
     * @throws IndexOutOfBoundsException if {@code refPanel >= shift.length}
     * @throws NullPointerException if {@code shift == null}
     */
    public void setFwdShift(int marker, int ancestry, double pAnc, double[] shift) {
        for (int j=0; j<shift.length; ++j) {
            shift[j] = pRecT[marker]*studyMu[ancestry]*q[ancestry][j]
                    + pNoRecTRecRho[ancestry][marker]*q[ancestry][j]*pAnc;
        }
    }

    /**
     * Returns the shift parameter for the HMM backward update.
     * @param marker the marker
     * @param ancestry the ancestry
     * @param bwdAncSum an ancestry-dependent sum
     * @param bwdSum the sum of {@code bwdAncSum} over all ancestries
     * @return the shift parameter for the HMM backward update
     * @throws IndexOutOfBoundsException if
     * {@code (marker < 0) || (marker >= this.nMarkers())}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().sampleData().nAnc())}
     */
    public double bwdShift(int marker, int ancestry, double bwdAncSum, double bwdSum) {
        return pRecT[marker]*bwdSum + pNoRecTRecRho[ancestry][marker]*bwdAncSum;
    }

    /**
     * Returns the inverse of the probability of no recombination between
     * the preceding marker and the current marker since admixture.
     * @param marker a marker index
     * @return the inverse of the probability of no recombination between
     * the preceding marker and the current marker since admixture
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public double invPNoRecT(int marker) {
        return invPNoRecT[marker];
    }

    /**
     * Returns the probability of no recombination between the preceding marker
     * and the current marker
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @return the probability of no recombination between the preceding marker
     * and the current marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().sampleData().nAnc())}
     */
    public double pNoRecTNoRecRho(int marker, int ancestry) {
        return pNoRecTNoRecRho[ancestry][marker];
    }

    /**
     * Returns a factor used to estimate the intensity parameter for the
     * exponential distribution of the probability of recombination
     * prior to admixture.
     * @param marker a marker index
     * @param ancestry an ancestry index
     * @param refPanel a reference panel index
     * @return a factor used to estimate the intensity parameter for the
     * exponential distribution of the probability of recombination
     * prior to admixture
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (ancestry < 0) || (ancestry >= this.params().sampleData().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) ||
     * (refPanel >= this.params().sampleData().nRefPanels())}
     */
    public double rhoFactor(int marker, int ancestry, int refPanel) {
        return  pRecT[marker]*studyMu[ancestry]*q[ancestry][refPanel] + 1.0 - pRecT[marker];
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
     * {@code (ancestry < 0) || (ancestry >= this.params().sampleData().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || 
     * (refPanel >= this.params().sampleData().nRefPanels())}
     */
    public double pHapChange(int marker, int ancestry, int refPanel) {
        return pRecT[marker]*studyMu[ancestry]*q[ancestry][refPanel]
                + pNoRecTRecRho[ancestry][marker]*q[ancestry][refPanel];
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
     * {@code (ancestry < 0) || (ancestry >= this.params().sampleData().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) ||
     * (refPanel >= this.params().sampleData().nRefPanels())}
     */
    public double pNoChange(int marker, int ancestry, int refPanel) {
        return pRecT[marker]*studyMu[ancestry]*q[ancestry][refPanel]
                + pNoRecTRecRho[ancestry][marker]*q[ancestry][refPanel]
                + pNoRecTNoRecRho[ancestry][marker];
    }
}

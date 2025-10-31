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

import blbutil.Const;
import blbutil.DoubleArray;
import blbutil.FloatArray;
import ints.IntArray;
import ints.IntList;
import java.io.PrintWriter;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.stream.IntStream;
import vcf.Marker;
import vcf.RefGT;

/**
 * <p>Class {@code EstimatedPanels} stores estimated reference panel
 * probabilities and mean recombination intensities for the target
 * haplotypes in a set of marker windows.</p>
 *
 * <p>Instances of class {@code EstimatedPanels} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class EstimatedPanels {

    private final SampleData sampleData;
    private final IntArray selectedTargHaps;
    private final RefGT targRefGT;
    private final int nRefPanels;
    private final int refPanelProbsLength;
    private final int nSelectedTargHaps;
    private final int[] windowEnds;
    private final double[] windowToProbsScaleFactor;
    private final double[] windowToNeScaleFactor;

    private final AdmixRecBuilder.ProbFormatter probFormatter;
    private final AtomicReferenceArray<FloatArray> hapToWindowPanelProbs;
    private final AtomicReferenceArray<FloatArray> hapToWindowScaledNe;

    /**
     * Constructs a new {@code EstimatedPanels} instance from the specified
     * data.
     * @param chromData immutable input data for local ancestry inference
     * on a chromosome
     * @throws NullPointerException if {@code (chromData == null)}
     */
    public EstimatedPanels(AdmixChromData chromData) {
        this.sampleData = chromData.sampleData();
        this.selectedTargHaps = selectedTargHaps(chromData);
        this.nSelectedTargHaps = selectedTargHaps.size();
        this.targRefGT = chromData.targRefGT();
        this.nRefPanels = sampleData.nRefPanels();
        this.refPanelProbsLength = targRefGT.nMarkers()*nRefPanels;
        DoubleArray genMap = chromData.map().genPos();
        this.windowEnds = windowEnds(genMap, chromData.par().panel_cm());
        this.windowToProbsScaleFactor = IntStream.range(0, windowEnds.length)
                .parallel()
                .mapToDouble(w -> 1.0/(windowEnds[w] - (w==0 ? 0 : windowEnds[w-1])))
                .toArray();
        this.windowToNeScaleFactor = IntStream.range(0, windowEnds.length)
                .parallel()
                .mapToDouble(w -> 100.0/
                        (genMap.get(windowEnds[w]-1) - genMap.get(w==0 ? 0 : windowEnds[w-1])))
                .toArray();
        this.probFormatter = new AdmixRecBuilder.ProbFormatter(3);
        this.hapToWindowPanelProbs = new AtomicReferenceArray<>(nSelectedTargHaps);
        this.hapToWindowScaledNe = new AtomicReferenceArray<>(nSelectedTargHaps);
    }

    private static IntArray selectedTargHaps(AdmixChromData chromData) {
        int nTargHaps = chromData.nTargHaps();
        AdmixPar par = chromData.par();
        int maxHaps = par.panel_haps();
        long seed = par.seed() + 103;
        IntList hapList = new IntList(nTargHaps);
        for (int h=0; h<nTargHaps; ++h) {
            hapList.add(h);
        }
        return ObservedHaps.randomSubset(hapList, maxHaps, seed);
    }

    private static int[] windowEnds(DoubleArray genPos, double windowCm) {
        IntList ends = new IntList(0);
        double nextEndGenPos = genPos.get(0) + windowCm;
        for (int m=0, n=genPos.size(); m<n; ++m) {
            if (genPos.get(m) > nextEndGenPos) {
                ends.add(m);
                nextEndGenPos += windowCm;
            }
        }
        ends.add(genPos.size());
        return ends.toArray();
    }

    /**
     * Returns the specified target haplotype with reference panel
     * probability and recombination intensity data
     * @param index an index in the list of target haplotypes with reference
     * panel and recombination intensity data
     * @return the specified target haplotype
     * @throws IndexOutOfBoundsException if
     * {@code ((index < 0) || (index >= this.nSelectedTargHaps()))}
     */
    public int selectedTargHap(int index) {
        return selectedTargHaps.get(index);
    }

    /**
     * Returns the number of target haplotypes with reference panel
     * probability and recombination intensity data.
     * @return the number of target haplotypes with reference panel
     * probability and recombination intensity data
     */
    public int nSelectedTargHaps() {
        return selectedTargHaps.size();
    }

    /**
     * Returns the number of reference panels.
     * @return the number of reference panels
     */
    public int nRefPanels() {
        return nRefPanels;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return targRefGT.nMarkers();
    }

    /**
     * Stores the specified reference panel probabilities and recombination
     * intensities for the specified target haplotype. The estimated
     * probability of reference panel {@code j} at marker {@code m} is required
     * to be stored in element {@code (m*this.nRefPanels() + j)} of the
     * {@code refPanelProbs} parameter.
     * @param index an index in the list of target haplotypes with
     * reference panel probability and recombination intensity data
     * @param refPanelProbs the reference panel probabilities
     * @param recombIntensities an array whose {@code k}-th element
     * is the recombination intensity between the {@code k}-th marker and
     * the preceding marker
     * @throws IllegalArgumentException if
     * {@code refPanelProbs.length != (this.nRefPanels() * this.nMarkers())}
     * @throws IllegalArgumentException if
     * {@code recombIntensities.length != this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || targHap >= this.nSelectedTargHaps()}
     * @throws NullPointerException if
     * {@code ((refPanelProbs == null) || (recombIntensities == null))}
     */
    public void set(int index, double[] refPanelProbs, double[] recombIntensities) {
        if (refPanelProbs.length != refPanelProbsLength) {
            throw new IllegalArgumentException(String.valueOf(refPanelProbs.length));
        }
        if (recombIntensities.length != targRefGT.nMarkers()) {
            throw new IllegalArgumentException(String.valueOf(recombIntensities.length));
        }
        hapToWindowPanelProbs.getAndSet(index, meanWindowPanelProbs(refPanelProbs));
        hapToWindowScaledNe.getAndSet(index, windowScaledNe(recombIntensities));
    }

    private FloatArray meanWindowPanelProbs(double[] panelProbs) {
        double[] meanProbs = new double[windowEnds.length*nRefPanels];
        int mStart = 0;
        int jStart = 0;
        int index = 0;
        for (int w=0; w<windowEnds.length; ++w) {
            int mEnd = windowEnds[w];
            int jEnd = jStart + nRefPanels;
            for (int m=mStart; m<mEnd; ++m) {
                for (int j=jStart; j<jEnd; ++j) {
                    meanProbs[j] += panelProbs[index++];
                }
            }
            for (int j=jStart; j<jEnd; ++j) {
                meanProbs[j] *= windowToProbsScaleFactor[w];
            }
            jStart = jEnd;
            mStart = mEnd;
        }
        return new FloatArray(meanProbs);
    }

    private FloatArray windowScaledNe(double[] recombIntensities) {
        double[] windowScaledNe = new double[windowEnds.length];
        int mStart = 0;
        for (int w=0; w<windowEnds.length; ++w) {
            int mEnd = windowEnds[w];
            for (int m=mStart; m<mEnd; ++m) {
                windowScaledNe[w] += recombIntensities[m];
            }
            windowScaledNe[w] *= windowToNeScaleFactor[w];
            mStart = mEnd;
        }
        return new FloatArray(windowScaledNe);
    }

    /**
     * Writes the reference panel probabilities for each marker window and
     * reference haplotype to the specified {@code PrintWriter}.
     * @param out the {@code PrintWriter} to which output will be printed
     * @throws NullPointerException if {@code (out == null)}
     */
    public void writePanelProbabilities(PrintWriter out) {
        int minWindowMarkers = sampleData.par().panel_markers();
        String[] output = IntStream.range(0, windowEnds.length)
                .parallel()
                .filter(w -> (windowEnds[w] - (w==0 ? 0 : windowEnds[w-1])) >= minWindowMarkers)
                .mapToObj(w -> windowReport(w))
                .toArray(String[]::new);
        for (String s : output) {
            out.print(s);
        }
    }

    private String windowReport(int window) {
        StringBuilder sb = new StringBuilder();
        int start = window*nRefPanels;
        int end = start + nRefPanels;
        String prefix = panelProbPrefix(window);
        for (int j=0; j<nSelectedTargHaps; ++j) {
            sb.append(prefix);
            sb.append(selectedTargHaps.get(j));
            sb.append(Const.tab);
            sb.append(hapToWindowScaledNe.get(j).get(window));
            FloatArray windowPanelProbs = hapToWindowPanelProbs.get(j);
            for (int k=start; k<end; ++k) {
                sb.append(Const.tab);
                sb.append(probFormatter.format(windowPanelProbs.get(k)));
            }
            sb.append(Const.nl);
        }
        return sb.toString();
    }

    private String panelProbPrefix(int window) {
        int start = window==0 ? 0 : windowEnds[window-1];
        int end = windowEnds[window];
        Marker startMarker = targRefGT.marker(start);
        Marker endMarker = targRefGT.marker(end-1);
        StringBuilder prefix = new StringBuilder(30);
        prefix.append(window);
        prefix.append(Const.tab);
        prefix.append(startMarker.chromID());
        prefix.append(Const.colon);
        prefix.append(startMarker.pos());
        prefix.append(Const.hyphen);
        prefix.append(endMarker.pos());
        prefix.append(Const.tab);
        return prefix.toString();
    }
}

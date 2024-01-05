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
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.stream.IntStream;

/**
 * <p>Class {@code GlobalAncestryProbs} stores the estimated global ancestry
 * proportion for each target sample.</p>
 *
 * <p>Instances of class {@code EstimatedGlobalAncProportions} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class EstimatedGlobalAncProportions {

    private final String[] sampleIds;
    private final String[] ancIds;
    private final AtomicReferenceArray<DoubleArray> sampAncProbs;

    /**
     * Constructs a new {@code GlobalAncProbs} instance from the specified
     * data.
     * @param sampleIds the sample identifiers
     * @param ancIds the list of ancestry identifiers
     * @throws NullPointerException if {@code (sampleIds == null) || (ancIds == null)}
     */
    public EstimatedGlobalAncProportions(String[] sampleIds, String[] ancIds) {
        this.sampleIds = sampleIds.clone();
        this.ancIds = ancIds.clone();
        this.sampAncProbs = initialSampAncProbs(sampleIds.length, ancIds.length);
    }

    private static AtomicReferenceArray<DoubleArray> initialSampAncProbs(
            int nSamples, int nAnc) {
        DoubleArray ancProbs = new DoubleArray(new double[nAnc]);
        DoubleArray[] fa2 = IntStream.range(0, nSamples)
                .parallel()
                .mapToObj(j -> ancProbs)
                .toArray(DoubleArray[]::new);
        return new AtomicReferenceArray<>(fa2);
    }

    /**
     * Returns the list of sample identifiers.
     * @return the list of sample identifiers
     */
    public String[] sampleIds() {
        return sampleIds.clone();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return sampleIds.length;
    }

    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return ancIds.length;
    }

    /**
     * Stores the specified ancestry probabilities for the specified
     * haplotype.
     * @param hap a haplotype index
     * @param nMarkers the number of markers
     * @param ancProbs the ancestry probabilities
     * @throws IllegalArgumentException if
     * {@code ancProbs.length != (this.nAnc() * nMarkers)}
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= 2*this.nSamples()}
     * @throws NullPointerException if {@code ancProbs == null}
     */
    public void add(int hap, int nMarkers, double[] ancProbs) {
        if (nMarkers*ancIds.length != ancProbs.length) {
            throw new IllegalArgumentException(String.valueOf(nMarkers));
        }
        int sample = hap >> 1;
        sampAncProbs.getAndUpdate(sample, probs -> update(probs, nMarkers, ancProbs));
    }

    private DoubleArray update(DoubleArray sampAncProbs, int nMarkers,
            double[] ancProbs) {
        double[] da = sampAncProbs.toArray();
        int index = 0;
        for (int m=0; m<nMarkers; ++m) {
            for (int a=0; a<ancIds.length; ++a) {
                da[a] += ancProbs[index++];
            }
        }
        return new DoubleArray(da);
    }

    /**
     * Prints one tab-delimited line per sample to the specified
     * {@code PrintWriter}. The first field is the sample identifier.
     * The remaining {@code this.nAnc()} fields are the posterior ancestry
     * probabilities for the sample. The posterior probability of an ancestry
     * for an individual is the mean ancestry probability calculated from
     * all markers and both haplotypes.
     * @param out the {@code PrintWriter} to which global ancestry proportions
     * will be written
     * @throws NullPointerException if {@code out == null}
     */
    public void writeGlobalAncestry(PrintWriter out) {
        String[] df3 = df3();
        int step = 1<<16;
        printHeaderLine(out);
        for (int start=0; start<sampleIds.length; start+=step) {
            int end = Math.min(start + step, sampleIds.length);
            printGlobalProbs(start, end, df3, out);
        }
    }

    private void printHeaderLine(PrintWriter out) {
        out.print("SAMPLE");
        for (String anc : ancIds) {
            out.print(Const.tab);
            out.print(anc);
        }
        out.println();
    }

    private static String[] df3() {
        DecimalFormat df3 = new DecimalFormat("#.###");
        return IntStream.range(0, 1001)
                .parallel()
                .mapToObj(j -> df3.format(0.001*j))
                .toArray(String[]::new);
    }

    private void printGlobalProbs(int start, int end,
            String[] df3, PrintWriter out) {
        String[] sa = IntStream.range(start, end)
                .parallel()
                .mapToObj(s -> globalProbsString(s, df3))
                .toArray(String[]::new);
        for (String s : sa) {
            out.println(s);
        }
    }

    private String globalProbsString(int sample, String[] df3) {
        DoubleArray fa = sampAncProbs.get(sample);
        StringBuilder sb = new StringBuilder(20 + 5*ancIds.length);
        sb.append(sampleIds[sample]);
        float sum = 0f;
        for (int j=0; j<ancIds.length; ++j) {
            sum += fa.get(j);
        }
        float factor = 1000f/sum;
        for (int j=0; j<ancIds.length; ++j) {
            sb.append(Const.tab);
            sb.append(df3[ (int) Math.rint(factor*fa.get(j)) ]);
        }
        return sb.toString();
    }
}

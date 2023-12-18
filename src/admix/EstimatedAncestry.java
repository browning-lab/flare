/*
 * Copyright 2021 Brian L. Browning
 * 
 * Copyright 2023 Genomics plc
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
import ints.IntArray;
import ints.IntList;
import ints.UnsignedByteArray;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.stream.IntStream;
import vcf.RefGT;

/**
 * <p>Class {@code EstimatedAncestry} stores estimated local ancestry of
 * target haplotypes.</p>
 *
 * <p>Instances of class {@code EstimatedAncestry} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class EstimatedAncestry {

    private final FixedParams fixedParams;
    private final RefGT targRefGT;
    private final int nAnc;
    private final int ancProbsLength;
    private final int nTargHaps;

    private final boolean storeAncProbs;
    private final AdmixRecBuilder.ProbFormatter probFormatter;
    private final AtomicReferenceArray<FloatArray> hapToAncProbs;
    private final AtomicReferenceArray<IntArray> hapToAnc;
    private final GlobalAncProbs globalAncProbs;

    /**
     * Constructs a new {@code EstimatedAncestry} instance from the specified
     * data.
     * @param admixData immutable input data for a chromosome and analysis
     * parameters for local ancestry inference
     * @param globalAncProbs an object for storing global ancestry
     * probabilities
     * @throws NullPointerException if
     * {@code (admixData == null) || (globalAncProbs == null)}
     */
    public EstimatedAncestry(AdmixChromData chromData, FixedParams fixedParams, GlobalAncProbs globalAncProbs) {
        this.fixedParams = fixedParams;
        AdmixPar par = fixedParams.par();
        this.targRefGT = chromData.targRefGT();
        this.nAnc = fixedParams.nAnc();
        this.ancProbsLength = targRefGT.nMarkers()*nAnc;
        this.nTargHaps = chromData.nTargHaps();

        this.storeAncProbs = par.probs();
        if (storeAncProbs) {
            probFormatter = new AdmixRecBuilder.ProbFormatter();
            hapToAncProbs = new AtomicReferenceArray<>(nTargHaps);
            hapToAnc = null;
        }
        else {
            probFormatter = null;
            hapToAncProbs = null;
            hapToAnc = new AtomicReferenceArray<>(nTargHaps);
        }
        this.globalAncProbs = globalAncProbs;
    }

    /**
     * Returns the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return nTargHaps;
    }

    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return nAnc;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return targRefGT.nMarkers();
    }

    /**
     * Stores the specified ancestry probabilities.
     * If the number of ancestries is {@code nAnc}, then the estimated
     * probability of ancestry at marker {@code m} and ancestry {@code i}
     * is stored in element {@code m*nAnc + i} of the returned list.
     * @param targHap a target haplotype index
     * @param ancProbs the ancestry probabilities
     * @throws IllegalArgumentException if
     * {@code ancProbs.length != (this.nAnc() * this.nMarkers())}
     * @throws IndexOutOfBoundsException if
     * {@code targHap < 0 || targHap >= this.nTargHaps()}
     * @throws NullPointerException if {@code ancProbs == null}
     *
     */
    public void set(int targHap, double[] ancProbs) {
        if (ancProbs.length != ancProbsLength) {
            throw new IllegalArgumentException(String.valueOf(ancProbs.length));
        }
        if (storeAncProbs) {
            hapToAncProbs.getAndUpdate(targHap, x -> new FloatArray(ancProbs));
        }
        else {
            hapToAnc.getAndUpdate(targHap, x -> mostProbableAncestry(ancProbs));
        }
        globalAncProbs.add(targHap, targRefGT.nMarkers(), ancProbs);
    }

    private IntArray mostProbableAncestry(double[] ancProbs) {
        int[] ia = new int[targRefGT.nMarkers()];
        for (int m=0, offset=0; offset<ancProbs.length; ++m, offset+=nAnc) {
            int nextOffset = offset + nAnc;
            int maxIndex = offset;
            for (int k=(offset+1); k<nextOffset; ++k) {
                if (ancProbs[k]>ancProbs[maxIndex]) {
                    maxIndex = k;
                }
            }
            ia[m] = maxIndex-offset;
        }
        return IntArray.packedCreate(ia, nAnc);
    }

    /**
     * Returns a list of byte arrays with GZIP compressed VCF records
     * for the specified markers that containing phased genotypes and
     * inferred local ancestry.
     * @param start the first marker (inclusive)
     * @param end the last marker (exclusive)
     * @return a list of byte arrays with GZIP compressed VCF records
     * @throws IllegalArgumentException if
     * {@code start < 0 || end >= this.nMarkers()}
     * @throws IllegalArgumentException if {@code end < start}
     */
    public UnsignedByteArray[] writeAncestry(int start, int end) {
        if (start<0) {
            throw new IndexOutOfBoundsException(String.valueOf(start));
        }
        if (end<start || end>targRefGT.nMarkers()) {
            throw new IndexOutOfBoundsException(String.valueOf(end));
        }
        int nThreads = fixedParams.par().nthreads();
        int stepSize = (end - start + nThreads - 1)/nThreads;
        IntList partEnds = new IntList(1<<8);
        for (int j=start; j<end; j+=stepSize) {
            partEnds.add(j);
        }
        partEnds.add(end);
        return  IntStream.range(1, partEnds.size())
                    .parallel()
                    .mapToObj(j -> partCompressedOutput(partEnds.get(j-1), partEnds.get(j)))
                    .toArray(UnsignedByteArray[]::new);
    }

    private UnsignedByteArray partCompressedOutput(int start, int end) {
        AdmixRecBuilder recBuilder = new AdmixRecBuilder(fixedParams, targRefGT, start, end);
        int nTargSamples = nTargHaps>>1;
        if (storeAncProbs) {
            for (int s=0; s<nTargSamples; ++s) {
                int h1 = s<<1;
                int h2 = h1 | 0b1;
                recBuilder.addSampleData(hapToAncProbs.get(h1), hapToAncProbs.get(h2), probFormatter);
            }
        }
        else {
            for (int s=0; s<nTargSamples; ++s) {
                int h1 = s<<1;
                int h2 = h1 | 0b1;
                recBuilder.addSampleData(hapToAnc.get(h1), hapToAnc.get(h2));
            }
        }
        return recBuilder.toUnsignedByteArray();
    }
}

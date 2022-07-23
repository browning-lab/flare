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
package vcf;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code LowMafGTRc} stores genotypes for a list of samples
 * at a marker. Instances of class {@code LowMafGTRec} store lists of
 * haplotypes carrying each non-major or missing allele. All genotypes are
 * considered to be unphased if any sample has an unphased or missing
 * genotype.</p>
 *
 * <p>Instances of class {@code LowMafGTRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LowMafGTRec implements GTRec {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int majorAllele;
    private final int[][] hapIndices;
    private final int[] missingSamples;
    private final boolean isPhased;

    @Override
    public long estBytes() {
        int overhead = (2 + (hapIndices.length - 1))*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + (4 + hapIndices.length)*8;  // assume 8 bytes per reference
        estBytes += (1 + hapIndices.length)*4; // 4 bytes per array to store length;
        for (int j=0; j<hapIndices.length; ++j) {
            if (j!=majorAllele) {
                estBytes += 4*hapIndices[j].length;
            }
            else {
                assert hapIndices[j]==null;
            }
        }
        estBytes += 4*missingSamples.length;
        return estBytes;
    }

    /**
     * Constructs a new {@code LowMafGTRec} representing
     * the specified VCF record's GT format field data.
     * @param listRep the VCF record genotype data
     * @throws NullPointerException if {@code listRep == null}
     */
    public LowMafGTRec(VcfRecGTParser.HapListRep listRep) {
        this.marker = listRep.marker();
        this.samples = listRep.samples();
        this.nHaps = samples.size()<<1;
        this.majorAllele = listRep.majorAllele();
        this.hapIndices = listRep.hapLists(true);
        this.missingSamples = listRep.missingSamples();
        this.isPhased = listRep.isPhased();
    }

    @Override
    public boolean isPhased(int sample) {
        if (sample < 0 || sample >= this.samples().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        return isPhased;
    }

    @Override
    public boolean isPhased() {
        return isPhased;
    }

    @Override
    public Samples samples() {
        return samples;
    }


    @Override
    public int size() {
        return nHaps;
    }

    @Override
    public Marker marker() {
        return marker;
    }


    @Override
    public int allele1(int sample) {
        return get(sample<<1);
    }

    @Override
    public int allele2(int sample) {
        return get((sample<<1) | 0b1);
    }

    @Override
    public int get(int hap) {
        if (hap < 0 || hap >= nHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        for (int j=0; j<hapIndices.length; ++j) {
            if (j != majorAllele) {
                if (Arrays.binarySearch(hapIndices[j], hap) >= 0) {
                    return j;
                }
            }
        }
        return Arrays.binarySearch(missingSamples, (hap>>1))>=0 ? -1 : majorAllele;
    }

    @Override
    public int[] alleles() {
        int[] ia = IntStream.range(0, nHaps)
                .map(h -> majorAllele)
                .toArray();
        for (int al=0; al<hapIndices.length; ++al) {
            if (al != majorAllele) {
                for (int h : hapIndices[al]) {
                    ia[h] = al;
                }
            }
        }
        for (int s : missingSamples) {
            int h1 = s<<1;
            int h2 = h1 | 0b1;
            ia[h1] = -1;
            ia[h2] = -1;
        }
        return ia;
    }

    /**
     * Returns the major allele.  If there are two major alleles with the same
     * allele count, the smaller allele is returned.
     * @return the major allele
     */
    public int majorAllele() {
        return majorAllele;
    }

    /**
     * Returns the number of copies of the specified allele.
     * @param allele an allele
     * @return the number of copies of the specified allele
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker().nAlleles()}
     */
    public int alleleCount(int allele) {
        if (allele==majorAllele) {
            int n = nHaps - (missingSamples.length<<1);
            for (int al=0; al<hapIndices.length; ++al) {
                if (al!=majorAllele) {
                    n -= hapIndices.length;
                }
            }
            return n;
        }
        else {
            return hapIndices[allele].length;
        }
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field. The returned VCF record
     * will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field
     */
    @Override
    public String toString() {
        return GTRec.toVcfRec(this);
    }
}

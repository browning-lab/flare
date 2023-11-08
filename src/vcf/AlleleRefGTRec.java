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

import bref4.BrefRec;
import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code AlleleRefGTRec} represent represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.</p>
 *
 * <p>Class {@code AlleleRefGTRec} stores the haplotypes that carry each
 * non-major allele.</p>
 *
 * <p>Instances of class {@code AlleleRefGTRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AlleleRefGTRec implements RefGTRec {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int majorAllele;
    private final int[][] alleleToHaps;

    @Override
    public long estBytes() {
        int overhead = (1 + (alleleToHaps.length - 1))*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + (3 + alleleToHaps.length)*8;  // assume 8 bytes per reference
        estBytes += (1 + alleleToHaps.length)*4; // 4 bytes per array to store length;
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (j!=majorAllele) {
                estBytes += 4*(1 + alleleToHaps[j].length); // 4 bytes per array to store length;
            }
            else {
                assert alleleToHaps[j]==null;
            }
        }
        estBytes += 8;  // to store nHaps and majorAllele
        return estBytes;
    }

    /**
     * Constructs a new {@code AlleleRefGTRec} instance from the specified data.
     *
     * @param rec the phased, non-missing genotype data
     * @throws NullPointerException if {@code rec == null}
     */
    public AlleleRefGTRec(RefGTRec rec) {
        this.alleleToHaps = rec.alleleToHaps();
        int majAllele = 0;
        while (alleleToHaps[majAllele]!=null) {
            ++majAllele;
        }
        this.marker = rec.marker();
        this.samples = rec.samples();
        this.nHaps = rec.size();
        this.majorAllele = majAllele;
    }

    /**
     * Constructs a new {@code AlleleRefGTRc} instance from the specified
     * data.
     *
     * @param gtp a VCF record parser that extracts sample genotypes
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased genotype or missing allele
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws NullPointerException if {@code gtp == null}
     */
    public AlleleRefGTRec(VcfRecGTParser gtp) {
        this.marker = gtp.marker();
        this.samples = gtp.samples();
        this.nHaps = 2*gtp.nSamples();
        this.alleleToHaps = gtp.nonMajRefIndices();
        int majAl = -1;
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (alleleToHaps[j]==null) {
                majAl = j;
                break;
            }
        }
        this.majorAllele = majAl;
    }

    /**
     * Constructs a new {@code AlleleRefGTRec} instance from the specified data.
     * The specified {@code hapIndices} array is required to contain exactly one
     * {@code null} element. The {@code null} element should be the major
     * allele because this is most memory-efficient, but this requirement is not
     * enforced. A haplotype index should be an element of only one array in
     * {@code hapIndices}. If a haplotype index is an element of more than
     * one array in {@code hapIndices}, it is assigned to the array with
     * smallest index.

     * @param marker the marker
     * @param samples the samples
     * @param hapIndices whose {@code j}-th element is a list of haplotypes
     * sorted in increasing order that carry the {@code j}-th allele, or is
     * {@code null}
     *
     * @throws IllegalArgumentException if the {@code (hapIndices} does
     * not contain exactly one {@code null} element
     * @throws IllegalArgumentException if any non-null element of
     * {@code hapIndices} is not a sorted list of distinct haplotype indices
     * between 0 (inclusive) and {@code 2*samples.size()} (exclusive)
     * @throws IllegalArgumentException if
     * {@code marker.nAlleles() != hapIndices.length}
     * @throws NullPointerException if
     * {@code marker == null || samples == null || hapIndices == null}
     */
    public AlleleRefGTRec(Marker marker, Samples samples, int[][] hapIndices) {
        this.marker = marker;
        this.samples = samples;
        this.nHaps = 2*samples.size();
        this.majorAllele = checkIndicesAndReturnNullIndex(hapIndices, nHaps);
        this.alleleToHaps = deepCopy(hapIndices);
    }

    static int checkIndicesAndReturnNullIndex(int[][] hapIndices, int nHaps) {
        int majAllele = -1;
        for (int j=0; j<hapIndices.length; ++j) {
            if (hapIndices[j]==null) {
                if (majAllele == -1) {
                    majAllele = j;
                }
                else {
                    throwArrayError();
                }
            }
            else {
                checkSorted(hapIndices[j], nHaps);
            }
        }
        if (majAllele == -1) {
            throwArrayError();
        }
        return majAllele;
    }

    private static void throwArrayError() {
        throw new IllegalArgumentException("invalid array");
    }

    private static void checkSorted(int[] ia, int nHaps) {
        if (ia.length>0 && (ia[0] < 0 || ia[ia.length - 1] >= nHaps)) {
            throwArrayError();
        }
        for (int k=1; k<ia.length; ++k) {
            if (ia[k-1] >= ia[k]) {
                throwArrayError();
            }
        }
    }

    static int[][] deepCopy(int[][] ia) {
        int[][] copy = new int[ia.length][];
        for (int j=0; j<ia.length; ++j) {
            if (ia[j]!=null) {
                copy[j] = ia[j].clone();
            }
        }
        return ia;
    }

    @Override
    public int[][] alleleToHaps() {
        return deepCopy(alleleToHaps);
    }

    @Override
    public IndexArray hapToAllele() {
        return new IndexArray(toIntArray(), alleleToHaps.length);
    }

    private IntArray toIntArray() {
        int[] ia = IntStream.range(0, nHaps)
                .map(h -> majorAllele)
                .toArray();
        for (int al=0; al<alleleToHaps.length; ++al) {
            if (alleleToHaps[al]!=null) {
                for (int h : alleleToHaps[al]) {
                    ia[h] = al;
                }
            }
        }
        return IntArray.packedCreate(ia, alleleToHaps.length);
    }

    @Override
    public int nAlleleCodedHaps() {
        return BrefRec.nonNullCnt(alleleToHaps);
    }

    @Override
    public boolean isPhased(int sample) {
        if (sample < 0 || sample >= this.samples().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        return true;
    }

    /**
     * Returns {@code true}.
     * @return {@code true}
     */
    @Override
    public boolean isPhased() {
        return true;
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
    public int get(int hap) {
        if (hap < 0 || hap >= nHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (j != majorAllele) {
                if (Arrays.binarySearch(alleleToHaps[j], hap) >= 0) {
                    return j;
                }
            }
        }
        return majorAllele;
    }

    @Override
    public boolean isAlleleCoded() {
        return true;
    }

    @Override
    public int majorAllele() {
        return majorAllele;
    }

    @Override
    public int[] alleleCounts() {
        int[] alCnts = new int[marker.nAlleles()];
        alCnts[majorAllele] = nHaps;
        for (int j=0; j<alCnts.length; ++j) {
            if (j!=majorAllele) {
                int alCnt = alleleToHaps[j].length;
                alCnts[j] = alCnt;
                alCnts[majorAllele] -= alCnt;
            }
        }
        return alCnts;
    }

    @Override
    public int alleleCount(int allele) {
        if (alleleToHaps[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return alleleToHaps[allele].length;
        }
    }

    @Override
    public int hapIndex(int allele, int copy) {
        if (alleleToHaps[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return alleleToHaps[allele][copy];
        }
    }

    @Override
    public boolean isCarrier(int allele, int hap) {
        return get(hap)==allele;
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

    @Override
    public int nMaps() {
        return 1;
    }

    @Override
    public IntArray[] maps() {
        return new IntArray[] {toIntArray()};
    }

    @Override
    public IntArray map(int index) {
        if (index!=0) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return toIntArray();
    }
}

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
package bref4;

import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.VcfRecGTParser;

/**
 * <p>Class {@code AllelicBrefRec} is a {@code BrefRec} instance.</p>
 *
 * <p>Instances of class {@code AllelicBrefRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AllelicBrefRec implements BrefRec {

    private final Marker marker;
    private final int size;
    private final int[][] alleleToHaps;
    private final int nullAllele;
    private final int nonNullAlleleCnt;

    /**
     * Constructs a new {@code AllelicBrefRec} instance from the specified
     * data.
     * @param rec a VCF record with phased, non-missing alleles.
     * @throws NullPointerException if {@code rec == null}
     */
    public AllelicBrefRec(RefGTRec rec) {
        this.marker = rec.marker();
        this.size = rec.size();
        this.alleleToHaps = rec.alleleToHaps();
        this.nullAllele = BrefRec.nullRow(alleleToHaps);
        this.nonNullAlleleCnt = BrefRec.nonNullCnt(alleleToHaps);
    }

    /**
     * Constructs a new {@code AllelicBrefRec} instance from the specified
     * data.
     * @param gtp a VCF record with phased, non-missing alleles.
     * @throws NullPointerException if {@code gtp == null}
     */
    public AllelicBrefRec(VcfRecGTParser gtp) {
        this.marker = gtp.marker();
        this.size = gtp.nSamples()<<1;
        this.alleleToHaps = gtp.nonMajRefIndices();
        this.nullAllele = BrefRec.nullRow(alleleToHaps);
        this.nonNullAlleleCnt = BrefRec.nonNullCnt(alleleToHaps);
    }

    /**
     * Constructs a new {@code AllelicBrefRec} instance from the specified
     * data. The contract for this constructor is unspecified if there exists
     * haplotypes {@code h1} and {@code h2} such that
     * {@code (0 <= h1 && h1 < h2 && h2 < this.size())} and
     * {@code (rec.get(h1) != rec.get(h2))} and
     * {@code (hapToSeq.get(h1) == hapToSeq.get(h2))}
     * @param hapToSeq a map from haplotype to sequence
     * @param rec the {@code AllelicBrefRec} whose haplotype indices will be
     * mapped
     * @throws IllegalArgumentException if {@code hapToSeq.size() > rec.size()}
     * @throws NullPointerException if {@code hapToSeq == null || rec == null}
     */
    private AllelicBrefRec(IndexArray hapToSeq, AllelicBrefRec rec) {
        if (hapToSeq.size() != rec.size()) {
            throw new IllegalArgumentException(String.valueOf(hapToSeq.size()));
        }
        int[][] newAlleleToHaps = new int[rec.alleleToHaps.length][];
        for (int j=0; j<rec.alleleToHaps.length; ++j) {
            if (rec.alleleToHaps[j]!=null) {
                newAlleleToHaps[j] = Arrays.stream(rec.alleleToHaps[j])
                        .map(i -> hapToSeq.get(i))
                        .distinct()
                        .sorted()
                        .toArray();
            }
        }
        this.marker = rec.marker();
        this.size = hapToSeq.valueSize();
        this.nullAllele = rec.nullAllele;
        this.alleleToHaps = newAlleleToHaps;
        this.nonNullAlleleCnt = BrefRec.nonNullCnt(alleleToHaps);
    }

    @Override
    public long estBytes() {
        int nNonNullAlleles = alleleToHaps.length - 1;
        long estBytes= 12*(nNonNullAlleles + 1);    // 12 bytes overhead per "owned" object
        estBytes += 12;                             // 4 bytes per int field
        estBytes += 8*(2 + alleleToHaps.length);    // 8 bytes per reference
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (j!=nullAllele) {
                estBytes += 4*alleleToHaps[j].length;  // 4 bytes per integer
                estBytes += 4;                         // 4 bytes to store inner array length
            }
        }
        estBytes += 4; // 4 bytes to store outer array length
        estBytes += 4; // assume 4 bytes for alignment
        return estBytes;
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public int get(int hap) {
        if (hap < 0 || hap >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (j != nullAllele) {
                if (Arrays.binarySearch(alleleToHaps[j], hap) >= 0) {
                    return j;
                }
            }
        }
        return nullAllele;
    }

    @Override
    public int[][] alleleToHaps() {
        int[][] ia = new int[alleleToHaps.length][];
        for (int j=0; j<ia.length; ++j) {
            if (j!=nullAllele) {
                ia[j] = alleleToHaps[j].clone();
            }
        }
        return ia;
    }

    @Override
    public IndexArray hapToAllele() {
        int[] ia = new int[size];
        Arrays.fill(ia, nullAllele);
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (j!=nullAllele) {
                for (int k : alleleToHaps[j]) {
                    ia[k] = j;
                }
            }
        }
        int nAlleles = marker.nAlleles();
        IntArray intArray = IntArray.packedCreate(ia, nAlleles);
        return new IndexArray(intArray, nAlleles);
    }

    @Override
    public int nullRow() {
        return nullAllele;
    }

    @Override
    public int nStoredHapIndices() {
        return nonNullAlleleCnt;
    }

    @Override
    public BrefRec applyMap(IndexArray hapToSeq) {
        return new AllelicBrefRec(hapToSeq, this);
    }
}

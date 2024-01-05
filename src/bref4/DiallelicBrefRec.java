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
package bref4;

import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.VcfRecGTParser;

/**
 * <p>Class {@code DiallelicBrefRec} is a {@code BrefRec} instance
 * that has two alleles.</p>
 *
 * <p>Instances of class {@code DiallelicBrefRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class DiallelicBrefRec implements BrefRec {

    private final Marker marker;
    private final int size;
    private final int nullAllele;
    private final int nonNullAllele;
    private final int[] nonNullHapIndices;

    /**
     * Constructs a new {@code DiallelicAlleleCodedBrefRec} instance from
     * the specified data.
     * @param rec a VCF record with phased, non-missing alleles
     * @throws IllegalArgumentException if {@code rec.marker().nAlleles() != 2}
     * @throws NullPointerException if {@code rec == null}
     */
    public DiallelicBrefRec(RefGTRec rec) {
        if (rec.marker().nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(rec.marker()));
        }
        this.marker = rec.marker();
        this.size = rec.size();
        int[][] alleleToHaps = rec.alleleToHaps();
        this.nullAllele = alleleToHaps[0]==null ? 0 : 1;
        this.nonNullAllele = (1 - nullAllele);
        this.nonNullHapIndices = alleleToHaps[nonNullAllele];
    }

    /**
     * Constructs a new {@code DiallelicAlleleCodedBrefRec} instance from
     * the specified data.
     * @param gtp a VCF record with phased, non-missing alleles
     * @throws IllegalArgumentException if {@code gtp.nAlleles() != 2}
     * @throws NullPointerException if {@code gtp == null}
     */
    public DiallelicBrefRec(VcfRecGTParser gtp) {
        if (gtp.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(gtp.nAlleles()));
        }
        this.marker = gtp.marker();
        this.size = 2*gtp.nSamples();
        int[][] alleleToHaps = gtp.nonMajRefIndices();
        this.nullAllele = alleleToHaps[0]==null ? 0 : 1;
        this.nonNullAllele = (1 - nullAllele);
        this.nonNullHapIndices = alleleToHaps[nonNullAllele];
    }

    /**
     * Constructs a new {@code DiallelicBrefRec} instance from the specified
     * data. The contract for this constructor is unspecified if there exists
     * haplotypes {@code h1} and {@code h2} such that
     * {@code (0 <= h1 && h1 < h2 && h2 < this.size())} and
     * {@code (rec.get(h1) != rec.get(h2))} and
     * {@code (hapToSeq.get(h1) == hapToSeq.get(h2))}.
     * @param hapToSeq a map from haplotype to sequence
     * @param rec the {@code DiallelicBrefRec} whose haplotype indices 
     * will be mapped
     * @throws IllegalArgumentException if {@code hapToSeq.size() > rec.size()}
     * @throws NullPointerException if {@code hapToSeq == null || rec == null}
     */    
    private DiallelicBrefRec(IndexArray hapToSeq,
            DiallelicBrefRec rec) {
        if (hapToSeq.size() != rec.size()) {
            throw new IllegalArgumentException(String.valueOf(hapToSeq.size()));
        }
        this.marker = rec.marker();
        this.size = hapToSeq.valueSize();
        this.nullAllele = rec.nullAllele;
        this.nonNullAllele = rec.nonNullAllele;
        this.nonNullHapIndices = Arrays.stream(rec.nonNullHapIndices)
                        .map(i -> hapToSeq.get(i))
                        .distinct()
                        .sorted()
                        .toArray();
    }

    @Override
    public long estBytes() {
        long estBytes = 12;                             // 12 bytes overhead per "owned" object
        estBytes += 12;                                 // 4 bytes per int field
        estBytes += 16;                                 // 8 bytes per reference
        estBytes += 4*nonNullHapIndices.length;         // 4 bytes per integer
        estBytes += 4;                                  // 4 bytes to store array length
        estBytes += 4;                                  // assume 4 bytes for alignment
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
        if (Arrays.binarySearch(nonNullHapIndices, hap) >= 0) {
            return nonNullAllele;
        }
        else {
            return nullAllele;
        }
    }

    @Override
    public int[][] alleleToHaps() {
        int[][] ia = new int[2][];
        ia[nonNullAllele] = nonNullHapIndices.clone();
        return ia;
    }

    @Override
    public IndexArray hapToAllele() {
        int[] ia = new int[size];
        Arrays.fill(ia, nullAllele);
        for (int j : nonNullHapIndices) {
            ia[j] = nonNullAllele;
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
        return nonNullHapIndices.length;
    }

    @Override
    public BrefRec applyMap(IndexArray hapToSeq) {
        return new DiallelicBrefRec(hapToSeq, this);
    }
}

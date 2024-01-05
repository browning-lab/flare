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
package vcf;

import bref4.BrefRec;
import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;

/**
 * <p>Class {@code IntArrayRefGT} represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.
 * Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>Instances of class {@code IntArrayRefGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IntArrayRefGTRec implements RefGTRec {

    private final Marker marker;
    private final Samples samples;
    private final IntArray alleles;

    @Override
    public long estBytes() {
        int overhead = 12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 3*8;  // assume 8 bytes per reference
        estBytes += 4*(1 + alleles.size()); // 4 bytes to store array length
        return estBytes;
    }

    /**
     * Constructs a new {@code IntArrayRefGT} instance from the
     * specified phased genotype data.
     * @param marker the marker
     * @param samples the list of samples
     * @param alleles the list of phased alleles
     * @throws IllegalArgumentException if
     * {@code alleles.length != 2*samples.size()}
     * @throws IllegalArgumentException if
     * {@code (alleles[j] < 0 || alleles[j] >= marker.nAlleles())} for any
     * {@code j} satisfying {@code 0 <= j && j < allele1.length}
     * @throws NullPointerException if
     * {@code marker==null || samples==null || alleles==null}
     */
    public IntArrayRefGTRec(Marker marker, Samples samples, int[] alleles) {
        if (alleles.length != 2*samples.size()) {
            throw new IllegalArgumentException(String.valueOf(alleles.length));
        }
        int nAlleles = marker.nAlleles();
        for (int a : alleles) {
            if (a<0 || a>=nAlleles) {
                throw new IllegalArgumentException(String.valueOf(a));
            }
        }
        this.marker = marker;
        this.samples = samples;
        this.alleles = IntArray.packedCreate(alleles, marker.nAlleles());
    }

    /**
     * Constructs a new {@code IntArrayRefGT} instance from the
     * specified phased genotype data.
     * @param marker the marker
     * @param samples the list of samples
     * @param alleles the list of phased alleles
     * @throws IllegalArgumentException if
     * {@code alleles.size() != 2*samples.size()}
     * @throws IllegalArgumentException if
     * {@code alleles.valueSize() >= marker.nAlleles()}
     * @throws NullPointerException if
     * {@code marker==null || samples==null || alleles==null}
     */
    public IntArrayRefGTRec(Marker marker, Samples samples, IndexArray alleles) {
        if (alleles.size() != 2*samples.size()) {
            throw new IllegalArgumentException(String.valueOf(alleles.size()));
        }
        if (alleles.valueSize() > marker.nAlleles()) {
            throw new IllegalArgumentException(String.valueOf(alleles.valueSize()));
        }
        this.marker = marker;
        this.samples = samples;
        this.alleles = alleles;
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
    public boolean isAlleleCoded() {
        return false;
    }

    @Override
    public int majorAllele() {
        return majorAllele(alleleCounts());
    }

    private int majorAllele(int[] alCnts) {
        int majAl = 0;
        for (int al=1; al<alCnts.length; ++al) {
            if (alCnts[al]>alCnts[majAl]) {
                majAl = al;
            }
        }
        return majAl;
    }

    @Override
    public int[] alleleCounts() {
        int[] alCnts = new int[marker.nAlleles()];
        for (int h=0, n=size(); h<n; ++h) {
            ++alCnts[get(h)];
        }
        return alCnts;
    }

    @Override
    public int alleleCount(int allele) {
        int[] alCnts = alleleCounts();
        if (allele==majorAllele(alCnts)) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return alCnts[allele];
        }
    }

    /**
     * Returns the samples. The returned samples are the filtered samples
     * after all sample exclusions.
     *
     * @return the samples.
     */
    @Override
    public Samples samples() {
        return samples;
    }


    @Override
    public int size() {
        return 2*samples.size();
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public int[][] alleleToHaps() {
        int[] alCnts = alleleCounts();
        int majAllele = 0;
        for (int al=1; al<alCnts.length; ++al) {
            if (alCnts[al] > alCnts[majAllele]) {
                majAllele = al;
            }
        }
        int[][] hapIndices = new int[alCnts.length][];
        for (int al=0; al<alCnts.length; ++al) {
            if (al!=majAllele) {
                hapIndices[al] = new int[alCnts[al]];
            }
        }
        Arrays.fill(alCnts, 0);
        for (int h=0, n=size(); h<n; ++h) {
            int al = get(h);
            if (al!=majAllele) {
                hapIndices[al][alCnts[al]++] = h;
            }
        }
        return hapIndices;
    }

    @Override
    public IndexArray hapToAllele() {
        return new IndexArray(alleles, marker.nAlleles());
    }

    @Override
    public int nAlleleCodedHaps() {
        return BrefRec.nonNullCnt(alleleToHaps());
    }

    @Override
    public int get(int hap) {
        return alleles.get(hap);
    }

    @Override
    public int hapIndex(int allele, int copy) {
        int[][] hapIndices = alleleToHaps();
        if (hapIndices[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return hapIndices[allele][copy];
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
        return new IntArray[] {alleles};
    }

    @Override
    public IntArray map(int index) {
        if (index!=0) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return alleles;
    }
}

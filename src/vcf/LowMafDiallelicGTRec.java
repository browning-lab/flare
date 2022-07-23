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
 * <p>Class {@code LowMafDiallelicGTRc} stores genotypes for a list of samples
 * at a diallelic marker. Instances of class {@code LowMafGTRec} store lists of
 * haplotypes carrying the minor or missing allele. All genotypes are
 * considered to be unphased if any sample has an unphased or missing
 * genotype.</p>
 *
 * <p>Instances of class {@code LowMafDiallelicGTRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LowMafDiallelicGTRec implements GTRec {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int majorAllele;
    private final int minorAllele;
    private final int[] missingSamples;
    private final int[] minorAlleles;
    private final boolean isPhased;

    @Override
    public long estBytes() {
        int overhead = 2*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 4*8;  // assume 8 bytes per reference
        estBytes += 16; // 4 bytes per int, assume 4 bytes per boolean due to alignment
        estBytes += 4*(missingSamples.length + minorAlleles.length);
        return estBytes;
    }

    /**
     * Constructs a new {@code LowMafDiallelicGTRec} representing the specified
     * VCF record's GT format field data.
     * @param listRep the VCF record genotype data
     * @throws IllegalArgumentException if
     * {@code listRep.marker().nAlleles() != 2}
     * @throws NullPointerException if {@code listRep == null}
     */
    public LowMafDiallelicGTRec(VcfRecGTParser.HapListRep listRep) {
        if (listRep.marker().nAlleles()!=2) {
            throw new IllegalArgumentException(
                    String.valueOf(listRep.marker().nAlleles()));
        }
        this.marker = listRep.marker();
        this.samples = listRep.samples();
        this.nHaps = samples.size()<<1;
        this.majorAllele = listRep.majorAllele();
        this.minorAllele = 1 - listRep.majorAllele();
        this.minorAlleles = listRep.hapLists(true)[minorAllele];
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
        if (Arrays.binarySearch(minorAlleles, hap) >= 0) {
            return minorAllele;
        }
        else if (Arrays.binarySearch(missingSamples, (hap>>1)) >= 0) {
            return -1;
        }
        else {
            return majorAllele;
        }
    }

    @Override
    public int[] alleles() {
        int[] ia = IntStream.range(0, nHaps)
                .map(h -> majorAllele)
                .toArray();
        for (int h : minorAlleles) {
            ia[h] = minorAllele;
        }
        for (int s : missingSamples) {
            int h1 = s <<1;
            int h2 = h1 | 0b1;
            ia[h1] = -1;
            ia[h2] = -1;
        }
        return ia;
    }

    /**
     * Returns the major allele.  If both alleles have the same allele count,
     * the reference allele is returned.
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
            return nHaps - minorAlleles.length - (missingSamples.length<<1);
        }
        else {
            return minorAlleles.length;
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

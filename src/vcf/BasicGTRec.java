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

/**
 * <p>Class {@code BasicGTRec} stores genotypes for a list of samples
 * at a single marker. The phased or unphased status of each genotype is
 * stored.</p>
 *
 * <p>Instances of class {@code BasicGTRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicGTRec implements GTRec {

    private final Marker marker;
    private final Samples samples;
    private final int[] alleles;
    private final boolean[] isPhased;
    private final boolean allPhased;

    @Override
    public long estBytes() {
        int overhead = 2*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 4*8;  // assume 8 bytes per reference
        estBytes += 4*alleles.length + isPhased.length;
        estBytes += 8; // 4 bytes per array to store length
        estBytes += 4; // assume 4 bytes for this.allPhased due to alignment
        return estBytes;
    }

    /**
     * Constructs a new {@code BasicGTRec} instance representing
     * the specified VCF record's GT format field data.
     *
     * @param recParser the VCF record genotype data
     * @throws IllegalArgumentException if a format error is detected
     * in the VCF record
     * @throws NullPointerException if {@code recParser == null}
     */
    public BasicGTRec(VcfRecGTParser recParser) {
        int nSamples = recParser.samples().size();
        int nHaps = nSamples<<1;
        int[] alleles0 = new int[nHaps];
        boolean[] isPhased0 = new boolean[nSamples];
        this.allPhased = recParser.storeAlleles(alleles0, isPhased0);
        this.marker = recParser.marker();
        this.samples = recParser.samples();
        this.alleles = alleles0;
        this.isPhased = isPhased0;
    }


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
    public boolean isPhased() {
        return allPhased;
    }

    @Override
    public boolean isPhased(int sample) {
        return isPhased[sample];
    }


    @Override
    public int allele1(int sample) {
        return alleles[sample<<1];
    }

    @Override
    public int allele2(int sample) {
        return alleles[(sample<<1) | 0b1];
    }

    @Override
    public int get(int hap) {
        return alleles[hap];
    }

    @Override
    public int[] alleles() {
        return alleles.clone();
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

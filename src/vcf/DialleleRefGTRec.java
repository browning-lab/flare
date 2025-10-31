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

import blbutil.BooleanArray;
import blbutil.Const;
import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code DialleleRefGTRec} represent represents phased,
 * non-missing genotypes for a list of reference samples at a single diallelic
 * marker.</p>
 *
 * <p>Class {@code DialleleRefGTRec} stores haplotypes that carry
 * the minor allele.</p>
 *
 * <p>Instances of class {@code DialleleRefGTRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class DialleleRefGTRec implements RefGTRec {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int nullRow;
    private final int nonNullRow;
    private final int[] nonNullRowHaps;

    @Override
    public long estBytes() {
        int overhead = 12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 3*8;  // assume 8 bytes per reference
        estBytes += 12; // 4 bytes per int
        estBytes += 4*(1 + nonNullRowHaps.length); // 4 bytes to store array length
        return estBytes;
    }

    /**
     * Constructs a new {@code TwoAlleleRefGTRec} instance from the
     * specified data.
     *
     * @param rec the phased, non-missing genotype data
     * @throws IllegalArgumentException if {@code rec.marker().nAlleles() != 2}
     * @throws NullPointerException if {@code rec == null}
     */
    public DialleleRefGTRec(RefGTRec rec) {
        if (rec.marker().nAlleles()!=2) {
            throw new IllegalArgumentException(
                    String.valueOf(rec.marker().nAlleles()!=2));
        }
        int[][] hapIndices = rec.alleleToHaps();
        int majAllele = 0;
        while (hapIndices[majAllele]!=null) {
            ++majAllele;
        }
        this.marker = rec.marker();
        this.samples = rec.samples();
        this.nHaps = rec.size();
        this.nullRow = majAllele;
        this.nonNullRow = 1 - majAllele;
        this.nonNullRowHaps = hapIndices[nonNullRow];
    }

    /**
     * Constructs a new {@code TwoAlleleRefGTRec} instance from the
     * specified data.
     *
     * @param gtp a VCF record parser that extracts sample genotypes
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased genotype or missing allele
     * @throws IllegalArgumentException if {@code gtp.nAlleles() != 2}
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws NullPointerException if {@code gtp == null}
     */
    public DialleleRefGTRec(VcfRecGTParser gtp) {
        if (gtp.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(gtp.nAlleles()));
        }
        int[][] nonMajIndices = gtp.nonMajAlleleIndices();
        this.marker = gtp.marker();
        this.samples = gtp.samples();
        this.nHaps = 2*gtp.nSamples();
        this.nullRow = nonMajIndices[0]==null ? 0 : 1;
        this.nonNullRow = 1 - nullRow;
        this.nonNullRowHaps = nonMajIndices[nonNullRow];
    }

    /**
     * Constructs a new {@code AlleleRefGTRec} instance from the specified
     * data. The {@code isAlleleRecord()} method of the returned
     * {@code RefGTRec} instance will return {@code true}. The contract for
     * this method is undefined if any two integers in
     * {@code alleleToHaps} are equal.
     * @param marker the marker
     * @param samples the samples
     * @param alleleToHaps an array of length {@code marker.nAlleles()} with
     * a unique {@code null} element and whose {@code j}-th element
     * is either {@code null} or an {@code int[]} whose elements are an
     * increasing list of the haplotypes that carry the {@code j}-th allele
     *
     * @throws IllegalArgumentException if {@code marker.nAlleles() != 2}
     * @throws IllegalArgumentException if
     * {@code marker.nAlleles() != alleleToHaps.length}
     * @throws IllegalArgumentException {@code alleleToHaps} does not have
     * a unique {@code null} element
     * @throws IllegalArgumentException if a non-null element of
     * {@code alleleToHaps} is not an increasing list of distinct haplotype
     * indices between 0 (inclusive) and {@code 2*samples.size()} (exclusive)
     * @throws NullPointerException if
     * {@code ((marker == null) || (samples == null) || (alleleToHaps == null))}
     */
    public DialleleRefGTRec(Marker marker, Samples samples,
            int[][] alleleToHaps) {
        if (marker.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(marker.nAlleles()));
        }
        if (alleleToHaps.length!=2) {
            throw new IllegalArgumentException(String.valueOf(alleleToHaps.length));
        }
        this.marker = marker;
        this.samples = samples;
        this.nHaps = 2*samples.size();
        this.nullRow = AlleleRefGTRec.checkAlleleToHaps(alleleToHaps, nHaps);
        this.nonNullRow = 1 - nullRow;
        this.nonNullRowHaps = alleleToHaps[nonNullRow].clone();
    }

    @Override
    public int[][] alleleToHaps() {
        int[][] hapIndices = new int[2][];
        hapIndices[nonNullRow] = nonNullRowHaps.clone();
        return hapIndices;
    }

    @Override
    public IndexArray hapToAllele() {
        return new IndexArray(toIntArray(), 2);
    }

    private IntArray toIntArray() {
        int[] ia = IntStream.range(0, nHaps)
                .map(h -> nullRow)
                .toArray();
        for (int h : nonNullRowHaps) {
            ia[h] = nonNullRow;
        }
        return IntArray.packedCreate(ia, 2);
    }

    @Override
    public int nullRowAlleleCnt() {
        return nonNullRowHaps.length;
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
        if (Arrays.binarySearch(nonNullRowHaps, hap) >= 0) {
            return nonNullRow;
        }
        else {
            return nullRow;
        }
    }

    @Override
    public boolean isAlleleRecord() {
        return true;
    }

    @Override
    public int nullRow() {
        return nullRow;
    }

    @Override
    public int[] alleleCounts() {
        int[] alCnts = new int[2];
        alCnts[nullRow] = nHaps - nonNullRowHaps.length;
        alCnts[nonNullRow] = nonNullRowHaps.length;
        return alCnts;
    }

    @Override
    public int alleleCount(int allele) {
        if (allele==nullRow) {
            return nHaps - nonNullRowHaps.length;
        }
        else if (allele==nonNullRow) {
            return nonNullRowHaps.length;
        }
        else {
            throw new IndexOutOfBoundsException(String.valueOf(allele));
        }
    }

    @Override
    public int nonNullRowHap(int allele, int copy) {
        if (allele==nullRow) {
            throw new IllegalArgumentException("null row");
        }
        else {
            return nonNullRowHaps[copy];
        }
    }

    @Override
    public boolean isCarrier(int allele, int hap) {
        return get(hap)==allele;
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

    @Override
    public String toString() {
        return toVcfRecord();
    }

    @Override
    public String toVcfRecord() {
        StringBuilder sb = new StringBuilder();
        MarkerUtils.appendFirst8Fields(marker, sb);  // no INFO/{AN,AC} update
//        MarkerUtils.appendFirst8Fields(marker, nHaps, alleleCounts(), sb);
        sb.append("\tGT");
        int index = 0;
        int nextMinorHap = (index<nonNullRowHaps.length) ? nonNullRowHaps[index++] : nHaps;
        for (int h=0; h<nHaps; ++h) {
            sb.append((h & 0b1)==0 ? Const.tab : Const.phasedSep);
            if (h==nextMinorHap) {
                sb.append(nonNullRow);
                nextMinorHap = (index<nonNullRowHaps.length) ? nonNullRowHaps[index++] : nHaps;
            }
            else {
                sb.append(nullRow);
            }
        }
        return sb.toString();
    }

    @Override
    public String toVcfRecord(BooleanArray isHaploid) {
        if (isHaploid==null) {
            return this.toString();
        }
        else {
            return toString0(isHaploid);
        }
    }

    private String toString0(BooleanArray isHaploid) {
        if (isHaploid.size()!=this.samples().size()) {
            throw new IllegalArgumentException(String.valueOf(isHaploid));
        }
        int[] indices = nonNullRowHaps;
        StringBuilder sb = new StringBuilder();
        MarkerUtils.appendFirst8Fields(marker, sb);  // no INFO/{AN,AC} update
        sb.append("\tGT");
        int index = 0;
        int nextMinorHap = (index<indices.length) ? indices[index++] : nHaps;
        for (int h=0; h<nHaps; ++h) {
            if ((h & 0b1)==0 || isHaploid.get(h>>1)==false) {
                sb.append((h & 0b1)==0 ? Const.tab : Const.phasedSep);
                if (h==nextMinorHap) {
                    sb.append(nonNullRow);
                    nextMinorHap = (index<indices.length) ? indices[index++] : nHaps;
                }
                else {
                    sb.append(nullRow);
                }
            }
        }
        return sb.toString();
    }
}

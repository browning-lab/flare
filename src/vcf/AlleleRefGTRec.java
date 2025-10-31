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
 * <p>Class {@code AlleleRefGTRec} represent represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.</p>
 *
 * <p>Class {@code AlleleRefGTRec} stores the haplotypes that carry each
 * for all alleles except one.</p>
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
        estBytes += 8;  // to store nHaps and nullRow
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
        this.alleleToHaps = gtp.nonMajAlleleIndices();
        int majAllele = 0;
        while (alleleToHaps[majAllele]!=null) {
            ++majAllele;
        }
        this.majorAllele = majAllele;
    }

    /**
     * Constructs a new {@code AlleleRefGTRec} instance from the specified
     * data. The {@code isAlleleRecord()} method of the returned
     * {@code RefGTRec} instance will return {@code true}. The contract for
     * this method is undefined if any two integers in {@code alleleToHaps}
     * are equal.
     * @param marker the marker
     * @param samples the samples
     * @param alleleToHaps an array of length {@code marker.nAlleles()} with
     * a unique {@code null} element and whose {@code j}-th element
     * is either {@code null} or an {@code int[]} whose elements are an
     * increasing list of the haplotypes that carry the {@code j}-th allele
     *
     * @throws IllegalArgumentException if
     * {@code marker.nAlleles() != alleleToHaps.length}
     * @throws IllegalArgumentException {@code alleleToHaps} does not
     * have a unique {@code null} element
     * @throws IllegalArgumentException if any non-null element of
     * {@code alleleToHaps} is not a sorted list of distinct haplotype indices
     * between 0 (inclusive) and {@code 2*samples.size()} (exclusive)
     * @throws NullPointerException if
     * {@code ((marker == null) || (samples == null) || (alleleToHaps == null))}
     */
    public AlleleRefGTRec(Marker marker, Samples samples, int[][] alleleToHaps) {
        if (marker.nAlleles()!=alleleToHaps.length) {
            throw new IllegalArgumentException(String.valueOf(alleleToHaps.length));
        }
        this.marker = marker;
        this.samples = samples;
        this.nHaps = 2*samples.size();
        this.majorAllele = checkAlleleToHaps(alleleToHaps, nHaps);
        this.alleleToHaps = deepCopy(alleleToHaps);
    }

    /**
     * Checks that the {@code alleleToHaps} array has a unique {@code null}
     * element and that each non-{@code null} element of {@code alleleToHaps}
     * is a list of increasing integers bounded between 0 (inclusive) and
     * {@code nHaps} (exclusive). The contract for this method is undefined
     * if any two integers in {@code alleleToHaps} are equal.
     *
     * @param alleleToHaps an array of length {@code marker.nAlleles()} with
     * a unique {@code null} element and whose {@code j}-th element
     * is either {@code null} or an {@code int[]} whose elements are an
     * increasing list of the haplotypes that carry the {@code j}-th allele
     * @param nHaps the number of haplotypes
     * @return the index of the {@code null} element of {@code alleleToHaps}
     * @throws IllegalArgumentException {@code alleleToHaps} does not
     * have a unique {@code null} element
     * @throws IllegalArgumentException if any non-{@code null}
     * element of {@code alleleToHaps} is not sorted in increasing order
     * @throws IllegalArgumentException if any integer in
     * {@code alleleToHaps} is negative or greater than or equal to
     * {@code nHaps}
     * @throws NullPointerException if {@code (alleleToHaps == null)}
     */
    static int checkAlleleToHaps(int[][] alleleToHaps, int nHaps) {
        int nullRow = -1;
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (alleleToHaps[j]==null) {
                if (nullRow == -1) {
                    nullRow = j;
                }
                else {
                    throwArrayError();
                }
            }
            else {
                checkSorted(alleleToHaps[j], nHaps);
            }
        }
        if (nullRow == -1) {
            throwArrayError();
        }
        return nullRow;
    }

    private static void throwArrayError() {
        throw new IllegalArgumentException("invalid array");
    }

    private static void checkSorted(int[] ia, int nHaps) {
        if (ia.length>0 && (ia[0] < 0 || ia[ia.length - 1] >= nHaps)) {
            throwArrayError();
        }
        for (int k=1; k<ia.length; ++k) {
            if (ia[k-1] > ia[k]) {
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
    public int nullRowAlleleCnt() {
        return RefGTRec.nonNullCnt(alleleToHaps);
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
    public boolean isAlleleRecord() {
        return true;
    }

    @Override
    public int nullRow() {
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
            return alleleCounts()[allele];
        }
        else {
            return alleleToHaps[allele].length;
        }
    }

    @Override
    public int nonNullRowHap(int allele, int copy) {
        if (alleleToHaps[allele]==null) {
            throw new IllegalArgumentException("null row");
        }
        else {
            return alleleToHaps[allele][copy];
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
        MarkerUtils.appendFirst8Fields(marker, sb); // no INFO/{AN,AC} update
//        MarkerUtils.appendFirst8Fields(marker, nHaps, alleleCounts(), sb);
        sb.append("\tGT");
        int[] minorIndices = sortedMinorIndices();
        int[] minorAlleles = minorAlleles(minorIndices);
        int index = 0;
        int nextMinorHap = (index<minorIndices.length) ? minorIndices[index] : nHaps;
        for (int h=0; h<nHaps; ++h) {
            sb.append((h & 0b1)==0 ? Const.tab : Const.phasedSep);
            if (h==nextMinorHap) {
                sb.append(minorAlleles[index++]);
                nextMinorHap = (index<minorAlleles.length) ? minorIndices[index] : nHaps;
            }
            else {
                sb.append(majorAllele);
            }
        }
        return sb.toString();
    }

    @Override
    public String toVcfRecord(BooleanArray isHaploid) {
        if (isHaploid==null) {
            return toVcfRecord();
        }
        else {
            return toString0(isHaploid);
        }
    }

    private String toString0(BooleanArray isHaploid) {
        if (isHaploid.size()!=this.samples().size()) {
            throw new IllegalArgumentException(String.valueOf(isHaploid));
        }
        StringBuilder sb = new StringBuilder();
        MarkerUtils.appendFirst8Fields(marker, sb); // no INFO/{AN,AC} update
        sb.append("\tGT");
        int[] minorIndices = sortedMinorIndices();
        int[] minorAlleles = minorAlleles(minorIndices);
        int index = 0;
        int nextMinorHap = (index<minorIndices.length) ? minorIndices[index] : nHaps;
        for (int h=0; h<nHaps; ++h) {
            if ((h & 0b1)==0 || isHaploid.get(h>>1)==false) {
                sb.append((h & 0b1)==0 ? Const.tab : Const.phasedSep);
                if (h==nextMinorHap) {
                    sb.append(minorAlleles[index++]);
                    nextMinorHap = (index<minorAlleles.length) ? minorIndices[index] : nHaps;
                }
                else {
                    sb.append(majorAllele);
                }
            }
        }
        return sb.toString();
    }

    private int[] sortedMinorIndices() {
        int cnt = 0;
        for (int[] haps : alleleToHaps) {
            if (haps!=null) {
                cnt += haps.length;
            }
        }
        int[] minorIndices = new int[cnt];
        int start = 0;
        for (int[] haps : alleleToHaps) {
            if (haps!=null) {
                System.arraycopy(haps, 0, minorIndices, start, haps.length);
                start += haps.length;
            }
        }
        assert start==minorIndices.length;
        Arrays.sort(minorIndices);
        return minorIndices;
    }

    private int[] minorAlleles(int[] minorIndices) {
        int[] minorAlleles = new int[minorIndices.length];
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (j!=majorAllele) {
                int[] haps = alleleToHaps[j];
                int index = 0;
                for (int hap : haps) {
                    index = Arrays.binarySearch(minorIndices, index, minorIndices.length, hap);
                    assert index>=0;
                    minorAlleles[index] = j;
                }
            }
        }
        return minorAlleles;
    }
}

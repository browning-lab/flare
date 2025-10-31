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

/**
 * <p>Interface {@code RefGTRec} represents represents phased genotype data
 * for one VCF record. For implementations of this interface, unless otherwise
 * specified in the implementation documentation, if the
 * {@code isAlleleRecord()} method returns {@code false}, the {@code nullRow()},
 * {@code alleleCount()}, and {@code nonNullRowHap()} methods will be
 * computationally expensive with compute time proportional to the number of
 * haplotypes. Alternatively, if the {@code isAlleleRecord()} method
 * returns {@code true}, the {@code maps()} and {@code map()} methods will be
 * computationally expensive with compute time proportional to the number
 * of haplotypes.
 * </p>
 * <p>Memory efficiency is greatest if {@code this.nullRow()} returns the
 * major allele, but it is not required that {@code this.nullRow()} return
 * the major allele because sample filtering can change the major allele.
 * </p>
 * <p>Instances of {@code RefGTRec} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface RefGTRec extends GTRec {

    /**
     * Constructs and returns a new {@code RefGTRec} instance from the specified
     * data. The {@code isNonMajorAlleleRecord()} method of the returned
     * {@code RefGTRec} instance will return {@code true}.
     * @param rec the phased, non-missing genotype data
     * @return an allele-coded {@code RefGTRec} instance for the
     * specified data
     * @throws NullPointerException if {@code rec == null}
     */
    static RefGTRec alleleRefGTRec(RefGTRec rec) {
        if (rec.isAlleleRecord()) {
            return rec;
        }
        if (rec.marker().nAlleles()==2) {
            return new DialleleRefGTRec(rec);
        }
        else {
            return new AlleleRefGTRec(rec);
        }
    }

    /**
     * Constructs and returns a new {@code RefGTRec} instance from the
     * specified data. The {@code isAlleleRecord()} method of the returned
     * {@code RefGTRec} instance will return {@code true}.
     *
     * @param gtp a VCF record parser that extracts sample genotypes
     * @return a new {@code RefGTRec} instance
     *
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased genotype or missing allele
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws NullPointerException if {@code (gtp == null)}
     */
    static RefGTRec alleleRefGTRec(VcfRecGTParser gtp) {
        if (gtp.nAlleles()==2) {
            return new DialleleRefGTRec(gtp);
        }
        else {
            return new AlleleRefGTRec(gtp);
        }
    }

    /**
     * Constructs and returns a new {@code RefGTRec} instance from the specified
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
     * @return an allele-coded {@code RefGTRec} instance
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
    static RefGTRec alleleRefGTRec(Marker marker, Samples samples,
            int[][] alleleToHaps) {
        if (marker.nAlleles()==2) {
            return new DialleleRefGTRec(marker, samples, alleleToHaps);
        }
        else {
            return new AlleleRefGTRec(marker, samples, alleleToHaps);
        }
    }

    /**
     * Returns {@code true} if this instance stores the indices of haplotypes
     * that carry each allele {@code j} where {@code (j != this.nullRow())},
     * and returns {@code false} otherwise.
     *
     * @return {@code true} if this instance stores the indices of haplotypes
     * that carry each allele {@code j} where {@code (j != this.nullRow())}
     */
    boolean isAlleleRecord();

    /**
     * Returns the sum of the lengths of non-null rows of the specified
     * two-dimensional array.
     * @param alleleToHaps a two-dimensional array
     * @return the sum of the lengths of non-null rows
     * @throws NullPointerException if {@code IalleleToHaps == null)}
     */
    public static int nonNullCnt(int[][] alleleToHaps) {
        int nonNullCnt=0;
        for (int[] ia : alleleToHaps) {
            if (ia!=null) {
                nonNullCnt += ia.length;
            }
        }
        return nonNullCnt;
    }

    /**
     * Returns the smallest row index {@code j} such that
     * {@code (allelesToHap[j] == null)} or {@code -1}
     * if no such allele index exists.
     *
     * @param alleleToHaps a two-dimensional array
     * @return the smallest row index {@code j} such that
     * {@code (allelesToHap[j] == null)} or {@code -1}
     * if no such allele index exists
     * @throws NullPointerException if {@code alleleToHaps == null}
     */
    public static int nullRow(int[][] alleleToHaps) {
        for (int j=0; j<alleleToHaps.length; ++j) {
            if (alleleToHaps[j]==null) {
                return j;
            }
        }
        return -1;
    }

    /**
     * Returns the index of the unique {@code null} element of
     * {@code this.alleleToHaps()}.
     *
     * @return the index of the unique {@code null} element of
     * {@code this.alleleToHaps()}
     */
    int nullRow();

    /**
     * Returns an array of length {@code marker.nAlleles()} with
     * a unique {@code null} element and whose {@code j}-th element
     * is either {@code null} or an {@code int[]} whose elements are an
     * increasing list of the haplotypes that carry the {@code j}-th allele.
     *
     * @return an array of length {@code marker.nAlleles()} with
     * a unique {@code null} element and whose {@code j}-th element
     * is either {@code null} or an {@code int[]} whose elements are an
     * increasing list of the haplotypes that carry the {@code j}-th allele
     */
    int[][] alleleToHaps();

    /**
     * Returns the number of haplotypes carrying an allele {@code j} where
     * {@code (j == this.nullRow())}.
     *
     * @return the number of haplotypes carrying an allele {@code j} where
     * {@code (j == this.nullRow())}
     */
    int nullRowAlleleCnt();

    /**
     * Returns an array of length {@code this.nAlleles()} whose {@code j}-th
     * element is the allele count of the {@code j}-th allele.
     * @return an array of allele counts
     */
    int[] alleleCounts();

    /**
     * Returns the number of haplotypes that carry the specified allele.
     * @param allele an allele index
     * @return the number of haplotypes that carry the specified allele
     * @throws IndexOutOfBoundsException if
     * {@code ( (allele < 0) ||  (allele >= this.nAlleles()) )}
     */
    int alleleCount(int allele);

    /**
     * Returns an {@code IndexArray} with {@code this.size()} elements that maps
     * haplotype to allele.
     *
     * @return an {@code IndexArray} with {@code this.size()} elements that maps
     * haplotype to allele
     */
    IndexArray hapToAllele();

    /**
     * Returns {@code true}.
     * @param sample the sample index
     * @return {@code true}
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    @Override
    boolean isPhased(int sample);

    /**
     * Returns {@code true}.
     * @return {@code true}
     */
    @Override
    boolean isPhased();

    /**
     * Returns index of the haplotype that carries the specified copy of the
     * specified allele.
     * @param allele an allele index
     * @param copy a copy index
     * @return index of the haplotype that carries the specified allele
     * @throws IllegalArgumentException if
     * {@code (allele == this.nullRow())}
     * @throws IndexOutOfBoundsException if
     * {@code ((allele < 0) || (allele >= this.nAlleles()))}
     * @throws IndexOutOfBoundsException if
     * {@code ((copy < 0) || (copy >= this.alleleCount(allele)))}
     */
    int nonNullRowHap(int allele, int copy);

    /**
     * Returns {@code true} if the specified haplotype carries the specified
     * allele and return {@code false} otherwise.
     * @param allele an allele index
     * @param hap a haplotype index
     * @return {@code true} if the specified haplotype carries the specified
     * allele
     * @throws IndexOutOfBoundsException if
     * {@code ((hap < 0) || (hap >= this.size()))}
     * @throws IndexOutOfBoundsException if
     * {@code ((allele < 0) || (allele >= this.nAlleles()))}
     */
    boolean isCarrier(int allele, int hap);

    /**
     * Returns {@code this.maps().length}
     * @return this.maps().length
     */
    int nMaps();

    /**
     * Returns an array of maps, whose composition maps haplotype indices
     * to alleles.  The allele on haplotype {@code h} is determined
     * by the following calculation:
     * <pre>
            IntArray[] maps = this.maps();
            int value = maps[0].get(h);
            for (int j=1; j&lt;maps.length; ++j) {
               value = indexArrays[j].get(value);
            }
            int allele = value
       </pre>
     * @return an array of maps, whose composition maps haplotype indices
     * to alleles
     */
    IntArray[] maps();

    /**
     * Returns {@code this.maps()[index]}.
     * @param index the index in {@code this.maps()}
     * @return {@code this.maps()[index]}
     * @throws IndexOutOfBoundsException if
     * {@code ((index < 0) || (index >= this.nMaps()))}
     */
    IntArray map(int index);

    /**
     * Returns the data in this {@code RefGTRec} as a string VCF record with
     * correct INFO/AN and INFO/AC fields and with FORMAT/GT  as the only
     * FORMAT field.
     * @return the data represented by {@code this} as a string VCF record
     */
    @Override
    public String toString();

    /**
     * Returns the data in this {@code RefGTRec} as a string VCF record with
     * correct INFO/AN and INFO/AC fields and with FORMAT/GT  as the only
     * FORMAT field.
     * @return the data represented by {@code this} as a string VCF record
     */
    public String toVcfRecord();

    /**
     * Returns the data in the specified {@code RefGTRec} object as a
     * string VCF record with with FORMAT/GT as the only FORMAT field.
     * If {@code (isHaploid == null)}, the returned string VCF record will
     * have correct INFO/AN and INFO/AC fields.  Otherwise, no changes
     * will be made to the VCF record INFO field.
     * If {@code (isHaploid != null)} and
     * {@code ( (0 <= j) && (j < isHaploid.size()) && (isHaploid.get(j) == true) )}
     * the second haplotype of the {@code j}-th sample will not printed.
     * @param isHaploid a boolean array with {@code this.samples().size()}
     * elements whose {@code j}-th value is {@code true} if the second
     * haplotype of the {@code j}-th sample should not be printed
     * @return the data represented by {@code this} as a string VCF record
     *
     * @throws IllegalArgumentException if {
     * {@code ((isHaploid !=null) && (isHaploid.size() != this.samples().size())}
     */
    public String toVcfRecord(BooleanArray isHaploid);

    /**
     * Returns the data in the specified {@code RefGTRec} object as a
     * string VCF record with correct INFO/AN and INFO/AC fields and
     * with FORMAT/GT  as the only FORMAT field. The implementation of
     * this method has suboptimal computational efficiency if
     * {@code (rec.isAlleleCoded() == true)}.
     * @param rec a string VCF record with phased, non-missing genotypes
     * @return the data represented by {@code rec} as a string VCF record
     * @throws NullPointerException if {@code (rec == null)}
     */
    public static String toString(RefGTRec rec) {
        int nHaps = rec.size();
        StringBuilder sb = new StringBuilder();
        MarkerUtils.appendFirst8Fields(rec.marker(), sb); // no INFO/{AN,AC} update
//        MarkerUtils.appendFirst8Fields(rec.marker(), nHaps, rec.alleleCounts(), sb);
        sb.append("\tGT");
        assert (nHaps & 0b1) == 0;
        for (int h=0; h<nHaps; h+=2) {
            sb.append(Const.tab);
            sb.append(rec.get(h));
            sb.append(Const.phasedSep);
            sb.append(rec.get(h | 0b1));
        }
        return sb.toString();
    }

    /**
     * Returns the data in the specified {@code RefGTRec} object as a
     * string VCF record with FORMAT/GT  as the only FORMAT field.
     * The implementation of this method has suboptimal computational
     * efficiency if {@code (rec.isAlleleCoded() == true)}. If
     * {@code (0 <= j) && (j < isHaploid.size()) && (isHaploid.get(j) == true)}
     * the second haplotype of the {@code j}-th sample will not printed.
     * @param rec a VCF record with phased, non-missing genotypes
     * @param isHaploid a boolean array with {@code this.samples().size()}
     * elements whose {@code j}-th value is {@code true} if the second
     * haplotype of the {@code j}-th sample should not be printed
     * @return the data represented by {@code rec} as a string VCF record
     *
     * @throws IllegalArgumentException if
     * {@code ((isHaploid.size() != rec.samples().size())}
     * @throws NullPointerException if
     * {@code ((rec == null) || (isHaploid == null))}
     */
    public static String toString(RefGTRec rec, BooleanArray isHaploid) {
        if (isHaploid.size()!=rec.samples().size()) {
            throw new IllegalArgumentException(String.valueOf(isHaploid.size()));
        }
        int nHaps = rec.size();
        StringBuilder sb = new StringBuilder();
        MarkerUtils.appendFirst8Fields(rec.marker(), sb); // no INFO/{AN,AC} update
        sb.append("\tGT");
        assert (nHaps & 0b1) == 0;
        for (int h=0; h<nHaps; h+=2) {
            sb.append(Const.tab);
            sb.append(rec.get(h));
            if (isHaploid.get(h>>1)==false) {
                sb.append(Const.phasedSep);
                sb.append(rec.get(h | 0b1));
            }
        }
        return sb.toString();
    }
}

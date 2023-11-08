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

import ints.IndexArray;
import ints.IntArray;

/**
 * <p>Interface {@code RefGTRec} represents represents phased genotype data
 * for one marker.  For implementations of this interface, unless otherwise
 * specified in the implementation documentation, if the {@code isAlleleCoded()}
 * method returns {@code false}, the {@code majorAllele()},
 * {@code alleleCount()}, and {@code hapIndex()} methods will be computationally
 * expensive with compute time proportional to the number of haplotypes.
 * Alternatively if the {@code isAlleleCoded()} method returns
 * {@code true}, the {@code maps()} and {@code map()} methods will be
 * computationally expensive with compute time proportional to the number
 * of haplotypes.
 * </p>
 * <p>All instances of {@code RefGTRec} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface RefGTRec extends GTRec {

    /**
     * Returns an allele-coded {@code RefGTRec} instance for the
     * specified data.
     * @param rec the phased, non-missing genotype data
     * @return an allele-coded {@code RefGTRec} instance for the
     * specified data
     * @throws NullPointerException if {@code rec == null}
     */
    static RefGTRec alleleRefGTRec(RefGTRec rec) {
        if (rec.isAlleleCoded()) {
            return rec;
        }
        if (rec.marker().nAlleles()==2) {
            return new TwoAlleleRefGTRec(rec);
        }
        else {
            return new AlleleRefGTRec(rec);
        }
    }

    /**
     * Constructs and returns a new allele-coded {@code RefGTRec} instance
     * from the specified data.
     *
     * @param gtp a VCF record parser that extracts sample genotypes
     * @return an allele-coded {@code RefGTRec} instance
     *
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased genotype or missing allele
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws NullPointerException if {@code gtp == null}
     */
    static RefGTRec alleleRefGTRec(VcfRecGTParser gtp) {
        if (gtp.nAlleles()==2) {
            return new TwoAlleleRefGTRec(gtp);
        }
        else {
            return new AlleleRefGTRec(gtp);
        }
    }

    /**
     * Constructs and returns a new allele-coded {@code RefGTRec} instance
     * from the specified data.
     *
     * @param marker the marker
     * @param samples the samples
     * @param alleleToHaps an array whose {@code j}-th element is {@code null}
     * if {@code j} is the major allele and otherwise is a list of haplotypes
     * sorted in increasing order that carry the {@code j}-th allele.
     * If there is more than one allele with a maximal allele count, the
     * major allele is the smallest allele with maximal allele count.
     * If a haplotype is contained in a list for more than one non-major allele,
     * the haplotype will be assumed to carry the smallest allele.
     * @return an allele-coded {@code RefGTRec} instance
     *
     * @throws IllegalArgumentException if the {@code (hapIndices[j] == null)}
     * and {@code j} is not the major allele or if
     * {@code (hapIndices[j] != null)} and {@code j} is the major allele
     * @throws IllegalArgumentException if any non-null element of
     * {@code hapIndices} is not a sorted list of distinct haplotype indices
     * between 0 (inclusive) and {@code 2*samples.size()} (exclusive)
     * @throws IllegalArgumentException if
     * {@code marker.nAlleles() != hapIndices.length}
     * @throws NullPointerException if
     * {@code marker == null || samples == null || hapIndices == null}
     */
    static RefGTRec alleleRefGTRec(Marker marker, Samples samples,
            int[][] alleleToHaps) {
        if (marker.nAlleles()==2) {
            return new TwoAlleleRefGTRec(marker, samples, alleleToHaps);
        }
        else {
            return new AlleleRefGTRec(marker, samples, alleleToHaps);
        }
    }

    /**
     * Returns an array of length {@code this.marker().nAlleles()} whose
     * {@code j}-th element is {@code null} or is a list of the
     * haplotypes in increasing order that carry allele {@code j}. Exactly
     * one element of the returned array must be {@code null}. It is
     * recommended that the {@code null} element correspond to the major allele
     * with lowest index, but this is not a requirement for implementations
     * of this method.
     *
     * @return Returns an array of length {@code this..marker().nAlleles()} whose
     * {@code j}-th element is {@code null} or is a list of the
     * haplotypes in increasing order that carry allele {@code j}
     */
    int[][] alleleToHaps();

    /**
     * Returns an {@code IndexArray} with {@code this.size()} elements that maps
     * haplotype to allele.
     *
     * @return an {@code IndexArray} with {@code this.size()} elements that maps
     * haplotype to allele
     */
    IndexArray hapToAllele();

    /**
     * Returns the sum of the lengths of non-null rows of
     * {@code this.alleleToHaps()}.
     *
     * @return the sum of the lengths of non-null rows of
     * {@code this.alleleToHaps()}
     */
    int nAlleleCodedHaps();

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
     * Returns {@code true} if this instance stores the indices of haplotypes
     * that carry non-major alleles, and returns {@code false} otherwise.
     *
     * @return {@code true} if this instance stores the indices of haplotypes
     * that carry non-major alleles
     */
    boolean isAlleleCoded();

    /**
     * Returns the major allele with lowest index.
     *
     * @return the major allele with lowest index
     */
    int majorAllele();

    /**
     * Returns an array of length {@code this.nAlleles()} whose {@code j}-th
     * element is the allele count of the {@code j}-th allele.
     * @return an array of allele counts
     */
    int[] alleleCounts();

    /**
     * Returns the number of haplotypes that carry the specified
     * non-major allele.
     * @param allele an allele index
     * @return the number of haplotypes that carry the specified non-major
     * allele
     * @throws IllegalArgumentException if
     * {@code allele == this.majorAllele()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 ||  allele >= this.nAlleles()}
     */
    int alleleCount(int allele);

    /**
     * Returns index of the haplotype that carries the specified copy of the
     * specified allele.
     * @param allele an allele index
     * @param copy a copy index
     * @return index of the haplotype that carries the specified allele
     * @throws IllegalArgumentException if
     * {@code allele == this.majorAllele()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 ||  allele >= this.nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code copy < 0 ||  copy >= this.alleleCount(allele)}
     */
    int hapIndex(int allele, int copy);

    /**
     * Returns {@code true} if the specified haplotype carries the specified
     * allele and return {@code false} otherwise.
     * @param allele an allele index
     * @param hap a haplotype index
     * @return {@code true} if the specified haplotype carries the specified
     * allele
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.size()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.nAlleles()}
     */
    boolean isCarrier(int allele, int hap);

    /**
     * Returns {@code this.maps().length}
     * @return this.maps().length
     */
    int nMaps();

    /**
     * Returns an array of maps, which when composed map haplotype indices
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
     * @return an array of maps, which when composed map haplotype indices
     * to alleles
     */
    IntArray[] maps();

    /**
     * Returns {@code this.maps()[index]}.
     * @param index the index in {@code this.maps()}
     * @return {@code this.maps()[index]}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMaps()}
     */
    IntArray map(int index);
}
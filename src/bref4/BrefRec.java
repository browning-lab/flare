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
import vcf.Marker;
import vcf.RefGTRec;
import vcf.VcfRecGTParser;

/**
 * <p>Interface {@code BrefRec} stores a marker and the list of sequence
 * indices that carry each marker allele.</p>
 *
 * <p>Classes that implement {@code BrefRec} are required to be immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface BrefRec extends IntArray {

    /**
     * Returns a new {@code BrefRec} instance constructed from the specified
     * data.
     * @param rec a VCF record with phased, non-missing alleles.
     * @return a new {@code BrefRec} instance constructed from the specified
     * data
     * @throws NullPointerException if {@code rec == null}
     */
    public static BrefRec from(RefGTRec rec) {
        if (rec.marker().nAlleles()==2) {
            return new DiallelicBrefRec(rec);
        }
        else {
            return new AllelicBrefRec(rec);
        }
    }

    /**
     * Returns a new {@code BrefRec} instance constructed from the specified
     * data.
     * @param gtp a VCF record with phased, non-missing alleles.
     * @return a new {@code BrefRec} instance constructed from the specified
     * data
     * @throws NullPointerException if {@code gtp == null}
     */
    public static BrefRec from(VcfRecGTParser gtp) {
        if (gtp.marker().nAlleles()==2) {
            return new DiallelicBrefRec(gtp);
        }
        else {
            return new AllelicBrefRec(gtp);
        }
    }

    /**
     * Returns the sum of the lengths of non-null rows of the specified
     * two-dimensional array.
     * @param alleleToHaps a two-dimensional array
     * @return the sum of the lengths of non-null rows
     * @throws NullPointerException if {@code alleleToHaps == null}
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
     * {@code this.allelesToHap()[j] == null} or {@code -1}
     * if no such allele index exists.
     *
     * @param alleleToHaps a two-dimensional array
     * @return the smallest row index {@code j} such that
     * {@code this.allelesToHap[j] == null} or {@code -1}
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
     * Returns the estimated number of bytes consumed by this object,
     * excluding the overhead bytes required by {@code this}.
     * @return the estimated number of bytes required to store this object
     */
    long estBytes();

    /**
     * Returns the marker.
     * @return the marker
     */
    Marker marker();

    /**
     * Returns the number of haplotypes.
     * @return the number of haplotypes
     */
    @Override
    int size();

    /**
     * Returns an array of length {@code this.marker().nAlleles()}
     * whose {@code j}-th element is either {@code null} or is an increasing
     * list of indices of haplotypes that carry allele {@code j}. Exactly
     * one element of the returned array will be {@code null}.
     *
     * @return an array of length {@code this.marker().nAlleles()}
     * whose {@code j}-th element is either {@code null} or is an increasing
     * list of indices of haplotypes that carry allele {@code j}
     */
    int[][] alleleToHaps();

    /**
     * Returns an {@code IndexArray} with {@code this.size()} elements that maps
     * haplotype to allele.
     *
     * @return an {@code IndexArray} with {@code this.size()} elements that maps
     * sequence to allele
     */
    IndexArray hapToAllele();

    /**
     * Returns the index of the {@code null} row of {@code this.alleleToHaps()}.
     *
     * @return the index of the {@code null} row of {@code this.alleleToHaps()}
     */
    int nullRow();

    /**
     * Returns the sum of the lengths of non-null rows of
     * {@code this.alleleToHaps()}.
     *
     * @return the sum of the lengths of non-null rows of
     * {@code this.alleleToHaps()}
     */
    int nStoredHapIndices();

    /**
     * Returns an {@code BrefRec} that is obtained from {@code this} by
     * applying the specified map to the haplotype indices.  The
     * contract for this method is unspecified if there exists
     * haplotypes {@code h1} and {@code h2} such that
     * {@code (0 < h1 && h1 < h2 && h2 < this.size())} and
     * {@code (this.get(h1) != this.get(h2))} and
     * {@code (map.get(h1) == map.get(h2))}.
     *
     * @param map a map from haplotype to sequence index
     * @return an {@code BrefRec} that is obtained from {@code this} by
     * applying the specified map to the haplotype indices
     * @throws IllegalArgumentException if
     * {@code this.size() != map.size()}
     * @throws NullPointerException if {@code map == null}
     */
    BrefRec applyMap(IndexArray map);
}

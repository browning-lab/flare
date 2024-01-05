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

/**
 * <p>Interface {@code GT} represents genotype data
 * for a list of markers and a list of samples.
 * </p>
 * <p>All instances of {@code GT} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface GT {

    /**
     * Returns the estimated number of bytes consumed by this object,
     * excluding the overhead bytes required by {@code this}.
     * @return the estimated number of bytes required to store this object
     */
    long estBytes();

    /**
     * Returns {@code true} if the markers are ordered by decreasing chromosome
     * base position, and returns {@code false} otherwise.
     * @return {@code true} if the markers are ordered by decreasing chromosome
     * base position
     */
    boolean isReversed();

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    Marker marker(int marker);

    /**
     * Returns the list of markers in order of increasing chromosome position.
     * If {@code (this.isReversed() == false)} then
     * {@code (this.markers().marker(j).equals(this.marker(j)) == true)} for
     * all {@code (0 <= j && j < this.nMarkers())}.
     * If {@code (this.isReversed() == true)} then
     * {@code (this.markers().marker(this.nMarkers() - 1 - j).equals(this.marker(j)) == true)}
     * for all {@code (0 <= j && j < this.nMarkers())}
     * @return the list of markers in order of increasing chromosome position
     */
    Markers markers();

    /**
     * Returns the number of haplotypes.  The returned value is equal to
     * {@code 2*this.nSamples()}.
     * @return the number of haplotypes
     */
    int nHaps();

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    int nSamples();

   /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();

    /**
     * Returns {@code true} if the genotype for each marker and sample
     * has non-missing alleles and is either haploid or diploid with
     * a phased allele separator, and returns {@code false} otherwise.
     * @return {@code true} if the genotype for each marker and sample
     * is a phased, non-missing genotype
     */
    boolean isPhased();

    /**
     * Returns the allele on the specified haplotype for the specified marker
     * or return -1 if the allele is missing. The order of the two alleles
     * is unspecified if {@code this.isPhased() == false}.
     * @param marker the marker index
     * @param hap the haplotype index
     * @return the allele on the specified haplotype for the specified marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap  >= this.nHaps()}
     */
    int allele(int marker, int hap);

    /**
     * Returns a {@code GT} instance restricted to genotype data for
     * the specified markers.
     * @param markers the list of markers in the returned instance
     * @param indices a list of distinct marker indices (from
     * {@code this.markers())} in increasing order
     * @return a {@code GT} instance restricted to genotype data for
     * the specified markers
     *
     * @throws IndexOutOfBoundsException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (indices[j] < 0 || indices[j] >= this.nMarkers())}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (1 <= j && j < indices.length)} such that
     * {@code (indices[j] <= indice[j - 1])}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (this.marker(indices[j]).equals(markers.marker(j)) == false)}
     * @throws NullPointerException if {@code indices == null}
     * @throws UnsupportedOperationException if {@code this.isReversed() == true}
     */
    GT restrict(Markers markers, int[] indices);

    /**
     * Returns a new {@code GT} instance restricted to genotype data for
     * the specified markers.
     * @param start the start marker (inclusive)
     * @param end the end marker (exclusive)
     * @return a {@code GT} instance restricted to genotype data for
     * the specified markers
     * @throws IllegalArgumentException if {@code start >= end}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 ||  end > this.markers()}
     */
    GT restrict(int start, int end);
}

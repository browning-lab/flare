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

import beagleutil.ChromIds;
import blbutil.BitArray;
import blbutil.Const;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * <p>Class {@code Markers} represent a list of markers in chromosome order.
 * </p>
 * <p>Instances of class {@code Markers} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Markers {

    private final Marker[] markers;
    private final Set<Marker> markerSet;
    private final int[] sumAlleles;
    private final int[] sumHapBits;
    private final int hashCode;

    /**
     * Returns the estimated number of bytes consumed by this object,
     * excluding the overhead bytes required by {@code this}.
     * @return the estimated number of bytes required to store this object
     */
    public long estBytes() {
        int overhead = (5 + markers.length)*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + (5 + 2*markers.length)*8;  // assume 8 bytes per reference
        estBytes += markers.length*104; // assume diallelic SNV with shared String[] alleles and no String[] ids
        estBytes += (1 + sumAlleles.length)*4; // 4 bytes per array to store length;
//        estBytes += (1 + sumGenotypes.length)*4; // 4 bytes per array to store length;
        estBytes += (1 + sumHapBits.length)*4; // 4 bytes per array to store length;
        estBytes += 4;
        return estBytes;
    }

    /**
     * Returns a new {@code Markers} instance that is constructed from
     * the specified data.
     * @param markers a list of markers in chromosome order
     * @return a new {@code Markers} instance corresponding to the
     * specified list of markers
     *
     * @throws IllegalArgumentException if markers on a chromosome are not
     * in chromosome order
     * @throws IllegalArgumentException if there are duplicate markers
     * @throws IllegalArgumentException if the markers on a chromosome
     * do not form a contiguous set of entries within the array
     *
     * @throws NullPointerException if
     * {@code markers == null} or if {@code markers[j] == null}
     * for any {@code j} satisfying {@code (0 <= j && j < markers.length)}
     */
    public static Markers create(Marker[] markers) {
        return new Markers(markers);
    }

    /**
     * Construct a new {@code Markers} instance that represents the
     * specified list of markers.
     * @param markers a list of markers in chromosome order
     *
     * @throws IllegalArgumentException if markers on a chromosome are not
     * in chromosome order
     * @throws IllegalArgumentException if there are duplicate markers
     * @throws IllegalArgumentException if the markers on a chromosome
     * do not form a contiguous set of entries within the array
     *
     * @throws NullPointerException if
     * {@code markers == null} or if {@code markers[j] == null}
     * for any {@code j} satisfying {@code (0 <= j && j < markers.length)}
     */
    private Markers(Marker[] markers) {
        checkMarkerPosOrder(markers);
        this.markers = markers.clone();
        this.markerSet = markerSet(markers);

        this.sumAlleles = cumSumAlleles(markers);
        this.sumHapBits = cumSumHaplotypeBits(markers);
        this.hashCode = Arrays.deepHashCode(markers);
    }

    private static void checkMarkerPosOrder(Marker[] markers) {
        if (markers.length < 2) {
            return;
        }
        Set<Integer> chromIndices = new HashSet<>();
        chromIndices.add(markers[0].chromIndex());
        chromIndices.add(markers[1].chromIndex());
        for (int j=2; j<markers.length; ++j) {
            int chr0 = markers[j-2].chromIndex();
            int chr1 = markers[j-1].chromIndex();
            int chr2 = markers[j].chromIndex();
            if (chr0 == chr1 && chr1==chr2) {
                int pos0 = markers[j-2].pos();
                int pos1 = markers[j-1].pos();
                int pos2 = markers[j].pos();
                if ( (pos1<pos0 && pos1<pos2) || (pos1>pos0 && pos1>pos2) ) {
                    String s = "markers not in chromosomal order: "
                            + Const.nl + markers[j-2]
                            + Const.nl + markers[j-1]
                            + Const.nl + markers[j];
                    throw new IllegalArgumentException(s);
                }
            }
            else if (chr1!=chr2) {
                if (chromIndices.contains(chr2)) {
                    String s = "markers on chromosome are not contiguous: "
                            + ChromIds.instance().id(chr2);
                    throw new IllegalArgumentException(s);
                }
                chromIndices.add(chr2);
            }
        }
    }

    private static Set<Marker> markerSet(Marker[] markers) {
        Set<Marker> markerSet = new HashSet<>(markers.length);
        for (Marker m : markers) {
            if (markerSet.add(m)==false) {
                throw new IllegalArgumentException("Duplicate marker: " + m);
            }
        }
        return markerSet;
    }

    private static int[] cumSumAlleles(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            ia[j] = ia[j-1] + markers[j-1].nAlleles();
        }
        return ia;
    }

    /**
     * Return an array of length {@code this.size() + 1} whose
     * {@code k}-th value is the the sum of the number of possible genotypes
     * for markers with index less than {@code k}
     * @return  an array whose {@code k}-th value is the the sum of the number
     * of possible genotypes for markers with index less than {@code k}
     */
    public int[] cumSumGenotypes() {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            ia[j] = ia[j-1] + MarkerUtils.nGenotypes(markers[j-1].nAlleles());
        }
        return ia;
    }

    private static int[] cumSumHaplotypeBits(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            int nAllelesM1 = markers[j-1].nAlleles() - 1;
            int nStorageBits = Integer.SIZE
                    - Integer.numberOfLeadingZeros(nAllelesM1);
            ia[j] = ia[j-1] + nStorageBits;
        }
        return ia;
    }

    /**
     * Returns a hash code value for the object.
     * The returned hash code equals
     * {@code Arrays.deepHashCode(this.markers())}.
     * @return a hash code value for the object
     */
    @Override
    public int hashCode() {
        return hashCode;
    }

    /**
     * Returns {@code true} if the specified object is a {@code Markers}
     * instance which represents the same list of markers as {@code this},
     * and returns {@code false} otherwise. Two lists of markers are
     * the same if the lists have the same size and if markers with the
     * same index in the two lists are equal.
     *
     * @param obj the object to be tested for equality with {@code this}
     *
     * @return {@code true} if the specified object is a {@code Markers}
     * instance which represents the same list of markers as {@code this}
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Markers other = (Markers) obj;
        return Arrays.deepEquals(this.markers, other.markers);
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int size() {
        return markers.length;
    }

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public Marker marker(int marker) {
        return markers[marker];
    }

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    public Marker[] markers() {
        return markers.clone();
    }

    /**
     * Returns {@code true} if the specified marker is not {@code null}
     * and is an element in the list of markers represented by {@code this},
     * and returns {@code false} otherwise.
     *
     * @param marker a marker
     *
     * @return {@code true} if the specified marker is not {@code null} and
     * is an element in the list of markers represented by {@code this}
     */
    public boolean contains(Marker marker) {
        return markerSet.contains(marker);
    }

    /**
     * Returns a {@code Markers} instance that represents
     * the specified range of marker indices.
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @return a {@code Markers} instance that represents
     * the specified range of marker indices
     *
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end > this.nMarkers()}
     * @throws IllegalArgumentException if {@code start >= end}
     */
    public Markers restrict(int start, int end) {
        if (end > markers.length) {
            throw new IndexOutOfBoundsException("end > this.nMarkers(): " + end);
        }
        return new Markers(Arrays.copyOfRange(markers, start, end));
    }

    /**
     * Returns a {@code Markers} instance that represents
     * the specified markers.
     * @param indices a list of distinct marker indices in increasing order
     * @return a new {@code Markers} instance that represents the specified
     * markers
     *
     * @throws IndexOutOfBoundsException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (indices[j] < 0 || indices[j] >= this.nMarkers())}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (1 <= j && j < indices.length)} such that
     * {@code (indices[j] <= indice[j - 1])}
     * @throws NullPointerException if {@code indices == null}
     */
    public Markers restrict(int[] indices) {
        Marker[] ma = new Marker[indices.length];
        ma[0] = markers[indices[0]];
        for (int j=1; j<indices.length; ++j) {
            if (indices[j] <= indices[j-1]) {
                throw new IllegalArgumentException(String.valueOf(indices[j]));
            }
            ma[j] = markers[indices[j]];
        }
        return new Markers(ma);
    }

    /**
     * Returns the sum of the number of alleles for
     * the markers with index less than the specified index.
     * @param marker a marker index
     * @return the sum of the number of alleles for
     * the markers with index less than the specified index
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker > this.nMarkers()}
     */
    public int sumAlleles(int marker) {
        return sumAlleles[marker];
    }

    /**
     * Returns {@code this.sumAlleles(this.nMarkers())}.
     * @return {@code this.sumAlleles(this.nMarkers())}
     */
    public int sumAlleles() {
        return sumAlleles[markers.length];
    }

    /**
     * Returns the number of bits requires to store a haplotype for the
     * markers with index less than the specified index.
     * @param marker a marker index
     * @return the number of bits requires to store a haplotype for the
     * markers with index less than the specified index
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker > this.nMarkers()}
     */
    public int sumHapBits(int marker) {
        return sumHapBits[marker];
    }

    /**
     * Returns {@code this.sumHaplotypeBits(this.nMarkers())}.
     * @return {@code this.sumHaplotypeBits(this.nMarkers())}
     */
    public int sumHapBits() {
        return sumHapBits[markers.length];
    }

    /**
     * Returns the specified allele stored in the specified {@code hapBits}
     * array.  The contract for this method is undefined if the specified
     * {@code hapBits} array was not created with the
     * {@code this.allelesToBits()} method.
     * @param hapBits the bit array storing the haplotype alleles
     * @return the specified allele stored in the specified {@code hapBits}
     * array.
     * @throws NullPointerException if {@code hapBits == null}
     */
    public int[] bitsToAlleles(BitArray hapBits) {
        int[] alleles = new int[markers.length];
        for (int m=0; m<alleles.length; ++m) {
            int start = sumHapBits[m];
            int end = sumHapBits[m+1];
            if (end==(start+1)) {
                alleles[m] = hapBits.get(start) ? 1 : 0;
            }
            else {
                int allele = 0;
                int mask = 1;
                for (int j=start; j<end; ++j) {
                    if (hapBits.get(j)) {
                        allele |= mask;
                    }
                    mask <<= 1;
                }
                alleles[m] = allele;
            }
        }
        return alleles;
    }

    /**
     * Stores the specified alleles in the specified {@code bitList}
     * @param alleles a sequence of alleles
     * @param bitList a sequence of bits
     * @throws IllegalArgumentException if
     * {@code alleles.length != this.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code bitList.size() != this.sumHaplotypeBits()}
     * @throws IllegalArgumentException if there exists a {@code k}
     * such that {@code (0 < k && k < alleles.length)} and
     * {@code (alleles[k] < 0 || alleles[k] >= this.marker(k).nAlleles()})
     * @throws NullPointerException if
     * {@code alleles == null || bitList == null}
     */
    public void allelesToBits(int[] alleles, BitArray bitList) {
        if (alleles.length != markers.length) {
            throw new IllegalArgumentException(String.valueOf(alleles.length));
        }
        if (bitList.size() != this.sumHapBits()) {
            throw new IllegalArgumentException(String.valueOf(bitList.size()));
        }
        for (int k=0; k<alleles.length; ++k) {
            int allele = alleles[k];
            if (allele < 0 || allele >= markers[k].nAlleles()) {
                String s = "allele \"" + allele + "\" out of bounds for marker: "
                        + markers[k];
                throw new IllegalArgumentException(s);
            }
            int mask = 1;
            for (int j=sumHapBits[k]; j<sumHapBits[k+1]; ++j) {
                if ((allele & mask)==mask) {
                    bitList.set(j);
                }
                else {
                    bitList.clear(j);
                }
                mask <<= 1;
            }
        }
    }

    /**
     * Stores the specified alleles in the specified {@code bitList}
     * @param marker a marker index
     * @param allele an allele index
     * @param bitList a sequence of bits
     * @throws IllegalArgumentException if
     * {@code alleles.length != this.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code bitList.size() != this.sumHaplotypeBits()}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker(marker).nAlleles()}
     * @throws NullPointerException if {@code bitList == null}
     */
    public void setAllele(int marker, int allele, BitArray bitList) {
        if (bitList.size() != this.sumHapBits()) {
            throw new IllegalArgumentException(String.valueOf(bitList.size()));
        }
        if (marker < 0 || marker >= markers.length) {
            throw new IndexOutOfBoundsException(String.valueOf(marker));
        }
        if (allele < 0 || allele >= markers[marker].nAlleles()) {
            throw new IndexOutOfBoundsException(String.valueOf(allele));
        }
        int mask = 1;
        for (int j=sumHapBits[marker], n=sumHapBits[marker+1]; j<n; ++j) {
            if ((allele & mask)==mask) {
                bitList.set(j);
            }
            else {
                bitList.clear(j);
            }
            mask <<= 1;
        }
    }

    /**
     * Returns the specified allele
     * @param hapBits a haplotype encoded as bits with the
     * {@code this.allelesToBits() method}
     * @param marker a marker index
     * @return the specified allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hapBits.size() < this.sumHaplotypeBits(marker + 1)}
     * @throws NullPointerException if {@code hapBits == null}
     */
    public int allele(BitArray hapBits, int marker) {
        int start = sumHapBits[marker];
        int end = sumHapBits[marker+1];
        if (end==(start+1)) {
            return hapBits.get(start) ? 1 : 0;
        }
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            if (hapBits.get(j)) {
                allele |= mask;
            }
            mask <<= 1;
        }
        return allele;
    }

    /**
     * Returns a string representation of {@code this}.
     * The exact details of the representation are unspecified and
     * subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return Arrays.toString(markers);
    }
}

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

import beagleutil.ChromIds;
import blbutil.Const;
import ints.IntList;

/**
 * <p>Class {@code Marker} represents a VCF record's CHROM, POS, ID, REF,
 * ALT, QUAL, FILTER, and INFO fields.</p>
 *
 * <p>Instances of class {@code Marker} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Marker implements Comparable<Marker> {

    /**
     * The maximum number of distinct chromosome identifiers indexed by
     * the {@code ChromIds.instance().getIndex()} method.
     */
    public static final int MAX_CHROMOSOMES = 0x10000;

    /**
     * The maximum number of alleles in a VCF record's REF and ALT fields
     */
    public static final int MAX_ALLELES = 0x10000;

    private final short chromIndex;
    private final int pos;
    private final short nAlleles;
    private final String fields;

    /**
     * Constructs a {@code Marker} instance for the specified data.
     * The contract for constructed {@code Marker} is undefined
     * if {@code (chromIndex < 0)}, if {@code (nAlleles < 0)},
     * if {@code (fields == null)}, if the tab-delimited VCF ID, REf, ALT,
     * QUAL, FILTER, and INFO in {@code fields} do not conform
     * to the VCF specification, or if {@code nAlleles} does not equal
     * the number of alleles in the REF and ALT fields of {@code fields}.
     * @param chromIndex the chromosome index
     * @param pos the base coordinate
     * @param nAlleles the number of alleles
     * @param fields the tab-delimited VCF ID, REF, ALT, QUAL, FILTER, and
     * INFO fields
     * @throws IndexOutOfBoundsException if
     * {@code (chromIndex >= Marker.MAX_CHROMOSOMES)}
     * @throws IndexOutOfBoundsException if the number of alleles in the REF
     * and ALT fields of {@code fields} is greater than
     * {@code Marker.MAX_ALLELES}
     */
    public Marker(int chromIndex, int pos, int nAlleles, String fields) {
        if (chromIndex > MAX_CHROMOSOMES) {
            throw new IllegalArgumentException(String.valueOf(chromIndex));
        }
        if (nAlleles > MAX_ALLELES) {
            throw new IllegalArgumentException(String.valueOf(nAlleles));
        }
        this.chromIndex = (short) chromIndex;
        this.pos = pos;
        this.nAlleles = (short) nAlleles;
        this.fields = fields;
    }

    /**
     * Constructs a new {@code Marker} instance from the data from the first
     * eight tab-delimited fields of the specified VCF record.  The contract
     * for this method is undefined if the first 8 fields of the VCF record
     * do not conform to the VCF 4.5 specification.
     * @param vcfRec A string that begins with the first 8 tab-delimited fields
     * of a VCF record
     * @return a new {@code Marker} instance
     *
     * @throws IndexOutOfBoundsException if the VCF record CHROM field,
     * {@code chrom} satisfies
     * {@code ChromIds.instance().getIndex(chrom) >= Marker.MAX_CHROMOSOMES}
     * @throws IndexOutOfBoundsException if the number of alleles in the REF
     * and ALT fields of the VCF records is greater than
     * {@code Marker.MAX_ALLELES}
     * @throws NullPointerException if {@code (vcfRecord == null)}
     */
    public static Marker fromVcfRecord(String vcfRec) {
        return fromVcfRecord(vcfRec, false);
    }

    /**
     * Constructs a new {@code Marker} instance from the data from the first
     * eight tab-delimited fields of the specified VCF record.  The contract
     * for this method is undefined if the first 8 fields of the VCF record
     * do not conform to the VCF 4.5 specification.
     * @param vcfRec A string that begins with the first 8 tab-delimited fields
     * of a VCF record
     * @param stripOptionalFields {@code true} if the VCF record's ID,
     * QUAL, FILTER, and INFO subfields should be discarded
     * @return a new {@code Marker} instance
     *
     * @throws IndexOutOfBoundsException if the VCF record CHROM field,
     * {@code chrom} satisfies
     * {@code ChromIds.instance().getIndex(chrom) >= Marker.MAX_CHROMOSOMES}
     * @throws IndexOutOfBoundsException if the number of alleles in the REF
     * and ALT fields of the VCF records is greater than
     * {@code Marker.MAX_ALLELES}
     * @throws NullPointerException if {@code (vcfRec == null)}
     */
    public static Marker fromVcfRecord(String vcfRec, boolean stripOptionalFields) {
        IntList tabs = MarkerUtils.first8TabIndices(vcfRec);
        StringBuilder sb = new StringBuilder();
        short chromIndex = MarkerUtils.chromIndex(vcfRec, vcfRec.substring(0, tabs.get(0)));
        int pos = Integer.parseInt(vcfRec.substring(tabs.get(0)+1, tabs.get(1)));
        if (stripOptionalFields) {
            sb.append('.');                                 // ID
            sb.append(vcfRec, tabs.get(2), tabs.get(4));    // REF, ALT
            sb.append("\t.\t.\t.");                         // QUAL, FILTER, INFO
        }
        else {
            sb.append(vcfRec, tabs.get(1) + 1, tabs.get(7));
        }
        short nAlleles = (short) nAlleles(vcfRec, tabs.get(3)+1, tabs.get(4));
        return new Marker(chromIndex, pos, nAlleles, sb.toString());
    }

    private static int nAlleles(String rec, int altStart, int altEnd) {
        if ((altEnd-altStart)==1 && rec.charAt(altStart)==Const.MISSING_DATA_CHAR) {
            return 1;
        }
        else {
            int nAlleles = 2;
            for (int j=altStart; j<altEnd; ++j) {
                if (rec.charAt(j)==',') {
                    ++nAlleles;
                }
            }
            return nAlleles;
        }
    }

    /**
     * Returns the VCF CHROM field.
     * @return the VCF CHROM field
     */
    public String chromID() {
        return ChromIds.instance().id(chromIndex);
    }

    /**
     * Returns the chromosome index.
     * @return the chromosome index
     */
    public int chromIndex() {
        return chromIndex;
    }

    /**
     * Returns the VCF POS field
     * @return the VCF POS field
     */
    public int pos() {
        return pos;
    }

    /**
     * Returns the VCF ID field.
     * @return the VCF ID field
     */
    public String id() {
        return fields.substring(0, fields.indexOf('\t'));
    }

    /**
     * Returns the tab-separated VCF REF and ALT fields
     * @return the tab-separated VCF REF and ALT fields
     */
    public String alleles() {
        int refStart = fields.indexOf('\t') + 1;
        int altStart = fields.indexOf('\t', refStart) + 1;
        int altEnd = fields.indexOf('\t', altStart);
        return fields.substring(refStart, altEnd);
    }

    /**
     * Returns the number of nucleotides in the reference allele.
     * @return the number of nucleotides in the reference allele
     */
    public int nRefBases() {
        int refStart = fields.indexOf('\t') + 1;
        int refEnd = fields.indexOf('\t', refStart);
        return (refEnd - refStart);
    }

    /**
     * Returns the number of alleles in the VCF REF and ALT fields.

     * @return the number of alleles in the VCF REF and ALT fields
     */
    public int nAlleles() {
        return nAlleles & 0xffff;
    }

    /**
     * Returns the minimum number of bits required to store a non-missing
     * allele.
     * @return the minimum number of bits required to store a non-missing
     * allele
     */
    public int bitsPerAllele() {
        return Integer.SIZE - Integer.numberOfLeadingZeros(nAlleles()-1);
    }

    /**
     * Returns the VCF QUAL field.
     * @return the VCF QUAL field
     */
    public String qual() {
        int filterEnd = fields.lastIndexOf('\t');
        int qualEnd = fields.lastIndexOf('\t', (filterEnd-1));
        int qualStart = fields.lastIndexOf('\t', (qualEnd-1)) + 1;
        return fields.substring(qualStart, qualEnd);
    }

    /**
     * Returns the VCF FILTER field.
     * @return the VCF FILTER field.
     */
    public String filter() {
        int filterEnd = fields.lastIndexOf('\t');
        int filterStart = fields.lastIndexOf('\t', (filterEnd-1)) + 1;
        return fields.substring(filterStart, filterEnd);
    }

    /**
     * Returns the VCF INFO field.
     * @return the VCF INFO field.
     */
    public String info() {
        return fields.substring(fields.lastIndexOf('\t') + 1);
    }

    /**
     * Returns the tab-delimited ID, REF, ALT, QUAL, FILTER, and INFO fields
     * @return the tab-delimited ID, REF, ALT, QUAL, FILTER, and INFO fields
     */
    public String fields() {
        return fields;
    }

    /**
     * <p>Returns the hash code value for this object.
     * The hash code is defined by the following calculation:
     * </p>
     *   <pre>
     *   int hash = 5;
     *   hash = 29 * hash + this.chromIndex();
     *   hash = 29 * hash + this.pos();
     *   hash = 29 * hash + this.alleles().hashCode();
    *   </pre>
     *
     * @return the hash code value for this marker
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 29 * hash + chromIndex;
        hash = 29 * hash + pos;
        hash = 29 * hash + alleles().hashCode();
        return hash;
    }

//    /**
//     * Returns the string that is obtained by reordering the marker's
//     * alleles in {@code this.fields()} so that the original {@code j}-th
//     * allele is the {@code newIndices[j]}-th allele in the returned string.
//     *
//     * @param indexMap an array whose {@code j}-th element is the
//     * index of the {@code}-th allele after re-ordering alleles
//     * @param flipStrand {@code true} if each allele in {@code marker.fields()}
//     * should be replaced with the allele on the complementary DNA strand
//     * @return the {@code this.fields()} string with re-ordered alleles
//     * @throws IllegalArgumentException if {@code indexMap} is not a
//     * permutation of this marker's allele indices
//     * @throws IllegalArgumentException if any allele is a symbolic
//     * allele
//     * @throws IllegalArgumentException if {@code (flipStrand == true)} and
//     * any character of a marker allele is not {@code 'A'}, {@code 'C'},
//     * {@code 'G'}, {@code 'T'}, {@code 'N'}, or {@code '*'}
//     * @throws NullPointerException if {@code (indexMap == null)}
//     */
//    public Marker reorderAlleles(int[] indexMap, boolean flipStrand) {
//        String newFields = AlleleMapper.reorderAltAlleles(this, indexMap, flipStrand);
//        return new Marker(chromIndex, pos, nAlleles, newFields);
//    }

    /**
     * Returns {@code true} if the specified object is a
     * {@code Marker} with the same chromosome, position, and alleles,
     * and returns {@code false} otherwise.
     * <p>
     * The return value is defined by the following calculation:
     * </p>
     *   <pre>
        if (this==obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Marker other = (Marker) obj;
        if (this.chromIndex() != other.chromIndex()) {
            return false;
        }
        if (this.pos() != other.pos()) {
            return false;
        }
        return this.alleles().equals(other.alleles());
     *  </pre>
     *
     * @param obj object to be compared with {@code this} for equality
     *
     * @return {@code true} if the specified object is a {@code Marker} with
     * the same chromosome, position, and alleles
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
        final Marker other = (Marker) obj;
        if (this.chromIndex != other.chromIndex) {
            return false;
        }
        if (this.pos != other.pos) {
            return false;
        }
        return this.alleles().equals(other.alleles());
    }

    /**
     * Compares this marker with the specified marker
     * for order, and returns a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal to,
     * or greater than the specified marker. Markers are compared
     * using the values returned by the {@code chromIndex()}, {@code pos()},
     * and {@code alleles()} methods. The returned value is defined
     * by the following calculation:
     *   <pre>
     *   if (this.chromIndex() != other.chromIndex()) {
     *       return (this.chromIndex &lt; other.chromIndex()) ? -1 : 1;
     *   }
     *   if (this.pos() != other.pos()) {
     *       return (this.pos &lt; other.pos()) ? -1 : 1;
     *   }
     *   return this.alleles().compareTo(other.alleles());
     *   </pre>
     *
     * @param other the {@code Marker} to be compared
     * @return a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal,
     * or greater than the specified marker
     */
    @Override
    public int compareTo(Marker other) {
        if (this.chromIndex != other.chromIndex) {
            return (this.chromIndex < other.chromIndex) ? -1 : 1;
        }
        if (this.pos != other.pos) {
            return (this.pos < other.pos) ? -1 : 1;
        }
        return this.alleles().compareTo(other.alleles());
    }

     /**
      * Returns a string equal to the first eight tab-delimited fields
      * of a VCF record corresponding to this marker (the CHROM, POS, ID,
      * REF, ALT, QUAL, FILTER, and INFO fields).
      *
      * @return a string equal to the first eight tab-delimited fields
      * of a VCF record corresponding to this marker
      */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(50);
        MarkerUtils.appendFirst8Fields(this, sb);
        return sb.toString();
    }
}

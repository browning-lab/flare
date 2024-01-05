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
import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Arrays;

/**
 * <p>Class {@code Marker} represents a VCF record's CHROM, POS, ID, REF,
 * ALT, QUAL, FILTER, and INFO fields. The number of alleles in the VCF
 * record must be less than or equal to {@code Marker.MAX_N_ALLELES}.</p>
 *
 * <p>Instances of class {@code Marker} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Marker implements Comparable<Marker> {

    public  static final short MAX_N_ALLELES = 0xff;
    private static final short STORED_N_ALLELES_MASK = 0xff;
    private static final short INDEXED_N_ALLELES_MASK = 0b111;

    private static final String[] SNV_PERMS = MarkerUtils.snvPerms();

    private static final short ID_STORED = (short) (1<<15);
    private static final short ALLELES_STORED = (short) (1<<14);
    private static final short QUAL_STORED = (short) (1<<13);
    private static final short FILTER_STORED = (short) (1<<12);
    private static final short INFO_STORED = (short) (1<<11);
    private static final short FLAG_BITS =
            ID_STORED | ALLELES_STORED | QUAL_STORED | FILTER_STORED | INFO_STORED;

    private final short chromIndex;
    private final int pos;
    private final short fieldInfo;
    private final String fields;

    private Marker(short chromIndex, int pos, short fieldInfo, String fields) {
        this.chromIndex = chromIndex;
        this.pos = pos;
        this.fieldInfo = fieldInfo;
        this.fields = fields;
    }

    /**
     * Constructs a new {@code BasicMarker} instance from the specified
     * string VCF record.  If {@code storeInfo == false} and a VCF record
     * INFO/END field is present, the INFO/END field will be stored.
     * @param vcfRec a VCF record
     * @param filter a filter for the VCF record's ID, QUAL, FILTER, and INFO
     * subfields
     * @throws IllegalArgumentException if the specified VCF record does not
     * contain at least 8 tab characters
     * @throws IllegalArgumentException if the VCF CHROM field contains
     * whitespace
     * @throws IllegalArgumentException if the specified VCF record has more
     * than 255 alleles
     * @throws IndexOutOfBoundsException if the index of the VCF CHROM field
     * exceeds {@code Short.MAX_VALUE}
     * @throws NullPointerException if
     * {@code (vcfRecord == null) || (filter==null)}
     * @throws NumberFormatException if the VCF record POS field is not a
     * parsable integer
     */
    public Marker(String vcfRec, VcfFieldFilter filter) {
        IntList tabs = MarkerUtils.first8TabIndices(vcfRec);
        short fInfo = (short) 0;
        StringBuilder sb = new StringBuilder();
        this.chromIndex = MarkerUtils.chromIndex(vcfRec, vcfRec.substring(0, tabs.get(0)));
        this.pos = Integer.parseInt(vcfRec.substring(tabs.get(0)+1, tabs.get(1)));
        if (filter.storeId(vcfRec, tabs.get(1)+1, tabs.get(2), sb)) {
             fInfo |= ID_STORED;
        }
        fInfo = storeAlleles(vcfRec, tabs.get(2)+1, tabs.get(4), fInfo, ALLELES_STORED, sb);
        if (filter.storeQual(vcfRec, tabs.get(4)+1, tabs.get(5),
                ((fInfo & FLAG_BITS)!=0), sb)) {
             fInfo |= QUAL_STORED;
        }
        if (filter.storeFilter(vcfRec, tabs.get(5)+1, tabs.get(6),
                ((fInfo & FLAG_BITS)!=0), sb)) {
             fInfo |= FILTER_STORED;
        }
        if (filter.storeInfo(vcfRec, tabs.get(6)+1, tabs.get(7),
                ((fInfo & FLAG_BITS)!=0), sb)) {
             fInfo |= INFO_STORED;
        }
        this.fields = (fInfo & FLAG_BITS)==0 ? null : sb.toString();
        this.fieldInfo = fInfo;
    }

    private static short storeAlleles(String vcfRec, int start, int end,
            short fieldInfo, short fieldInfoMask, StringBuilder sb) {
        String refAltAlleles = new String(vcfRec.substring(start, end));
        int snvIndex = snvIndex(refAltAlleles);
        if (snvIndex>=0) {
            int nAlleles = refAltAlleles.endsWith("\t.") ? 1 : ((end - start + 1) >> 1);
            fieldInfo |= (snvIndex<<3);
            fieldInfo |= nAlleles;
        }
        else {
            int tabIndex = refAltAlleles.indexOf(Const.tab);
            if ((tabIndex+1)==refAltAlleles.length()) {
                String s = "ERROR: missing ALT field: "
                        + MarkerUtils.truncate(vcfRec, 80);
                throw new IllegalArgumentException(s);
            }
            int nAlleles = nAlleles(refAltAlleles, tabIndex);
            if (nAlleles > STORED_N_ALLELES_MASK) {
                throw new IndexOutOfBoundsException(String.valueOf(nAlleles));
            }
            if ((fieldInfo & FLAG_BITS)!=0) {
                sb.append(Const.tab);
            }
            sb.append(refAltAlleles);
            fieldInfo |= nAlleles;
            fieldInfo |= fieldInfoMask;
        }
        return fieldInfo;
    }

    private static int snvIndex(String refAndAlt) {
        int index = Arrays.binarySearch(SNV_PERMS, refAndAlt);
        if (index<0) {
            index = -index-1;
        }
        if (index==SNV_PERMS.length) {
            return -1;
        }
        else {
            return SNV_PERMS[index].startsWith(refAndAlt) ? index : -1;
        }
    }

    private static int nAlleles(String refAltAlleles, int tabIndex) {
        int startAllele = tabIndex + 1;
        if (startAllele==(refAltAlleles.length()-1)
                && refAltAlleles.charAt(startAllele)==Const.MISSING_DATA_CHAR) {
            return 1;
        }
        else {
            int nAlleles = 2;
            startAllele = refAltAlleles.indexOf(Const.comma, startAllele) + 1;
            while (startAllele>0) {
                ++nAlleles;
                startAllele = refAltAlleles.indexOf(Const.comma, startAllele) + 1;
            }
            return nAlleles;
        }
    }

    /**
     * Returns the VCF CHROM field.
     * @return the VCF CHROM field
     */
    public String chrom() {
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
     * Returns {@code true} if the VCF ID field has non-missing
     * data, and returns {@code false} otherwise.
     *
     * @return {@code true} if the VCF ID field has non-missing
     * data, and returns {@code false} otherwise
     */
    public boolean hasIdData() {
        return (fieldInfo & ID_STORED)==ID_STORED;
    }

    /**
     * Returns {@code true} if the VCF QUAL field has non-missing
     * data, and returns {@code false} otherwise.
     *
     * @return {@code true} if the VCF QUAL field has non-missing
     * data, and returns {@code false} otherwise
     */
    public boolean hasQualData() {
        return (fieldInfo & QUAL_STORED)==QUAL_STORED;
    }

    /**
     * Returns {@code true} if the VCF FILTER field has non-missing
     * data, and returns {@code false} otherwise.
     *
     * @return {@code true} if the VCF FILTER field has non-missing
     * data, and returns {@code false} otherwise
     */
    public boolean hasFilterData() {
        return (fieldInfo & FILTER_STORED)==FILTER_STORED;
    }

    /**
     * Returns {@code true} if the VCF INFO field has non-missing
     * data, and returns {@code false} otherwise.
     *
     * @return {@code true} if the VCF INFO field has non-missing
     * data, and returns {@code false} otherwise
     */
    public boolean hasInfoData() {
        return (fieldInfo & INFO_STORED)==INFO_STORED;
    }

    /**
     * Returns the VCF ID field.
     * @return the VCF ID field
     */
    public String id() {
        if ((fieldInfo & ID_STORED)==ID_STORED) {
            int endIndex = fields.indexOf(Const.tab);
            return endIndex<0 ? fields : fields.substring(0, endIndex);
        }
        else {
            return Const.MISSING_DATA_STRING;
        }
    }

    /**
     * Returns the tab-separated VCF REF and ALT fields
     * @return the tab-separated VCF REF and ALT fields
     */
    public String refAlt() {
        if ((fieldInfo & ALLELES_STORED)==ALLELES_STORED) {
            int start = 0;
            if ((fieldInfo & ID_STORED)==ID_STORED) {
                start = fields.indexOf(Const.tab) + 1;      // start of REF
            }
            int end = fields.indexOf(Const.tab, start);     // end of REF
            end = fields.indexOf(Const.tab, end+1);         // end of ALT
            return end<0 ? fields.substring(start) : fields.substring(start, end);
        }
        else {
            int snvIndex = (fieldInfo>>>3) & 0xff;
            int nAlleles = fieldInfo & INDEXED_N_ALLELES_MASK;
            if (nAlleles==1) {
                return SNV_PERMS[snvIndex];
            }
            else {
                int nChars =(nAlleles<<1) - 1;
                return SNV_PERMS[snvIndex].substring(0, nChars);
            }
        }
    }

    /**
     * Returns the number of alleles for the marker, including the REF
     * allele.
     * @return the number of alleles for the marker, including the REF
     * allele
     */
    public int nAlleles() {
        if ((fieldInfo & ALLELES_STORED)==ALLELES_STORED) {
            return (fieldInfo & STORED_N_ALLELES_MASK);
        }
        else {
            return fieldInfo & INDEXED_N_ALLELES_MASK;
        }
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

    private int qualStartIndex() {
        int start = 0;
        if ((fieldInfo & ID_STORED)==ID_STORED) {
            start = fields.indexOf(Const.tab) + 1;         // skip ID field
        }
        if ((fieldInfo & ALLELES_STORED)==ALLELES_STORED) {
            start = fields.indexOf(Const.tab, start) + 1;  // skip REF field
            start = fields.indexOf(Const.tab, start) + 1;  // skip ALT field
        }
        return start;
    }

    private String extractField(int start) {
        int end = fields.indexOf(Const.tab, start);
        return end<0 ? fields.substring(start) : fields.substring(start, end);
    }

    /**
     * Returns the VCF QUAL field.
     * @return the VCF QUAL field
     */
    public String qual() {
        if ((fieldInfo & QUAL_STORED)==QUAL_STORED) {
            return extractField(qualStartIndex());
        }
        else {
            return Const.MISSING_DATA_STRING;
        }
    }

    /**
     * Returns the VCF FILTER field.
     * @return the VCF FILTER field.
     */
    public String filter() {
        if ((fieldInfo & FILTER_STORED)==FILTER_STORED) {
            int start = qualStartIndex();
            if ((fieldInfo & QUAL_STORED)==QUAL_STORED) {
                start = fields.indexOf(Const.tab, start) + 1; // skip QUAL field
            }
            return extractField(start);
        }
        else {
            return Const.MISSING_DATA_STRING;
        }
    }

    /**
     * Returns the VCF INFO field.
     * @return the VCF INFO field.
     */
    public String info() {
        if ((fieldInfo & INFO_STORED)==INFO_STORED) {
            int start = qualStartIndex();
            if ((fieldInfo & QUAL_STORED)==QUAL_STORED) {
                start = fields.indexOf(Const.tab, start) + 1;  // skip QUAL field
            }
            if ((fieldInfo & FILTER_STORED)==FILTER_STORED) {
                start = fields.indexOf(Const.tab, start) + 1;  // skip FILTER field
            }
            return extractField(start);
        }
        else {
            return Const.MISSING_DATA_STRING;
        }
    }

    /**
     * <p>Returns the hash code value for this object.
     * The hash code is defined by the following calculation:
     * </p>
     * <pre>
     *   int hash = 5;
     *   hash = 29 * hash + chromIndex;
     *   hash = 29 * hash + this.pos;
     *   hash = 29 * hash + refAlt().hashCode();
     * </pre>
     *
     * @return the hash code value for this marker
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 29 * hash + chromIndex;
        hash = 29 * hash + pos;
        hash = 29 * hash + refAlt().hashCode();
        return hash;
    }

    /**
     * Returns {@code true} if the specified object is a
     * {@code Marker} with the same chromosome, position, and alleles,
     * and returns {@code false} otherwise.
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
        return this.refAlt().equals(other.refAlt());
    }

    /**
     * Compares this marker with the specified marker
     * for order, and returns a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal to,
     * or greater than the specified marker.  Markers are compared
     * in using the values returned by the {@code chromIndex()}, {@code pos()}
     * and {@code refAlt()} returned methods. The returned value is
     * defined by the following calculation:
     *   <pre>
     *   if (this.chromIndex() != other.chromIndex()) {
     *       return (this.chromIndex &lt; other.chromIndex()) ? -1 : 1;
     *   }
     *   if (this.pos() != other.pos()) {
     *       return (this.pos &lt; other.pos()) ? -1 : 1;
     *   }
     *   return this.refAlt().compareTo(other.refAlt());
     *  </pre>
     *
     * @param other the {@code Marker} to be compared
     * @return a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal,
     * or greater than the specified marker
     */
    @Override
    public int compareTo(Marker other) {
        if (this.chromIndex != other.chromIndex()) {
            return (this.chromIndex < other.chromIndex) ? -1 : 1;
        }
        if (this.pos != other.pos) {
            return (this.pos < other.pos) ? -1 : 1;
        }
        return this.refAlt().compareTo(other.refAlt());
    }

    /**
     * Writes a representation of the VCF record ID, REF, ALT, QUAL, FILTER,
     * and INFO fields to the specified output. The exact details of the
     * representation are unspecified and subject to change. The written data
     * can be read with the {@code Marker.readNonPosFields()} method.
     * @param out the output destination
     * @throws IOException if an I/O error occurs
     * @throws NullPointerException if {@code out == null}
     */
    public void writeNonPosFields(DataOutput out) throws IOException {
        out.writeShort(fieldInfo);
        if ((fieldInfo & FLAG_BITS)!=0) {
            out.writeUTF(fields);
        }
    }

    /**
     * Reads the VCF record ID, REF, ALT, QUAL, FILTER, and INFO fields and
     * returns a marker with these fields and the specified CHROM and POS
     * fields.  The contract for this method is unspecified if
     * {@code chromIndex} is not a valid chromosome index.
     * @param chromIndex the chromosome index
     * @param pos the chromosome position
     * @param in the input source
     * @return a {@code Marker}
     * @throws IOException if an I/O error occurs
     * @throws NullPointerException if {@code in == null}
     */
    public static Marker readNonPosFields(short chromIndex, int pos,
            DataInput in) throws IOException {
        short fieldInfo = (short) in.readShort();
        String fields = (fieldInfo & Marker.FLAG_BITS) == 0 ? null : in.readUTF();
        return new Marker(chromIndex, pos, fieldInfo, fields);
    }

     /**
      * Returns a string equal to the first five tab-delimited fields
      * of a VCF record corresponding to this marker (the CHROM, POS, ID,
      * REF, and ALT fields).
      *
      * @return a string equal to the first five tab-delimited fields
      * of a VCF record corresponding to this marker
      */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(50);
        sb.append(chrom());
        sb.append(Const.tab);
        sb.append(pos);
        sb.append(Const.tab);
        sb.append(id());
        sb.append(Const.tab);
        sb.append(refAlt());
        return sb.toString();
    }
}

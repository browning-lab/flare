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
package bref;

import beagleutil.ChromIds;
import blbutil.Const;
import blbutil.Filter;
import blbutil.Utilities;
import ints.CharArray;
import ints.IntArray;
import ints.UnsignedByteArray;
import java.io.DataInput;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.Samples;
import vcf.HapRefGTRec;
import vcf.VcfFieldFilter;

/**
 * <p>Class {@code Bref3Reader} contains methods for reading a bref3
 * (binary reference format version 3) file.
 * </p>
 * <p>Instances of class {@code Bref3Reader} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref3Reader {

    static final String READ_ERR = "Error reading file";
    private static final String[][] SNV_PERMS = snvPerms();

    private final Filter<Marker> markerFilter;
    private final String program;
    private final Bref3Header brefHeader;
    private final int[] includedHapIndices;
    private final int[] invIncludedHapIndices;
    private final Samples samples;

    private final byte[] byteBuffer;
    private final char[] hapToSeq;

    /**
     * Constructs a new {@code Bref3Reader} instance.
     * @param source the source bref3 file or {@code null} if the bref3 file
     * is read from stdin
     * @param dataIn a {@code DataInput} instance reading from a bref3 file
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of the specified bref3 file
     * @throws NullPointerException if {@code (dataIn == null)}
     */
    public Bref3Reader(File source, DataInput dataIn) {
        this(source, dataIn, Filter.acceptAllFilter(), Filter.acceptAllFilter());
    }

    /**
     * Constructs a new {@code Bref3Reader} instance.
     * @param source the source bref3 file or {@code null} if the bref3 file
     * is read from stdin
     * @param dataIn a {@code DataInput} instance reading from a bref3 file
     * @param sampleFilter a sample filter
     * @param markerFilter a marker filter
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of the specified bref3 file
     * @throws NullPointerException if
     * {@code (dataIn == null) || (sampleFilter == null) || (markerFilter == null)}
     */
    public Bref3Reader(File source, DataInput dataIn,
            Filter<String> sampleFilter, Filter<Marker> markerFilter) {
        if (markerFilter==null) {
            throw new NullPointerException("markerFilter==null");
        }
        this.brefHeader = new Bref3Header(source, dataIn, sampleFilter);
        this.program = brefHeader.program();
        this.includedHapIndices = brefHeader.filteredHapIndices();
        this.invIncludedHapIndices = brefHeader.invfilteredHapIndices();
        this.samples = brefHeader.samples();
        this.byteBuffer = new byte[2*invIncludedHapIndices.length];
        this.hapToSeq = new char[includedHapIndices.length];
        this.markerFilter = markerFilter;
    }

    /**
     * Returns the list of samples in the brefFile used to construct this
     * instance
     * @return the list of samples
     */
    Samples samples() {
        return samples;
    }

    /**
     * Returns the bref file program string in the brefFile used to
     * construct this instance
     * @return the program string
     */
    String program() {
        return program;
    }

    /**
     * Reads the next binary reference format data block.  The contract for
     * this method is undefined if the next byte read from the specified
     * {@code DataInput} object is not the first byte of a bref data block
     * or is not the first byte of the sentinal denoting no more
     * bref data blocks.
     * @param bref the bref file and associated file pointer
     * @param buffer the collection to which the next block of records
     * will be added
     * @throws NullPointerException if {@code bref == null || buffer == null}
     */
    void readBlock(DataInput bref, Collection<RefGTRec> buffer) {
        try {
            int nRecs = Integer.MAX_VALUE;
            while (buffer.isEmpty() && nRecs!=0) {
                nRecs = bref.readInt();
                if (nRecs!=0) {
                    readBlock(bref, buffer, nRecs);
                }
            }
        } catch (IOException ex) {
            Utilities.exit(ex, READ_ERR);
        }
    }

    private void readBlock(DataInput dataIn, Collection<RefGTRec> buffer,
            int nRecs) throws IOException {
        assert nRecs!= 0;
        String chrom = dataIn.readUTF();
        int chromIndex = ChromIds.instance().getIndex(chrom);
        int nSeq = dataIn.readUnsignedShort();
        dataIn.readFully(byteBuffer);
        IntArray hapToSeq = hapToSeq(byteBuffer);
        for (int j=0; j<nRecs; ++j) {
            Marker marker = readMarker(dataIn, chromIndex);
            byte flag = dataIn.readByte();
            if (flag==0) {
                RefGTRec rec = readHapRecord(dataIn, marker, samples,
                        hapToSeq, nSeq);
                if (markerFilter.accept(marker)) {
                    buffer.add(rec);
                }
            }
            else if (flag==1) {
                RefGTRec rec = readAlleleRecord(dataIn, marker, samples);
                if (markerFilter.accept(marker)) {
                    buffer.add(rec);
                }
            }
            else {
                Utilities.exit(READ_ERR);
            }
        }
    }

    private IntArray hapToSeq(byte[] byteBuffer) {
        for (int k=0; k<hapToSeq.length; ++k) {
            int offset = (includedHapIndices[k]<<1);
            int b1 = byteBuffer[offset] & 0xff;
            int b2 = byteBuffer[offset+1] & 0xff;
            hapToSeq[k] = (char) ((b1 << 8) + b2);
        }
        return new CharArray(hapToSeq);
    }

    private RefGTRec readHapRecord(DataInput bref, Marker marker,
            Samples samples, IntArray hapToSeq, int nSeq) throws IOException {
        bref.readFully(byteBuffer, 0, nSeq);
        IntArray seqToAllele = new UnsignedByteArray(byteBuffer, 0, nSeq);

//        // following code can be uncommented to check data consistency
//        if (seqToAllele.max() >= marker.nAlleles()) {
//            throw new IllegalArgumentException("inconsistent data");
//        }

        return new HapRefGTRec(marker, samples, hapToSeq, seqToAllele);
    }

    private static Marker readMarker(DataInput dataIn, int chromIndex)
            throws IOException {
        int pos = dataIn.readInt();
        String id = readByteLengthStringArrayAndJoin(dataIn, Const.semicolon);
        byte alleleCode = dataIn.readByte();
        if (alleleCode == -1) {
            String[] strAlleles = readStringArray(dataIn);
            String refAndAltFields = refAndAltFields(strAlleles);
            int end = dataIn.readInt();
            String vcfRecPrefix = vcfRecPrefix(chromIndex, pos, id,
                    refAndAltFields, end);
            boolean storeId = true;
            VcfFieldFilter filter = new VcfFieldFilter(storeId, false, false, false);
            return new Marker(vcfRecPrefix, filter);
        }
        else {
            int nAlleles = 1 + (alleleCode & 0b11);
            int permIndex = alleleCode >> 2;
            String[] strAlleles = alleleString(permIndex, nAlleles);
            String refAndAltFields = refAndAltFields(strAlleles);
            int end = -1;
            String vcfRecPrefix = vcfRecPrefix(chromIndex, pos, id,
                    refAndAltFields, end);
            boolean storeId = true;
            VcfFieldFilter filter = new VcfFieldFilter(storeId, false, false, false);
            return new Marker(vcfRecPrefix, filter);
        }
    }

    private static String refAndAltFields(String[] strAlleles) {
        StringBuilder sb = new StringBuilder();
        sb.append(strAlleles[0]);
        if (strAlleles.length==1) {
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=1; j<strAlleles.length; ++j) {
                sb.append(j==1 ? Const.tab : Const.comma);
                sb.append(strAlleles[j]);
            }
        }
        return sb.toString();
    }

    private static String vcfRecPrefix(int chrom, int pos, String id,
            String refAndAltFields, int end) {
        StringBuilder sb = new StringBuilder(64);
        sb.append(ChromIds.instance().id(chrom));
        sb.append(Const.tab);
        sb.append(pos);
        sb.append(Const.tab);
        sb.append(id);
        sb.append(Const.tab);
        sb.append(refAndAltFields);
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR); // QUAL
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR); // FILTER
        sb.append(Const.tab);
        if (end>=0) {
            sb.append("END=");                  // INFO
            sb.append(end);                     // INFO
        }
        else {
            sb.append(Const.MISSING_DATA_CHAR); // INFO
        }
        sb.append(Const.tab);
        sb.append("GT");
        return sb.toString();
    }

    private RefGTRec readAlleleRecord(DataInput di, Marker marker,
            Samples samples) throws IOException {
        int nAlleles = marker.nAlleles();
        int[][] hapIndices = new int[nAlleles][];
        for (int j=0; j<nAlleles; ++j) {
            hapIndices[j] = readAlleleCodedHapList(di);
        }
        return RefGTRec.alleleRefGTRec(marker, samples, hapIndices);
    }

    private int[] readAlleleCodedHapList(DataInput dataInput) throws IOException {
        int length = dataInput.readInt();
        if (length == -1) {
            return null;
        }
        else {
            int[] ia = new int[length];
            int index = 0;
            for (int j=0; j<ia.length; ++j) {
                int hap = dataInput.readInt();
                if (invIncludedHapIndices[hap]>=0) {
                    ia[index++] = invIncludedHapIndices[hap];
                }
            }
            return index<ia.length ? Arrays.copyOf(ia, index) : ia;
        }
    }

    private static String readByteLengthStringArrayAndJoin(DataInput dataIn,
            char delim) throws IOException {
        int length = dataIn.readUnsignedByte();
        if (length <= 0) {
            return Const.MISSING_DATA_STRING;
        }
        else {
            StringBuilder sb = new StringBuilder();
            for (int j=0; j<length; ++j) {
                if (j>0) {
                    sb.append(delim);
                }
                sb.append(dataIn.readUTF());
            }
            return sb.toString();
        }
    }

    static String[] readStringArray(DataInput di) throws IOException {
        int length = di.readInt();
        return readStringArray(di, length);
    }

    /* Returns null if length is negative */
    private static String[] readStringArray(DataInput dataIn, int length)
            throws IOException {
        if (length<0) {
            return null;
        } else if (length==0) {
            return new String[0];
        } else {
            String[] sa=new String[length];
            for (int j=0; j<sa.length; ++j) {
                sa[j]=dataIn.readUTF();
            }
            return sa;
        }
    }

    static String[][] snvPerms() {
        String[] bases=new String[]{"A", "C", "G", "T"};
        List<String[]> perms=new ArrayList<>(24);
        permute(new String[0], bases, perms);
        return perms.toArray(new String[0][]);
    }

    private static void permute(String[] start, String[] end,
            List<String[]> perms) {
        if (end.length==0) {
            perms.add(start);
        }
        else {
            for (int j=0; j<end.length; ++j) {
                String[] newStart = Arrays.copyOf(start, start.length + 1);
                newStart[start.length] = end[j];

                String[] newEnd = new String[end.length - 1];
                if (j > 0) {
                    System.arraycopy(end, 0, newEnd, 0, j);
                }
                if (j < newEnd.length) {
                    System.arraycopy(end, j+1, newEnd, j, (newEnd.length - j));
                }
                permute(newStart, newEnd, perms);
            }
        }
    }

    /**
     * Returns an array that is obtained by taking the first {@code length}
     * elements of the specified permutation of "A", "C", "G", and "T".
     * The list of 24 permutations of "A", "C", "G", and "T" are sorted
     * in lexicographic order.
     * @param permIndex an index of a permutation of the bases "A",
     * "C", "G", and "T"
     * @param length the number of elements in the returned array
     * @return an array that is obtained by taking the first {@code length}
     * elements of the specified permutation of "A", "C", "G", and "T"
     * @throws IndexOutOfBoundsException if
     * {@code permIndex < 0 || permIndex >= 24}
     * @throws IndexOutOfBoundsException if {@code length < 0 || length >= 4}
     */
    private static String[] alleleString(int permIndex, int length) {
        return Arrays.copyOf(SNV_PERMS[permIndex], length);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        return sb.toString();
    }
}

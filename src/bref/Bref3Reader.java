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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import vcf.BasicMarker;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.Samples;
import vcf.SeqCodedRefGTRec;

/**
 * <p>Class {@code Bref3Reader} contains methods for reading a bref 3
 * (binary reference format) file.
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
    private final Samples samples;
    private final int nHaps;
    private final byte[] byteBuffer;

    /**
     * Constructs a new {@code Bref3Reader} instance.
     * @param bref a {@code DataInput} instance reading from a bref v3 file
     * @param markerFilter a marker filter or {@code null}
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of the specified bref v3 file
     * @throws NullPointerException if {@code file == null}
     */
    public Bref3Reader(DataInput bref, Filter<Marker> markerFilter) {
        if (markerFilter == null) {
            markerFilter = Filter.acceptAllFilter();
        }
        String[] sampleIds = null;
        String programString = null;
        try {
            readAndCheckMagicNumber(bref);
            programString = readString(bref);
            sampleIds = readStringArray(bref);
        } catch (IOException ex) {
            Utilities.exit(ex, READ_ERR);
        }
        boolean[] isDiploid = new boolean[sampleIds.length];
        Arrays.fill(isDiploid, true);
        this.program = programString;
        this.samples = Samples.fromIds(sampleIds, isDiploid);
        this.markerFilter = markerFilter;
        this.nHaps = 2*samples.size();
        this.byteBuffer = new byte[2*nHaps];
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
     * Returns the marker filter
     * @return the marker filter
     */
    Filter<Marker> markerFilter() {
        return markerFilter;
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

    private void readBlock(DataInput bref, Collection<RefGTRec> buffer,
            int nRecs) throws IOException {
        assert nRecs!= 0;
        String chrom = readString(bref);
        int chromIndex = ChromIds.instance().getIndex(chrom);
        int nSeq = bref.readUnsignedShort();
        bref.readFully(byteBuffer);
        IntArray hapToSeq = new CharArray(byteBuffer);
        for (int j=0; j<nRecs; ++j) {
            Marker marker = readMarker(bref, chromIndex);
            byte flag = bref.readByte();
            if (flag==0) {
                RefGTRec rec = readSeqCodedRecord(bref, marker, samples,
                        hapToSeq, nSeq);
                if (markerFilter.accept(marker)) {
                    buffer.add(rec);
                }
            }
            else if (flag==1) {
                RefGTRec rec = readHapCodedRec(bref, marker, samples);
                if (markerFilter.accept(marker)) {
                    buffer.add(rec);
                }
            }
            else {
                Utilities.exit(READ_ERR);
            }
        }
    }

    private RefGTRec readSeqCodedRecord(DataInput bref, Marker marker,
            Samples samples, IntArray hapToSeq, int nSeq) throws IOException {
        bref.readFully(byteBuffer, 0, nSeq);
        IntArray seqToAllele = new UnsignedByteArray(byteBuffer, 0, nSeq);

//        // following code can be uncommented to check data consistency
//        if (seqToAllele.max() >= marker.nAlleles()) {
//            throw new IllegalArgumentException("inconsistent data");
//        }

        return new SeqCodedRefGTRec(marker, samples, hapToSeq, seqToAllele);
    }

    private static void readAndCheckMagicNumber(DataInput di) throws IOException {
        int magicNumber=di.readInt();
        if (magicNumber!=AsIsBref3Writer.MAGIC_NUMBER_V3) {
            String s = "ERROR: Unrecognized input file.  Was input file created "
                    + Const.nl + "with a different version of the bref program?";
            Utilities.exit(s);
        }
    }

    private static Marker readMarker(DataInput di, int chromIndex)
            throws IOException {
        int pos = di.readInt();
        String[] ids = readByteLengthStringArray(di);
        byte alleleCode = di.readByte();
        if (alleleCode == -1) {
            String[] strAlleles = readStringArray(di);
            int end = di.readInt();
            return new BasicMarker(chromIndex, pos, ids, strAlleles, end);
        }
        else {
            int nAlleles = 1 + (alleleCode & 0b11);
            int permIndex = alleleCode >> 2;
            String[] strAlleles = alleleString(permIndex, nAlleles);
            int end = -1;
            return new BasicMarker(chromIndex, pos, ids, strAlleles, end);
        }
    }

    private static RefGTRec readHapCodedRec(DataInput di, Marker marker,
            Samples samples) throws IOException {
        int nAlleles = marker.nAlleles();
        int[][] hapIndices = new int[nAlleles][];
        for (int j=0; j<nAlleles; ++j) {
            hapIndices[j] = readIntArray(di);
        }
        return RefGTRec.hapCodedInstance(marker, samples, hapIndices);
    }

    private static int[] readIntArray(DataInput di) throws IOException {
        int length = di.readInt();
        if (length == -1) {
            return null;
        }
        else {
            int[] ia = new int[length];
            byte[] ba = new byte[4*length]; // will overflow if 4*length >= 2^31
            di.readFully(ba);
            for (int j=0; j<ba.length; j+=4) {
                ia[j/4] = (((ba[j] & 0xff) << 24) | ((ba[j+1] & 0xff) << 16) |
                            ((ba[j+2] & 0xff) << 8) | (ba[j+3] & 0xff));
            }
            return ia;
        }
    }

    private static String readString(DataInput di) throws IOException {
        return di.readUTF();
    }

    private static String[] readByteLengthStringArray(DataInput di)
            throws IOException {
        int length = di.readUnsignedByte();
        return readStringArray(di, length);
    }

    static String[] readStringArray(DataInput di) throws IOException {
        int length = di.readInt();
        return readStringArray(di, length);
    }

    /* Returns null if length is negative */
    private static String[] readStringArray(DataInput di, int length)
            throws IOException {
        if (length<0) {
            return null;
        } else if (length==0) {
            return new String[0];
        } else {
            String[] sa=new String[length];
            for (int j=0; j<sa.length; ++j) {
                sa[j]=readString(di);
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

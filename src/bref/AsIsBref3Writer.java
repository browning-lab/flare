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
package bref;

import beagleutil.ChromIds;
import blbutil.Const;
import blbutil.FileUtil;
import blbutil.Utilities;
import ints.IntArray;
import ints.IntList;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import vcf.Marker;
import vcf.MarkerUtils;
import vcf.RefGTRec;
import vcf.Samples;

/**
 * <p>Class {@code AsIsBref3Writer} writes VCF data with phased, non-missing
 * genotypes to a binary reference format v3 (bref) file.  Each record that
 * is written will have the same internal representation (allele-coded or
 * sequence-coded) as the {@code RefGTRec} passed to the {@code write()} method.
 * The {@code close()} method must be called after the last invocation of the
 * {@code write()} method in order to ensure that any buffered data are
 * written to the output binary reference file.
 * </p>
 * <p>Instances of class {@code AsIsBref3Writer} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AsIsBref3Writer implements BrefWriter {

    /**
     * The end of file code for a bref file.
     */
    public static final int END_OF_DATA = 0;

    /**
     * The integer denoting denoting the end of the index in a bref file
     */
    public static final long END_OF_INDEX = -999_999_999_999_999L;

    /**
     * The initial integer in a bref version 3 file.
     */
    public static final int MAGIC_NUMBER_V3 = 2055763188;

    /**
     * The byte value denoting a sequence coded record
     */
    public static final byte SEQ_CODED = 0;

    /**
     * The byte value denoting an allele coded record
     */
    public static final byte ALLELE_CODED = 1;

    private final String WRITE_ERR = "Error writing file";
    private final String CONTIGUITY_ERR = "Error: chromosomes not contiguous";

    private final Set<String> BASES_SET = basesSet();
    private final String[][] SNV_PERMS = Bref3Reader.snvPerms();
    private final Comparator<String[]> ALLELES_COMP = allelesComparator();

    public final int MAX_SAMPLES = (1 << 29) - 1;

    private int lastChromIndex = -1;
    private IntArray hap2Seq = null;
    private long bytesWritten = 0;

    private final File bref;
    private final Samples samples;
    private final int nHaps;
    private final List<RefGTRec> recBuffer;
    private final List<BrefBlock> index;

    private final DataOutputStream brefOut;
    private final ByteArrayOutputStream baos;
    private final DataOutputStream utf8Buffer;

    /**
     * Constructs a new {@code AsIsBref4Writer} for the specified data.
     * The Java virtual machine will exit with an error message if an I/O
     * error occurs during object construction
     *
     * @param program the name of the program which is creating the
     * binary reference file.
     * @param samples the list of samples whose genotype data will
     * be written in binary reference format
     * @param brefFile name of the output binary reference file or
     * {@code null} if the output should be directed to standard output
     *
     * @throws IllegalArgumentException if {
     * {@code samples.size() > AsIsBref4Writer.MAX_SAMPLES}
     * @throws NullPointerException if {@code program == null || samples == null}
     */
    public AsIsBref3Writer(String program, Samples samples, File brefFile) {
        if (program==null) {
            throw new NullPointerException(String.class.toString());
        }
        if (samples.size() > MAX_SAMPLES) {
            throw new IllegalArgumentException(String.valueOf(samples.size()));
        }
        this.bref = brefFile;
        this.samples = samples;
        this.nHaps = 2*samples.size();
        this.recBuffer = new ArrayList<>(500);
        this.index = new ArrayList<>(500);
        this.brefOut = dataOutputStream(bref);
        this.baos = new ByteArrayOutputStream(100);
        this.utf8Buffer = new DataOutputStream(baos);
        try {
            brefOut.writeInt(MAGIC_NUMBER_V3);
            bytesWritten += Integer.BYTES;
            writeString(program, brefOut);
            writeStringArray(samples.ids(), brefOut);
        } catch (IOException ex) {
            Utilities.exit(ex, WRITE_ERR);
        }
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public void write(RefGTRec rec) {
        if (rec.samples().equals(samples)==false) {
            Utilities.exit("ERROR: inconsistent data");
        }
        if (startNewBlock(rec)) {
            writeAndClearBuffer();
        }
        recBuffer.add(rec);
    }

     private boolean startNewBlock(RefGTRec rec) {
        boolean startNewBlock = false;
        int chromIndex = rec.marker().chromIndex();
        if (chromIndex!=lastChromIndex) {
            lastChromIndex = chromIndex;
            hap2Seq = null;
            startNewBlock = true;
        }
        if (rec.isAlleleRecord()==false) {
            if (hap2Seq==null) {
                hap2Seq = rec.map(0);
            }
            else if (rec.map(0)!=hap2Seq) {
                hap2Seq = rec.map(0);
                startNewBlock = true;
            }
        }
        return startNewBlock;
     }

    @Override
    public void close() {
        try {
            writeAndClearBuffer();
            brefOut.writeInt(END_OF_DATA);
            bytesWritten += Integer.BYTES;

            long indexOffset = bytesWritten;
            writeIndex();
            brefOut.writeLong(indexOffset);
            bytesWritten += Long.BYTES;

            brefOut.close();

        } catch (IOException ex) {
            Utilities.exit(ex, "Error closing file");
        }
    }

    private void writeAndClearBuffer() {
        if (recBuffer.isEmpty()==false) {
            try {
                Marker m = recBuffer.get(0).marker();
                index.add(new BrefBlock(m.chromIndex(), m.pos(), bytesWritten));
                brefOut.writeInt(recBuffer.size());
                bytesWritten += Integer.BYTES;
                writeString(m.chromID(), brefOut);
                writeHapToSeq();
                for (int j=0, n=recBuffer.size(); j<n; ++j) {
                    RefGTRec rec = recBuffer.get(j);
                    if (rec.isAlleleRecord()) {
                        writeAlleleCodedRec(rec);
                    }
                    else {
                        writeSeqCodedRec(rec);
                    }
                }
                recBuffer.clear();
            } catch (IOException ex) {
                Utilities.exit(ex, WRITE_ERR);
            }
        }
    }

    private void writeHapToSeq() throws IOException {
        RefGTRec rec = null;
        for (int j=0, n=recBuffer.size(); j<n && rec==null; ++j) {
            RefGTRec candidate = recBuffer.get(j);
            if (candidate.isAlleleRecord()==false) {
                rec = candidate;
            }
        }
        if (rec==null) {
            brefOut.writeChar(0);
            for (int j=0; j<nHaps; ++j) {
                brefOut.writeChar(0);
            }
        }
        else {
            assert rec.nMaps()==2;
            IntArray hapToSeq = rec.map(0);
            IntArray seqToAllele = rec.map(1);
            brefOut.writeChar(seqToAllele.size());
            for (int j=0, n=hapToSeq.size(); j<n; ++j) {
                brefOut.writeChar(hapToSeq.get(j));
            }
        }
        bytesWritten += (nHaps + 1)*Character.BYTES;
    }

    private void writeSeqCodedRec(RefGTRec rec) throws IOException {
        if (rec.marker().nAlleles() >= 256) {
            Utilities.exit("ERROR: more than 256 alleles: " + rec.marker());
        }
        assert rec.nMaps()==2;
        IntArray seq2Allele = rec.map(1);
        writeMarker(rec.marker());
        brefOut.writeByte(SEQ_CODED);
        for (int j=0, n=seq2Allele.size(); j<n; ++j) {
            brefOut.writeByte(seq2Allele.get(j));
        }
        bytesWritten += (seq2Allele.size() + 1) * Byte.BYTES;
    }

    private void writeAlleleCodedRec(RefGTRec rec) throws IOException {
        assert rec.isAlleleRecord();
        int nAlleles = rec.marker().nAlleles();
        int nullRow = rec.nullRow();
        writeMarker(rec.marker());
        brefOut.writeByte(ALLELE_CODED);
        bytesWritten += Byte.BYTES;
        for (int a=0; a<nAlleles; ++a) {
            if (a == nullRow) {
                brefOut.writeInt(-1);
                bytesWritten += Integer.BYTES;
            }
            else {
                int alCnt = rec.alleleCount(a);
                brefOut.writeInt(rec.alleleCount(a));
                for (int c=0; c<alCnt; ++c) {
                    brefOut.writeInt(rec.nonNullRowHap(a, c));
                }
                bytesWritten += (alCnt + 1)*Integer.BYTES;
            }
        }
    }

    private void writeMarker(Marker marker) throws IOException {
        String[] ids = MarkerUtils.ids(marker);
        int nIds = Math.min(ids.length, 255);
        brefOut.writeInt(marker.pos());
        brefOut.writeByte(nIds);
        bytesWritten += (Integer.BYTES + Byte.BYTES);
        for (int j=0; j<nIds; ++j) {
            writeString(ids[j], brefOut);
        }
        String[] alleleList = MarkerUtils.alleles(marker);
        byte alleleCode = isSNV(alleleList) ? snvCode(alleleList) : -1;
        brefOut.writeByte(alleleCode);
        bytesWritten += Byte.BYTES;
        if (alleleCode == -1) {
            writeStringArray(alleleList, brefOut);
            brefOut.writeInt(extractEnd(marker));
            bytesWritten += Integer.BYTES;
        }
    }

    private static int extractEnd(Marker marker) {
        String info = marker.info();
        int index = 4;  // start of base coordinate if info.startsWith("END=")==true
        if (info.startsWith("END=")==false) {
            index = info.indexOf(";END=") + 5;
            if (index<5) {  // ";END=" not found
                return -1;
            }
        }
        int endIndex = info.indexOf(Const.semicolon, index);
        if (endIndex == -1) {
            endIndex = info.length();
        }
        return Integer.parseInt(info.substring(index, endIndex));
    }

    private byte snvCode(String[] alleles) {
        int x = Arrays.binarySearch(SNV_PERMS, alleles, ALLELES_COMP);
        if (x < 0) {
            x = (-x - 1);
        }
        int code = (x << 2) + (alleles.length - 1);
        return (byte) code;
    }

    private boolean isSNV(String[] alleleList) {
        for (int j=0; j<alleleList.length; ++j) {
            if (BASES_SET.contains(alleleList[j])==false) {
                return false;
            }
        }
        return true;
    }

    private static Comparator<String[]> allelesComparator() {
        return (String[] o1, String[] o2) -> {
            int n = Math.min(o1.length, o2.length);
            for (int k=0; k<n; ++k) {
                char c1 = o1[k].charAt(0);
                char c2 = o2[k].charAt(0);
                if (c1 != c2) {
                    return (c1 < c2) ? -1 : 1;
                }
            }
            if (o1.length != o2.length) {
                return o1.length < o2.length ? -1 : 1;
            }
            else {
                return 0;
            }
        };
    }

    private Set<String> basesSet() {
        Set<String> set = new HashSet<>(4);
        set.add("A");
        set.add("C");
        set.add("G");
        set.add("T");
        return Collections.unmodifiableSet(set);
    }

    private void writeIndex() throws IOException {
        writeIndexChroms();
        int lastChrIndex = -1;
        for (int j=0, n=index.size(); j<n; ++j) {
            BrefBlock bb = index.get(j);
            long offset = bb.offset();
            int ci = bb.chromIndex();
            if (ci!=lastChrIndex) {
                lastChrIndex = ci;
                offset = -offset;
            }
            brefOut.writeLong(offset);
            brefOut.writeInt(bb.pos());
        }
        brefOut.writeLong(END_OF_INDEX);
        bytesWritten += (Long.BYTES + Integer.BYTES)*index.size() + Long.BYTES;
    }

    private void writeIndexChroms() throws IOException {
        int lastChrIndex = -1;
        List<String> chromList = new ArrayList<>();
        IntList firstChromBlock = new IntList();
        for (int j=0, n=index.size(); j<n; ++j) {
            BrefBlock bb = index.get(j);
            int chromIndex = bb.chromIndex();
            if (chromIndex!=lastChrIndex) {
                chromList.add(ChromIds.instance().id(chromIndex));
                firstChromBlock.add(j);
                lastChrIndex = chromIndex;
            }
        }
        Set<String> set = new HashSet<>(chromList);
        if (chromList.size() != set.size()) {
            Utilities.exit(CONTIGUITY_ERR);
        }
        String[] chromIds = chromList.toArray(new String[0]);
        int[] firstBlocks = firstChromBlock.toArray();
        writeStringArray(chromIds, brefOut);
        for (int b : firstBlocks) {
            brefOut.writeInt(b);
        }
    }

    private DataOutputStream dataOutputStream(File file) {
        OutputStream dos;
        if (file==null) {
            dos = new DataOutputStream(System.out);
        }
        else {
            dos = FileUtil.bufferedOutputStream(file);
        }
        return new DataOutputStream(dos);
    }

    private void writeStringArray(String[] sa, DataOutputStream dos)
            throws IOException {
        dos.writeInt(sa.length);
        bytesWritten += Integer.BYTES;
        for (String s : sa) {
            writeString(s, dos);
        }
    }

    private void writeString(String s, DataOutputStream dos)
           throws IOException {
        baos.reset();
        utf8Buffer.writeUTF(s);
        bytesWritten += baos.size();
        baos.writeTo(dos);
    }
}

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

import blbutil.BlockLineReader;
import blbutil.FileIt;
import blbutil.Utilities;
import blbutil.VcfFileIt;
import bref.SeqCoder3;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.Function;
import java.util.function.Predicate;

/**
 * <p>Class {@code RefIt} represents  an iterator whose {@code next()}
 * method returns an object storing data from a VCF record with
 * phased, non-missing genotypes.
 * </p>
 * <p>Instances of class {@code RefIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefIt implements VcfFileIt<RefGTRec> {

    /**
     * The default number of {@code GTRec} objects that are
     * stored in a buffer.
     */
    private static final int DEFAULT_BUFFER_SIZE = 1<<10;

    private final VcfHeader vcfHeader;
    private final Function<String, RefGTRec> mapper;
    private final Predicate<Marker> markerFilter;

    private final BlockLineReader reader;
    private final List<RefGTRec> lowFreqBuffer;
    private final Deque<RefGTRec> recBuffer;
    private final SeqCoder3 seqCoder;
    private final int maxSeqCodedAlleles;
    private final int maxSeqCodingMajorCnt;

    private int lastChrom = -1;

    /**
     * Create and returns a new {@code RefIt} instance from the specified
     * iterator.
     * @param it an iterator that returns lines of a VCF file
     * @return a new {@code RefIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code it}
     * @throws NullPointerException if {@code it == null}
     */
    public static RefIt create(FileIt<String> it) {
        return RefIt.create(it, FilterUtil.acceptAllPredicate(),
                FilterUtil.acceptAllPredicate(), false, DEFAULT_BUFFER_SIZE);
    }

    /**
     * Create and returns a new {@code RefIt} instance from the specified
     * objects.
     * @param it an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @return a new {@code RefIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code it}
     * @throws IllegalArgumentException if {@code bufferSize < 1}
     * @throws NullPointerException if {@code it == null}
     */
    public static RefIt create(FileIt<String> it, Predicate<String> sampleFilter,
            Predicate<Marker> markerFilter) {
        return new RefIt(it, sampleFilter, markerFilter, false,
                DEFAULT_BUFFER_SIZE);
    }

    /**
     * Create and returns a new {@code RefIt} instance from the specified
     * objects.
     * @param it an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param stripOptionalFields {@code true} if the VCF record's ID,
     * QUAL, FILTER, and INFO subfields should be discarded
     * @param bufferSize the number of VCF records stored in a buffer
     * @return a new {@code RefIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code it}
     * @throws IllegalArgumentException if {@code bufferSize < 1}
     * @throws NullPointerException if {@code it == null}
     */
    public static RefIt create(FileIt<String> it, Predicate<String> sampleFilter,
            Predicate<Marker> markerFilter, boolean stripOptionalFields,
            int bufferSize) {
        return new RefIt(it, sampleFilter, markerFilter, stripOptionalFields,
                bufferSize);
    }

    private RefIt(FileIt<String> it, Predicate<String> sampleFilter,
            Predicate<Marker> markerFilter, boolean stripOptionalFields, int blockSize) {
        if (blockSize < 1) {
            throw new IllegalArgumentException(String.valueOf(blockSize));
        }
        if (markerFilter==null) {
            markerFilter = FilterUtil.acceptAllPredicate();
        }
        String[] vcfHeaderLines = VcfHeader.readVcfHeader(it).toArray(String[]::new);
        String firstDataLine = it.hasNext() ? it.next() : null;
        if (firstDataLine==null) {
            throw new IllegalArgumentException("Missing VCF data lines (" + it.source() + ")");
        }
        boolean[] isDiploid = VcfHeader.isDiploid(firstDataLine);
        this.vcfHeader = new VcfHeader(it.source(), vcfHeaderLines, isDiploid, sampleFilter);
        this.mapper = (String s) -> {
            return RefGTRec.alleleRefGTRec(new VcfRecGTParser(vcfHeader, s, stripOptionalFields));
        } ;
        this.markerFilter = markerFilter;
        this.seqCoder = new SeqCoder3(vcfHeader.samples());
        this.maxSeqCodedAlleles = Math.min(seqCoder.maxNSeq(), SeqCoder3.MAX_NALLELES);
        this.maxSeqCodingMajorCnt = maxSeqCodingMajorCnt(vcfHeader.samples());

        this.lowFreqBuffer = new ArrayList<>();
        this.recBuffer = new ArrayDeque<>(blockSize);

        int nBlocks = 1;
        this.reader = BlockLineReader.create(it, blockSize, nBlocks);
        fillRecBuffer(firstDataLine);
    }

    private int maxSeqCodingMajorCnt(Samples samples) {
        int nHaps = samples.size() << 1;
        return (int) Math.floor(nHaps*SeqCoder3.COMPRESS_FREQ_THRESHOLD - 1);
    }

    private void fillRecBuffer() {
        fillRecBuffer(null);
    }

    private void fillRecBuffer(String firstDataLine) {
        assert recBuffer.isEmpty();
        while (recBuffer.isEmpty()) {
            String[] lines = reader.next();
            if (firstDataLine!=null) {
                lines = combine(firstDataLine, lines);
                firstDataLine = null;
            }
            if (lines==BlockLineReader.SENTINAL) {
                flushCompressedRecords();
                return;
            }
            else {
                RefGTRec[] recs = parseLines(lines);
                for (int j=0; j<recs.length; ++j) {
                    RefGTRec rec = recs[j];
                    int chrom = rec.marker().chromIndex();
                    if (lastChrom == -1) {
                        lastChrom = chrom;
                    }
                    if (chrom!=lastChrom || lowFreqBuffer.size()==Integer.MAX_VALUE) {
                        flushCompressedRecords();
                        lastChrom = chrom;
                    }
                    if (applySeqCoding(rec)==false) {
                        lowFreqBuffer.add(rec);
                    }
                    else {
                        boolean success = seqCoder.add(rec);
                        if (success == false) {
                            flushCompressedRecords();
                            success = seqCoder.add(rec);
                            assert success;
                        }
                        lowFreqBuffer.add(null);
                    }
                }
            }
        }
    }

    private String[] combine(String firstDataLine, String[] lines) {
        String[] modLines = new String[lines.length + 1];
        modLines[0] = firstDataLine;
        System.arraycopy(lines, 0, modLines, 1, lines.length);
        return modLines;
    }

    private RefGTRec[] parseLines(String[] lines) {
        RefGTRec[] recs = null;
        try {
            recs = Arrays.stream(lines)
                    .parallel()
                    .map(mapper)
                    .filter(e -> markerFilter.test(e.marker()))
                    .toArray(RefGTRec[]::new);
        }
        catch (Throwable t) {
            Utilities.exit(t);
        }
        return recs;
    }

    private void flushCompressedRecords() {
        List<RefGTRec> list = seqCoder.getCompressedList();
        int index = 0;
        for (int j=0, n=lowFreqBuffer.size(); j<n; ++j) {
            GTRec rec = lowFreqBuffer.get(j);
            if (rec==null) {
                lowFreqBuffer.set(j, list.get(index++));
            }
        }
        recBuffer.addAll(lowFreqBuffer);
        lowFreqBuffer.clear();
    }

    @Override
    public void close() {
        reader.close();
        recBuffer.clear();
        lowFreqBuffer.clear();
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !recBuffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public RefGTRec next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        RefGTRec first = recBuffer.removeFirst();
        if (recBuffer.isEmpty()) {
            fillRecBuffer();
        }
        return first;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException(this.getClass().toString());
    }

    @Override
    public String source() {
        return reader.source();
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }

    @Override
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        sb.append(" : ");
        sb.append(reader.source());
        return sb.toString();
    }

    private boolean applySeqCoding(RefGTRec rec) {
        assert rec.isAlleleRecord();
        if (rec.marker().nAlleles() > maxSeqCodedAlleles) {
            return false;
        }
        int nHaps = rec.size();
        int nullRow = rec.nullRow();
        int nullRowCnt = nHaps;
        for (int a=0, n=rec.marker().nAlleles(); a<n; ++a) {
            if (a!=nullRow) {
                nullRowCnt -= rec.alleleCount(a);
            }
        }
        return nullRowCnt<=maxSeqCodingMajorCnt;
    }
}

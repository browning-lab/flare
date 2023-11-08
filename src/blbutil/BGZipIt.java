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
package blbutil;

import ints.IntList;
import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;

/**
 * <p>Class {@code BGZipIt} is a {@code blbutil.FileIt<String>} whose
 * {@code next()} method returns lines of a bgzip-compressed file.
 * </p>
 * <p>The GZIP file format specification is described
 * <a href="https://www.ietf.org/rfc/rfc1952.txt">RFC 1952</a>
 * and the BGZIP file format specification is described in the
 * <a href="https://samtools.github.io/hts-specs/SAMv1.pdf">
 * Sequence Alignment/Map Format Specification</a>
 * </p>
 * <p>Instances of class {@code BGZipIt} are not thread safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BGZipIt implements FileIt<String> {

    private static final byte CR = 0x0D;
    private static final byte LF = 0x0A;
    private static final byte[] EOF = new byte[0];

    private static final byte GZIP_ID1 = 31;
    private static final byte GZIP_ID2 = (byte) 139;
    private static final byte GZIP_CM = 8;
    private static final byte GZIP_FLG = (1 << 2); // only required set bit
    private static final byte GZIP_XLEN1 = 6;
    private static final byte GZIP_XLEN2 = 0;
    private static final byte BGZIP_SI1 = 66;
    private static final byte BGZIP_SI2 = 67;
    private static final byte BGZIP_SLEN1 = 2;
    private static final byte BGZIP_SLEN2 = 0;

    private final InputStream is;
    private final File source;
    private final int nBufferedBlocks;
    private byte[] leftOverBytes;
    private final ArrayDeque<String> lines;

    /**
     * Constructs a new {@code BGZipIt} instance from the specified data
     * @param is an input stream that reads from a gzip-compressed
     * VCF file
     * @param nBufferedBlocks the number of buffered gzip blocks
     * @throws IllegalArgumentException if {@code nBufferedBlocks < 1}
     * @throws NullPointerException if {@code is == null}
     */
    public BGZipIt(InputStream is, int nBufferedBlocks) {
        this(is, nBufferedBlocks, null);
    }

    /**
     * Constructs a new {@code BGZipIt} instance from the specified data
     * @param is an input stream that reads gzip-compressed
     * VCF data
     * @param nBufferedBlocks the number of buffered gzip blocks
     * @param source the gzip-compressed VCF file that is read
     * @throws IllegalArgumentException if {@code nBufferedBlocks < 1}
     * @throws NullPointerException if {@code is == null}
     */
    public BGZipIt(InputStream is, int nBufferedBlocks, File source) {
        if (nBufferedBlocks < 1) {
            throw new IllegalArgumentException(String.valueOf(nBufferedBlocks));
        }
        this.is = is;
        this.source = source;
        this.nBufferedBlocks = nBufferedBlocks;
        this.leftOverBytes = new byte[0];
        this.lines = new ArrayDeque<>();
        fillBuffer();
    }

    @Override
    public void close() {
        try {
            is.close();
        } catch (IOException ex) {
            Utilities.exit(ex);
        }
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return lines.isEmpty()==false;
    }

    /**
     * Returns the next line of the VCF file. End of line characters are
     * not included in the returned line.
     * @return the next line of the VCF file
     * @throws NoSuchElementException if the VCF file has no more lines
     */
    @Override
    public String next() {
        String s = lines.remove();
        if (lines.isEmpty()) {
            fillBuffer();
        }
        return s;
    }

    @Override
    public File file() {
        return source;
    }

    private void fillBuffer() {
        byte[][] blocks = readAndInflateBlocks(is, leftOverBytes, nBufferedBlocks);
        if (blocks.length>0) {
            int[] eolIndices = IntStream.range(0, blocks.length)
                    .parallel()
                    .flatMap(j -> eolIndices(j, blocks[j]))
                    .toArray();
            leftOverBytes = leftOverBytes(blocks, eolIndices);
            addToLines(blocks, eolIndices, lines);
        }
    }

    private static IntStream eolIndices(int index, byte[] bytes) {
        IntList il = new IntList();
        for (int j=0; j<bytes.length; ++j) {
            if (bytes[j]==LF) {
                il.add(index);
                il.add(j);
            }
        }
        return il.stream();
    }

    private static byte[] leftOverBytes(byte[][] blocks, int[] eolIndices) {
        if (blocks.length==0) {
            return new byte[0];
        }
        else {
            int lastBlock = blocks.length-1;
            int endIndex = blocks[lastBlock].length;
            if (eolIndices.length==0) {
                return merge(blocks, 0, 0, lastBlock, endIndex);
            }
            else {
                int startBlock = eolIndices[eolIndices.length-2];
                int startIndex = eolIndices[eolIndices.length-1] + 1;
                return merge(blocks, startBlock, startIndex, lastBlock, endIndex);
            }
        }
    }

    private static void addToLines(byte[][] blocks, int[] eolIndices,
            ArrayDeque<String> lines) {
        List<String> tmpList = IntStream.range(0, eolIndices.length)
                .parallel()
                .filter(j -> (j & 0b1)==0)
                .mapToObj(j -> toString(blocks, eolIndices, j))
                .collect(Collectors.toList());
        lines.addAll(tmpList);
    }

    private static String toString(byte[][] blocks, int[] eolIndices, int index) {
        int block = eolIndices[index];
        int endIndex = eolIndices[index + 1];
        byte[] merged;
        if (index==0) {
            merged = merge(blocks, 0, 0, block, endIndex);
        }
        else {
            int startBlock = eolIndices[index-2];
            int startIndex = eolIndices[index-1] + 1;
            merged = merge(blocks, startBlock, startIndex, block, endIndex);
        }
        int lengthM1 = merged.length-1;
        if (lengthM1>=0 && merged[lengthM1]==CR) {
            // Correct for CR LF line ending on Windows systems
            return new String(merged, 0, lengthM1, StandardCharsets.UTF_8);
        }
        else {
            return new String(merged, StandardCharsets.UTF_8);
        }
    }

    private static byte[] merge(byte[][] blocks, int startBlock, int startIndex,
            int lastBlock, int endIndex) {
        // merge correctly handles startIndex == blocks[startBlock].length
        if (lastBlock==startBlock) {
            return Arrays.copyOfRange(blocks[startBlock], startIndex, endIndex);
        }
        else {
            int size = 0;
            for (int j=startBlock; j<lastBlock; ++j) {
                size += blocks[j].length;
            }
            size -= startIndex;
            size += endIndex;
            byte[] merged = new byte[size];
            int len = (blocks[startBlock].length - startIndex);
            System.arraycopy(blocks[startBlock], startIndex, merged, 0, len);
            for (int j=(startBlock + 1); j<lastBlock; ++j) {
                System.arraycopy(blocks[j], 0, merged, len, blocks[j].length);
                len += blocks[j].length;
            }
            System.arraycopy(blocks[lastBlock], 0, merged, len, endIndex);
            assert merged.length == (len + endIndex);
            return merged;
        }
    }

    private static byte[] inflateBlock(byte[] ba) {
        ByteArrayOutputStream os = new ByteArrayOutputStream(ba.length);
        byte[] buffer = new byte[1<<13];
        try (ByteArrayInputStream bais = new ByteArrayInputStream(ba);
                GZIPInputStream gzis = new GZIPInputStream(bais)) {
            int bytesRead;
            while ((bytesRead = gzis.read(buffer)) != -1) {
                os.write(buffer, 0, bytesRead);
            }
        }
        catch (IOException e) {
            Utilities.exit(e);
        }
        return os.toByteArray();
    }

    private static byte[][] readAndInflateBlocks(InputStream is, byte[] initialBytes, int nBlocks) {
        ArrayList<byte[]> compressedBlocks = new ArrayList<>(nBlocks);
        for (int j=0; j<nBlocks; ++j) {
            byte[] ba = readCompressedBlock(is);
            if (ba.length>0) {
                compressedBlocks.add(ba);
            }
            else if (ba==EOF) {
                break;
            }
        }
        byte[][] blocks = compressedBlocks.stream()
                .parallel()
                .map(ba -> inflateBlock(ba))
                .toArray(byte[][]::new);
        if (initialBytes.length>0) {
            int newLength = initialBytes.length + blocks[0].length;
            byte[] prependedBlock = Arrays.copyOf(initialBytes, newLength);
            System.arraycopy(blocks[0], 0, prependedBlock, initialBytes.length,
                    blocks[0].length);
            blocks[0] = prependedBlock;
        }
        return blocks;
    }

    private static byte[] readCompressedBlock(InputStream is) {
        byte[] ba = new byte[18];
        try {
            int bytesRead = 0;
            int offset = 0;
            while (offset<ba.length
                    && (bytesRead = is.read(ba, offset, ba.length - offset)) != -1) {
                offset += bytesRead;
            }
            if (offset==0) {
                return EOF;
            }
            if (offset==ba.length && isStartOfBgzipBlock(ba)) {
                int blockSize = ((ba[16] & 0xff) | ((ba[17] & 0xff) << 8)) + 1;
                ba = Arrays.copyOf(ba, blockSize);
                while (offset<ba.length
                        && (bytesRead = is.read(ba, offset, ba.length - offset)) != -1) {
                    offset += bytesRead;
                }
                if (offset < ba.length) {
                    Utilities.exit("Premature end of BGZIP block");
                }
            }
            else {
                Utilities.exit("Invalid BGZIP block header");
            }
        }
        catch (IOException e) {
            Utilities.exit(e);
        }
        return ba;
    }

    /**
     * Returns {@code true} if the first 16 bytes of the specified input stream
     * are a gzip header that includes a 6 byte extra field containing
     * the block size as described in the bgzip specification, and returns
     * {@code false} otherwise. The method sets a mark before reading
     * the initial bytes from the stream, and resets the stream to the
     * mark position before returning.
     * @param bis a buffered input stream
     * @return {@code true} if the first 16 bytes of the specified input stream
     * are a gzip header that includes a 6 byte extra field containing
     * the block size as described in the bgzip specification
     */
    public static boolean beginsWithBgzipBlock(BufferedInputStream bis) {
        assert bis.markSupported();
        int maxBytes = 16;
        int bytesRead = 0;
        int offset = 0;
        byte[] ba = new byte[maxBytes];
        bis.mark(maxBytes);
        try {
            while (offset<ba.length
                    && (bytesRead = bis.read(ba, offset, ba.length - offset)) != -1) {
                offset += bytesRead;
            }
            bis.reset();
        }
        catch(IOException ex) {
            Utilities.exit(ex);
        }
        return offset==ba.length && isStartOfBgzipBlock(ba);
    }

    private static boolean isStartOfBgzipBlock(byte[] buffer) {
    // isStartOfBgzipBlock() returns false if additional non-bgzip
    // subfields are present
        return (buffer.length >= 16
                && buffer[0] == GZIP_ID1)
                && (buffer[1] == GZIP_ID2)
                && (buffer[2] == GZIP_CM)
                && ((buffer[3] & GZIP_FLG)!=0)
                && (buffer[10] == GZIP_XLEN1)
                && (buffer[11] == GZIP_XLEN2)
                && (buffer[12] == BGZIP_SI1)
                && (buffer[13] == BGZIP_SI2)
                && (buffer[14] == BGZIP_SLEN1)
                && (buffer[15] == BGZIP_SLEN2);
    }
}

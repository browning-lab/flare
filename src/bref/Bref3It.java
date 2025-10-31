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

import blbutil.FileUtil;
import blbutil.Utilities;
import blbutil.VcfFileIt;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.NoSuchElementException;
import java.util.function.Predicate;
import vcf.FilterUtil;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.Samples;
import vcf.VcfHeader;

/**
 * <p>Class {@code Bref3It} represents  an iterator whose {@code next()} which
 * returns records from a bref version 3 file.
 * </p>
 * <p>Instances of class {@code Bref3It} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Bref3It implements VcfFileIt<RefGTRec> {

    private final String source;
    private final DataInputStream dataIn;
    private final Bref3Reader bref3Reader;
    private final Deque<RefGTRec> buffer;
    private final VcfHeader vcfHeader;

    /**
     * Constructs a new {@code Bref3It} instance.
     * @param brefFile a bref version 3 file or {@code null} if the
     * bref version 3 file is to be read from standard input
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of the specified bref file
     */
    public Bref3It(File brefFile) {
        this(brefFile, FilterUtil.acceptAllPredicate(), FilterUtil.acceptAllPredicate());
    }

    /**
     * Constructs a new {@code Bref4It} instance.
     * @param brefFile a bref v3 file or {@code null} if the bref3 file
     * is to be read from stdin
     * @param sampleFilter a sample filter
     * @param markerFilter a marker filter
     *
     * @throws IllegalArgumentException if a format error is detected in
     * the specified bref v3 file
     * @throws NullPointerException if
     * {@code (sampleFilter == null) || (markerFilter == null)}
     */
    public Bref3It(File brefFile, Predicate<String> sampleFilter,
            Predicate<Marker> markerFilter) {
        InputStream dis;
        if (brefFile==null) {
            dis = new BufferedInputStream(System.in);
        }
        else {
            dis = FileUtil.bufferedInputStream(brefFile);
        }
        this.source = brefFile==null ? "stdin" : brefFile.toString();
        this.dataIn = new DataInputStream(dis);
        this.bref3Reader = new Bref3Reader(brefFile, dataIn, sampleFilter, markerFilter);
        this.buffer = new ArrayDeque<>(500);
        bref3Reader.readBlock(dataIn, buffer);
        this.vcfHeader = vcfHeader(bref3Reader.samples());
    }

    private VcfHeader vcfHeader(Samples samples) {
        String[] sampleIds = samples.ids();
        boolean[] isDiploid = new boolean[sampleIds.length];
        Arrays.fill(isDiploid, true);
        String[] lines = new String[2];
        lines[0] = "##fileformat=bref3";
        StringBuilder sb = new StringBuilder(80 + 8*sampleIds.length);
        sb.append(VcfHeader.HEADER_PREFIX);
        for (String id : sampleIds) {
            sb.append('\t');
            sb.append(id);
        }
        lines[1] = sb.toString();
        return new VcfHeader(source, lines, isDiploid);
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !buffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public RefGTRec next() {
        if (hasNext()==false) {
            throw new NoSuchElementException();
        }
        RefGTRec rec = buffer.removeFirst();
        if (buffer.isEmpty()) {
            bref3Reader.readBlock(dataIn, buffer);
        }
        return rec;
    }

    @Override
    public void close() {
        try {
            dataIn.close();
        } catch (IOException ex) {
            Utilities.exit(ex, "Error closing file");
        }
        buffer.clear();
    }

    @Override
    public String source() {
        return source;
    }

    @Override
    public Samples samples() {
        return bref3Reader.samples();
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
        sb.append(source);
        return sb.toString();
    }
}

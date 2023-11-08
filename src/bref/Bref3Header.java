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

import blbutil.Const;
import blbutil.Filter;
import blbutil.Utilities;
import java.io.DataInput;
import java.io.File;
import java.io.IOException;
import java.util.stream.IntStream;
import vcf.Samples;

/**
 * <p>Class {@code Bref3Header} represents the header of a bref3 file (binary
 * reference format version 3).</p>
 *
 * <p>Instances of class {@code Bref3Header} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref3Header {

    private final String program;
    private final String[] sampleIds;
    private final int[] filteredHapIndices;
    private final int[] invFilteredHapIndices;
    private final Samples samples;

    /**
     * Reads and constructs a {@code Bref3Header} instance from the specified
     * data input stream. The constructor reads an {@code int} bref3 magic
     * number, a string with the name of the encoding program, and an array
     * of string sample identifiers from the data input stream.  All strings
     * are coded in modified UTF-8 format.  The Java Virtual
     * Machine will print an error message and exit if an I/O error or a
     * file format error is detected.
     *
     * @param source the source bref3 file or {@code null} if the bref3 file
     * is read from stdin
     * @param dataIn a data input stream
     * @param sampleFilter a sampleFilter
     * @throws NullPointerException if
     * {@code (dataIn == null) || (sampleFilter == null)}
     */
    public Bref3Header(File source, DataInput dataIn, Filter<String> sampleFilter) {
        String tmpProgramString = null;
        String[] tmpSampleIds = null;
        try {
            checkMagicNumber(dataIn.readInt());
            tmpProgramString = dataIn.readUTF();
            tmpSampleIds = Bref3Reader.readStringArray(dataIn);
        } catch (IOException ex) {
            Utilities.exit(ex, "Error reading file");
        }
        this.program = tmpProgramString;
        this.sampleIds = tmpSampleIds;
        int[] filteredSampleIndices = filteredSampleIndices(sampleIds, sampleFilter);
        if (filteredSampleIndices.length==0) {
            noSamplesError(source==null ? "stdin" : source.toString());
        }
        this.filteredHapIndices = hapIndices(filteredSampleIndices);
        this.invFilteredHapIndices = invArray(filteredHapIndices, (sampleIds.length<<1));
        this.samples = samples(sampleIds, filteredSampleIndices);
    }

    private static void checkMagicNumber(int magicNumber) throws IOException {
        if (magicNumber!=AsIsBref3Writer.MAGIC_NUMBER_V3) {
            String s = "ERROR: Unrecognized input file.  Was input file created "
                    + Const.nl + "with a different version of the bref program?";
            Utilities.exit(s);
        }
    }

    private static int[] filteredSampleIndices(String[] ids, Filter<String> filter) {
        return IntStream.range(0, ids.length)
                .filter(j -> filter.accept(ids[j]))
                .toArray();
    }

    private static void noSamplesError(String source) {
        String err = "All samples in the bref3 file have been excluded";
        String info = Const.nl + "Error      :  " + err
                + Const.nl     + "Bref3 file :  " + source;
        Utilities.exit(new Throwable(err), info);
    }

    /**
     * Returns the array of haplotype indices corresponding to the
     * specified array of sample indices. The returned array
     * will have length {@code (2 * sampIndices.length)}.  The
     * {@code (2 * k)}-th element of the returned array is
     * {@code (2 * sampIndices[k])}. The {@code ((2 * k) + 1)}-th element
     * of the returned array is {@code ((2 * sampIndices[k]) + 1)}.
     * @param sampIndices a list of sample indices
     * @return an array of haplotype indices
     *
     * @throws IllegalArgumentException if {@code sampIndices.length >= (1<<30)}
     * @throws NullPointerException if {@code sampIndices == null}
     */
    private static int[] hapIndices(int[] sampIndices) {
        if (sampIndices.length >= (1<<30)) {
            // 2*sampIndices.length will wrap around to a negative value
            throw new IllegalArgumentException(String.valueOf(sampIndices.length));
        }
        int[] hapIndices = new int[sampIndices.length<<1];
        for (int j=0, index=0; j<sampIndices.length; ++j) {
            int hap1 = sampIndices[j]<<1;
            hapIndices[index++] = hap1;
            hapIndices[index++] = (hap1 | 0b1);
        }
        return hapIndices;
    }

    /**
     * Returns the specified list of samples
     * @param sampleIds A list of sample identifiers
     * @param sampleIndices a list of {@code sampleIds} array indices
     * @return the specified list of samples
     * @throws IllegalArgumentException if any two elements of
     * {@code sampleIndices} are equal
     * @throws IndexOutOfBoundsException if there is a {@code j} such that
     * {@code (0 &le j) && (j < sampleIndices.length)} and
     * {@code (sampleIndices[j] < 0) || (sampleIndices[j] >= sampleIds.length)}
     * @throws NullPointerException if
     * {@code (sampleIds = null) || (sampleIndices == null)}
     */
    private static Samples samples(String[] sampleIds, int[] sampleIndices) {
        String[] ids = new String[sampleIndices.length];
        boolean[] isDiploid = new boolean[sampleIndices.length];
        for (int j=0; j<ids.length; ++j) {
            ids[j] = sampleIds[sampleIndices[j]];
            isDiploid[j] = true;
        }
        return Samples.fromIds(ids, isDiploid);
    }

    /**
     * Returns an array with the specified length that maps
     * {@code hapIndices[j]} to {@code j}. All indices of the returned array
     * that are not {@code hapIndices} array values are assigned the value
     * {@code -1}.
     * @param size the length of the returned array
     * @param hapIndices an array with nonnegative values
     * @return an inverse array
     *
     * @throws IllegalArgumentException if {@code size < 0}
     * @throws IllegalArgumentException if
     * {@code (hapIndices[j] == hapIndices[k])} for some
     * {@code ((0 &le j) && (j < k) && (k < hapIndices.length))}
     * @throws IndexOutOfBoundsException if
     * {@code (hapIndices[j] < 0) || (hapIndices[j] &ge; size)}
     * @throws NullPointerException if {@code hapIndices == null}
     */
    private static int[] invArray(int[] hapIndices, int size) {
        if (size < 0) {
            throw new IllegalArgumentException(String.valueOf(size));
        }
        int[] inverseArray = IntStream.range(0, size)
                .map(j -> -1)
                .toArray();
        for (int j=0; j<hapIndices.length; ++j) {
            if (inverseArray[hapIndices[j]] >= 0) {
                String s = "duplicate array value: " + hapIndices[j];
                throw new IllegalArgumentException(s);
            }
            inverseArray[hapIndices[j]] = j;
        }
        return inverseArray;
    }

    /**
     * Returns the program used to encode the bref3 file.
     *
     * @return the program used to encode the bref3 file
     */
    public String program() {
        return program;
    }

    /**
     * Returns the list of unfiltered sample identifiers.
     * @return the list of unfiltered sample identifiers
     */
    public String[] unfilteredSampleIds() {
        return sampleIds.clone();
    }

    /**
     * Returns the list of filtered haplotype indices.
     * @return the list of filtered haplotype indices
     */
    public int[] filteredHapIndices() {
        return filteredHapIndices.clone();
    }

    /**
     * Returns an array with length {@code (2 * this.unfilteredSampleIds())}
     * that maps {@code this.includedHapIndices()[j]} to {@code j}. All
     * other indices of the returned array that are assigned the value
     * {@code -1}.
     * @return an inverse of the {@code this.filteredHapIndices()} array
     */
    public int[] invfilteredHapIndices() {
        return invFilteredHapIndices.clone();
    }

    /**
     * Returns the list of filtered samples.
     * @return the list of filtered samples
     */
    public Samples samples() {
        return samples;
    }
}

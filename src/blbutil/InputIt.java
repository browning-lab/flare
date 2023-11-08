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

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

/**
 * <p>Class {@code InputIt} is a buffered iterator whose {@code next()}
 * method returns lines of a text input stream.
 * </p>
 * <p>If an {@code IOException} is thrown when an {@code InputIt}
 * instance reads from the text input stream, the {@code IOException}
 * is trapped, an error message is written to standard out, and the
 * Java Virtual Machine is terminated.
 * </p>
 * <p>The GZIP file format specification is described
 * <a href="https://www.ietf.org/rfc/rfc1952.txt">RFC 1952</a>
 * and the BGZIP file format specification is described in the
 * <a href="https://samtools.github.io/hts-specs/SAMv1.pdf">
 * Sequence Alignment/Map Format Specification</a>
 * </p>
 *
 * <p>Instances of class {@code InputIt} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class InputIt implements FileIt<String> {

    private final File file;
    private final BufferedReader in;
    private String next = null;

    /**
     * Constructs a new {@code InputStreamIterator} with default buffer
     * size that will iterate through lines of the specified input stream.
     *
     * @param is input stream of text data
     * @param file the file that is the source of the input stream
     */
    private InputIt(InputStream is, File file) {
        BufferedReader br = null;
        try {
            InputStreamReader isr = new InputStreamReader(is);
            br = new BufferedReader(isr);
            next = br.readLine();
        }
        catch(IOException e) {
            Utilities.exit(e, "Error reading " + file);
        }
        this.in = br;
        this.file = file;
    }

    /**
     * Constructs a new {@code InputStreamIterator} with default buffer size
     * that will iterate through the lines of the specified input stream.
     *
     * @param is input stream of text data
     * @param file the file that is the source of the input stream
     * @param bufferSize the buffer size in bytes
     *
     * @throws IllegalArgumentException if {@code bufferSize < 0}
     */
    private InputIt(InputStream is, File file, int bufferSize) {
        BufferedReader br = null;
        try {
            InputStreamReader isr = new InputStreamReader(is);
            br = new BufferedReader(isr, bufferSize);
            next = br.readLine();
        }
        catch(IOException e) {
            Utilities.exit(e, "Error reading " + file);
        }
        this.in = br;
        this.file = file;
    }

    @Override
    public File file() {
        return file;
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return (next != null);
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public String next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        String current = next;
        try {
            next = in.readLine();
        }
        catch (IOException e) {
            Utilities.exit(e, "Error reading " + file);
        }
        return current;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        String s = this.getClass().toString() + ".remove()";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public void close() {
        try {
            in.close();
        }
        catch (IOException e) {
            Utilities.exit(e, "Error closing " + in);
        }
        next=null;
    }

    /**
     * Returns a string representation of this iterator.  The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of this iterator
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(200);
        sb.append("[file= ");
        sb.append(file);
        sb.append("; next=\"");
        sb.append(next);
        sb.append("\"]");
        return sb.toString();
    }

    /**
     * Constructs and returns an {@code InputIt} instance with the default
     * buffer size that iterates through lines of text read from standard input.
     *
     * @return a new {@code InputIt} instance that iterates
     * through lines of text read from standard input
     */
    public static InputIt fromStdIn() {
        File file = null;
        return new InputIt(System.in, file);
    }

    /**
     * Constructs and returns a buffered {@code FileIt<String>} instance
     * that iterates through lines of the specified compressed or
     * uncompressed text file. If the filename ends in ".gz" or ".bgz", the
     * file must be GZIP-compressed.
     * @param file a compressed or uncompressed text file
     * @param nBufferedBlocks the number buffered GZIP blocks if the
     * specified file is bgzip-compressed
     * @return {@code FileIt<String>} instance that iterates through the
     * lines of the specified file
     *
     * @throws NullPointerException if {@code file == null}
     */
    public static FileIt<String> fromBGZipFile(File file, int nBufferedBlocks) {
        String filename = file.getName();
        try {
            InputStream is = new FileInputStream(file);
            BufferedInputStream bis = new BufferedInputStream(is);
            if (filename.endsWith(".gz") || filename.endsWith(".bgz")) {
                if (BGZipIt.beginsWithBgzipBlock(bis)) {
                    return new BGZipIt(bis, nBufferedBlocks, file);
                }
                else {
                    return new InputIt(new GZIPInputStream(bis), file);
                }
            }
            else {
                return new InputIt(is, file);
            }
        }
        catch(FileNotFoundException e) {
            Utilities.exit(e, "Error opening " + file);
        }
        catch(IOException e) {
            Utilities.exit(e, "Error reading " + file);
        }
        assert false;
        return null;
    }

    /**
     * Constructs and returns a buffered {@code InputIt} instance that
     * iterates through lines of the specified compressed or uncompressed
     * text file. If the filename ends in ".gz", the file must be
     * tGZIP-compressed.
     *
     * @param file a compressed or uncompressed text file
     * @return  a buffered {@code InputIt} instance that iterates
     * through lines of the specified text file
     * @throws NullPointerException if {@code file == null}
     */
    public static InputIt fromGzipFile(File file) {
        String filename = file.getName();
        try {
            InputStream is = new FileInputStream(file);
            if (filename.endsWith(".gz") || filename.endsWith(".bgz")) {
                return new InputIt(new GZIPInputStream(is), file);
            }
            else {
                return new InputIt(is, file);
            }
        }
        catch(FileNotFoundException e) {
            Utilities.exit(e, "Error opening " + file);
        }
        catch(IOException e) {
            Utilities.exit(e, "Error reading " + file);
        }
        assert false;
        return null;
    }

    /**
     * Constructs and returns a buffered {@code InputIt} instance
     * that iterates through lines of the specified text file.
     *
     * @param file a text file
     * @return a buffered {@code InputIt} instance that iterates through
     * lines of the specified text file
     * @throws NullPointerException if {@code filename == null}
     */
    public static InputIt fromTextFile(File file) {
        try {
            return new InputIt(new FileInputStream(file), file);
        }
        catch(FileNotFoundException e) {
            Utilities.exit(e, "Error opening " + file);
        }
        assert false;
        return null;
    }
}

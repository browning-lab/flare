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
package ints;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;

/**
 * <p>Class {@code UnsignedByteIndexArray} represents an immutable
 * array of integer values between 0 and 255 inclusive that is stored
 * as a {@code byte[]} array.  A {@code 0xff} mask is used to convert
 * {@code int} values to {@code byte values}.
 * </p>
 * <p>
 * Instances of {@code UnsignedByteArray} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class UnsignedByteArray implements IntArray {

    private final byte[] ba;

    /**
     * Constructs a new {@code UnsignedByteArray} instance from
     * the specified data.
     * @param ba an array of bytes which are interpreted as unsigned byte
     * values between 0 and 255
     * @throws NullPointerException if {@code ba == null}
     */
    public UnsignedByteArray(byte[] ba) {
        this.ba = ba.clone();
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the
     * specified data.
     * @param ba an array of bytes which are interpreted as unsigned byte
     * values between 0 and 255
     * @param from the first element to be included (inclusive)
     * @param to the last element to be included (exclusive)
     * @throws IndexOutOfBoundsException if {@code (from < 0 || to > ia.length)}
     * @throws NegativeArraySizeException if {@code to > from}
     * @throws NullPointerException if {@code ba == null}
     */
    public UnsignedByteArray(byte[] ba, int from, int to) {
        this.ba = Arrays.copyOfRange(ba, from, to);
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from
     * the specified data.
     * @param ia an array of integer values between 0 and 255 inclusive
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > 255)} for any index {@code j}
     * satisfying {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code ia == null}
     */
    public UnsignedByteArray(int[] ia) {
        this(ia, 0, ia.length);
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the specified data.
     * @param il an list of integer values between 0 and 255 inclusive
     * @throws IllegalArgumentException if
     * {@code (il.get(j) < 0 || il.get(j) > 255)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < il.size())}
     * @throws NullPointerException if {@code il == null}
     */
    public UnsignedByteArray(IntList il) {
        this(il, 0, il.size());
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the
     * specified data.
     * @param ia an array of integer values between 0 and 255 inclusive
     * @param from the first element to be included (inclusive)
     * @param to the last element to be included (exclusive)
     * @throws IllegalArgumentException if
     * {@code ((ia[j] < 0) || (ia[j] > 255))} for any index {@code j}
     * satisfying {@code (j >= from && j < to)}
     * @throws IndexOutOfBoundsException if {@code ((from < 0) || (to > ia.length))}
     * @throws NegativeArraySizeException if {@code (to > from)}
     * @throws NullPointerException if {@code (ia == null)}
     */
    public UnsignedByteArray(int[] ia, int from, int to) {
        this.ba = new byte[to - from];
        for (int j=from; j<to; ++j) {
            if (ia[j] < 0 || ia[j] > 255) {
                throw new IllegalArgumentException(String.valueOf(ia[j]));
            }
            ba[j - from] = (byte) ia[j];
        }
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the
     * specified data.
     * @param il an list of integer values between 0 and 255 inclusive
     * @param from the first element to be included (inclusive)
     * @param to the last element to be included (exclusive)
     * @throws IllegalArgumentException if
     * {@code ((il.get(j) < 0) || (il.get(j) > 255))} for any index {@code j}
     * satisfying  {@code ((j >= from) && (j < to))}
     * @throws IndexOutOfBoundsException if
     * {@code ((from < 0) || (to > il.size()))}
     * @throws NegativeArraySizeException if {@code (to > from)}
     * @throws NullPointerException if {@code (il == null)}
     */
    public UnsignedByteArray(IntList il, int from, int to) {
        this.ba = new byte[to - from];
        for (int j=from; j<to; ++j) {
            int value = il.get(j);
            if (value < 0 || value > 255) {
                throw new IllegalArgumentException(String.valueOf(value));
            }
            ba[j - from] = (byte) value;
        }
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the
     * specified data.
     * @param baos a byte array output stream whose elements are interpreted
     * as unsigned byte values between 0 and 255
     * @throws NullPointerException if {@code baos == null}
     */
    public UnsignedByteArray(ByteArrayOutputStream baos) {
        this.ba = baos.toByteArray();
    }

    @Override
    public int size() {
        return ba.length;
    }

    @Override
    public int get(int index) {
        return ba[index] & 0xff;
    }

    @Override
    public String toString() {
        return IntArray.asString(this);
    }

    /**
     * Write {@code this} byte array to the specified output stream.
     * @param os an output stream
     * @throws IOException if an I/O error occurs
     */
    public void write(OutputStream os) throws IOException {
        os.write(ba);
    }
}

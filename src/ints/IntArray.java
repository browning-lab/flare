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
package ints;

import java.util.Arrays;

/**
 * <p>Interface {@code IntArray} represents an immutable {@code int[]} array.
 * </p>
 * Instances of class {@code IntArray} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface IntArray {

    /**
     * Returns the number of elements in this {@code IntArray}.
     * @return the number of elements in this {@code IntArray}
     */
    int size();

    /**
     * Returns the specified array element.
     * @param index an array index
     * @return the specified array element
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    int get(int index);

    /**
     * Returns a copy of the specified array.
     * @param ia a list of integers
     * @return a copy of the specified array
     * @throws NullPointerException if {@code ia == null}
     */
    static int[] toArray(IntArray ia) {
        int[] copy = new int[ia.size()];
        for (int j=0; j<copy.length; ++j) {
            copy[j] = ia.get(j);
        }
        return copy;
    }

    /**
     * Returns a string representation of this {@code IntArray} by applying
     * {@code java.utils.Arrays.toString()} to an equivalent {@code int[]}
     * object.
     *
     * @param ia a list of integers
     * @return a string representation of this {@code IntArray}.
     * @throws NullPointerException if {@code ia == null}
     */
    static String asString(IntArray ia) {
        return Arrays.toString(toArray(ia));
    }

    /**
     * Returns {@code true} if the specified {@code IntArray} objects
     * represent the same sequence of integer values, and returns {@code false}
     * otherwise.
     * @param a a sequence of integer values
     * @param b a sequence of integer values
     * @return {@code true} if the specified {@code IntArray} objects
     * represent the same sequence of integer values
     */
    static boolean equals(IntArray a, IntArray b) {
        if (a==b) {
            return true;
        }
        else if (a.size()!=b.size()) {
            return false;
        }
        else {
            for (int j=0, n=a.size(); j<n; ++j) {
                if (a.get(j)!=b.get(j)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Returns the maximum element, or {@code Integer.MIN_VALUE} if
     * {@code this.size() == 0}.
     * @param ia a list of integers
     * @return the maximum element
     * @throws NullPointerException if {@code ia == null}
     */
    static int max(IntArray ia) {
        int max = Integer.MIN_VALUE;
        for (int j=0, n=ia.size(); j<n; ++j) {
            if (ia.get(j) > max) {
                max = ia.get(j);
            }
        }
        return max;
    }

    /**
     * Returns the minimum element, or {@code Integer.MAX_VALUE} if
     * {@code this.size() == 0}.
     * @param ia a list of integers
     * @return the minimum element
     * @throws NullPointerException if {@code ia == null}
     */
    static int min(IntArray ia) {
        int min = Integer.MAX_VALUE;
        for (int j=0, n=ia.size(); j<n; ++j) {
            if (ia.get(j) < min) {
                min = ia.get(j);
            }
        }
        return min;
    }

    /**
     * Returns a new {@code IntArray} instance that has the same
     * sequence of integers as the specified array of non-negative integers.
     * @param ia the array of non-negative integers to be copied
     * @param valueSize the exclusive end of the range of
     * array values
     * @return a new {@code IntArray} instance that has
     * the same sequence of integers as the specified array
     * @throws IllegalArgumentException if {@code valueSize < 1}
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code ia == null}
     */
    static IntArray packedCreate(int[] ia, int valueSize) {
        if (valueSize < 1) {
            throw new IllegalArgumentException(String.valueOf(valueSize));
        }
        if (valueSize<=16) {
            return new PackedIntArray(ia, valueSize);
        }
        else if (valueSize<=256) {
            return new UnsignedByteArray(ia, valueSize);
        }
        else if (valueSize<=65536) {
            return new CharArray(ia, valueSize);
        }
        else {
            return new WrappedIntArray(ia, valueSize);
        }
    }

    /**
     * Returns a new {@code IntArray} instance that has the same
     * sequence of integers as the specified list of non-negative integers.
     * @param il the list of non-negative integers to be copied
     * @param valueSize the exclusive end of the range of list values
     * @return a new {@code IntArray} instance that has
     * the same sequence of integers as the specified list
     * @throws IllegalArgumentException if {@code valueSize < 1}
     * @throws IllegalArgumentException if
     * {@code (il.get(j) < 0 || il.get(j)> valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < il.size())}
     * @throws NullPointerException if {@code il == null}
     */
    static IntArray packedCreate(IntList il, int valueSize) {
        if (valueSize < 1) {
            throw new IllegalArgumentException(String.valueOf(valueSize));
        }
        if (valueSize<=16) {
            return new PackedIntArray(il, valueSize);
        }
        else if (valueSize<=256) {
            return new UnsignedByteArray(il, valueSize);
        }
        else if (valueSize<=65536) {
            return new CharArray(il, valueSize);
        }
        else {
            return new WrappedIntArray(il, valueSize);
        }
    }

    /**
     * Returns a new {@code IntArray} instance that has the same
     * sequence of integers as the specified array of non-negative integers.
     * Each integer of the returned object is stored in 1, 2, or 4 bytes.
     * @param ia the array of non-negative integers to be copied
     * @param valueSize the exclusive end of the range of
     * array values
     * @return a new {@code IntArray} instance that has
     * the same sequence of integers as the specified array
     * @throws IllegalArgumentException if {@code valueSize < 1}
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code ia == null}
     */
    static IntArray create(int[] ia, int valueSize) {
        if (valueSize < 1) {
            throw new IllegalArgumentException(String.valueOf(valueSize));
        }
        if (valueSize<=256) {
            return new UnsignedByteArray(ia, valueSize);
        }
        else if (valueSize<=65536) {
            return new CharArray(ia, valueSize);
        }
        else {
            return new WrappedIntArray(ia, valueSize);
        }
    }

    /**
     * Returns a new {@code IntArray} instance that has the same
     * sequence of integers as the specified list of non-negative integers.
     * Each integer of the returned object is stored in 1, 2, or 4 bytes.
     * @param il the list of non-negative integers to be copied
     * @param valueSize the exclusive end of the range of list values
     * @return a new {@code IntArray} instance that has
     * the same sequence of integers as the specified list
     * @throws IllegalArgumentException if {@code valueSize < 1}
     * @throws IllegalArgumentException if
     * {@code (il.get(j) < 0 || il.get(j)> valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < il.size())}
     * @throws NullPointerException if {@code il == null}
     */
    static IntArray create(IntList il, int valueSize) {
        if (valueSize < 1) {
            throw new IllegalArgumentException(String.valueOf(valueSize));
        }
        if (valueSize<=256) {
            return new UnsignedByteArray(il, valueSize);
        }
        else if (valueSize<=65536) {
            return new CharArray(il, valueSize);
        }
        else {
            return new WrappedIntArray(il, valueSize);
        }
    }
}

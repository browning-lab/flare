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

/**
 * <p>Class {@code IndexArray} stores an array whose entries
 * are elements of a bounded set of non-negative integers along with an upper
 * bound. Both the stored array and upper bound are specified at the time of
 * object construction. It is the client code that constructs a
 * {@code IndexArray} object's responsibility to ensure that the bound
 * specified at construction is correct. The contract for this class
 * is undefined if the bound is specified at object construction is incorrect.
 * It is recommended that the bound be the minimum integer that is greater
 * than all elements in the stored array, but not is not a requirement.
 * </p>
 * <p>Instances of class {@code IndexArray} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IndexArray implements IntArray {

    private final IntArray intArray;
    private final int valueSize;

    /**
     * Constructs a new {@code IndexArray} instance from the specified
     * data. The contract for this constructor and class is unspecified if
     * the sequence indices are not a subset of
     * {@code 0, 1, ..., (valueSize - 1)}.
     *
     * @param intArray an array with non-negative values
     * @param valueSize a value that is greater than all elements
     * of the specified array
     * @throws NullPointerException if {@code IndexArray == null}
     */
    public IndexArray(int[] intArray, int valueSize) {
        if (intArray==null) {
            throw new NullPointerException(IntArray.class.toString());
        }
        this.intArray = IntArray.packedCreate(intArray, valueSize);
        this.valueSize = valueSize;
    }

    /**
     * Constructs a new {@code IndexArray} instance from the specified
     * data. The contract for this constructor and class is unspecified if
     * the sequence indices are not a subset of
     * {@code 0, 1, ..., (valueSize - 1)}.
     *
     * @param intArray an array with non-negative values
     * @param valueSize a value that is greater than all elements
     * of the specified array
     * @throws NullPointerException if {@code IndexArray == null}
     */
    public IndexArray(IntArray intArray, int valueSize) {
        if (intArray==null) {
            throw new NullPointerException(IntArray.class.toString());
        }
        this.intArray = intArray;
        this.valueSize = valueSize;
    }

    @Override
    public int get(int index) {
        return intArray.get(index);
    }

    @Override
    public int size() {
        return intArray.size();
    }

    /**
     * Returns the value size that was specified at construction.
     * @return the value size that was specified at construction
     */
    public int valueSize() {
        return valueSize;
    }

    /**
     * Returns the wrapped {@code IntArray} object.
     * @return the wrapped {@code IntArray} object
     */
    public IntArray intArray() {
        return intArray;
    }

    /**
     * Returns the minimum integer that is greater than all elements in the
     * specified array of non-negative values.
     * @param ia an array of non-negative values
     * @return the minimum integer that is greater than all elements in the
     * specified array of non-negative values
     * @throws IllegalArgumentException if any element of the specified array is
     * negative
     */
    public static int valueSize(int[] ia) {
        int max = -1;
        for (int i : ia) {
            if (i<0) {
                throw new IllegalArgumentException(String.valueOf(i));
            }
            if (i > max) {
                max = i;
            }
        }
        return (max + 1);
    }

    /**
     * Returns the minimum integer that is greater than all elements in the
     * specified list of non-negative values.
     * @param ia an array of non-negative values
     * @return the minimum integer that is greater than all elements in the
     * specified array of non-negative values
     * @throws IllegalArgumentException if any element of the specified array
     * is negative
     */
    public static int valueSize(IntArray ia) {
        int max = -1;
        for (int j=0, n=ia.size(); j<n; ++j) {
            int i = ia.get(j);
            if (i<0) {
                throw new IllegalArgumentException(String.valueOf(i));
            }
            if (i > max) {
                max = i;
            }
        }
        return (max + 1);
    }
}

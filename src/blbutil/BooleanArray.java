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
package blbutil;

import java.util.Arrays;

/**
 * <p>Class {@code BooleanArray} represents an immutable array of
 * boolean values.</p>
 *
 * <p>Instances of {@code BooleanArray} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BooleanArray {

    private final boolean[] ba;

    /**
     * Returns the estimated number of bytes consumed by this object,
     * excluding the overhead bytes required by {@code this}.
     * @return the estimated number of bytes required to store this object
     */
    public long estBytes() {
        int overhead = 12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 8;  // assume 8 bytes per object reference
        estBytes += 4 + ba.length; // 4 bytes to store array length
        return estBytes;
    }

    /**
     * Constructs a {@code BooleanArray} instance from the specified data.
     * @param ba the array to be copied and stored
     * @throws NullPointerException if {@code (ba == null)}
     */
    public BooleanArray(boolean[] ba) {
        this.ba = ba.clone();
    }

    /**
     * Returns the number of values in this {@code BooleanArray}.
     * @return the number of values in this {@code BooleanArray}
     */
    public int size() {
        return ba.length;
    }

    /**
     * Returns the specified {@code boolean} value.
     * @param index an array index
     * @return the specified {@code boolean} value
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public boolean get(int index) {
        return ba[index];
    }

    /**
     * Returns {@code java.util.Arrays.hashCode(this.toArray())}.
     * @return {@code java.util.Arrays.hashCode(this.toArray())}
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(ba);
    }

    /**
     * Returns {@code true} if this {@code BooleanArray} and the specified
     * object represent identical sequences of boolean values and returns
     * {@code false} otherwise.
     * @param obj the object to be compared with {@code this}
     * @return {@code true} if this {@code BooleanArray} and the specified
     * object represent identical sequences of boolean values
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final BooleanArray other = (BooleanArray) obj;
        return Arrays.equals(this.ba, other.ba);
    }

    /**
     * Returns the sequence of boolean values as a {@code boolean[]} array.
     * @return the sequence of boolean values as a {@code boolean[]} array
     */
    public boolean[] toArray() {
        return ba.clone();
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.toArray)}.
     *
     * @return {@code java.util.Arrays.toString(this.toArray)}
     */
    @Override
    public String toString() {
        return Arrays.toString(ba);
    }
}

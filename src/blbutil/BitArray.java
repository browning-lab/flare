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
 * <p>Interface {@code BitArray} represents a mutable sequence of bits
 * with a fixed length.</p>
 *
 * <p>Instances of {@code BitArray} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BitArray {

    private static final int LOG2_BITS_PER_WORD = 6;
    private static final long WORD_MASK = 0xffffffffffffffffL;
    private static final int BIT_INDEX_MASK = (1 << LOG2_BITS_PER_WORD) - 1;

    private final long[] words;
    private final int size;

    /**
     * Returns the estimated number of bytes consumed by this object,
     * excluding the overhead bytes required by {@code this}.
     * @return the estimated number of bytes required to store this object
     */
    public long estBytes() {
        int overhead = 12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 8;  // assume 8 bytes per reference
        estBytes += 4 + (8*words.length); // 4 bytes to store array length
        estBytes += 4;
        return estBytes;
    }

    /**
     * Constructs a {@code BitArray} instance with the specified {@code size}
     * and having all bits set to 0 (unset).
     * @param size the number of bits
     * @throws IllegalArgumentException if {@code size < 0}
     */
    public BitArray(int size) {
        if (size<0) {
            throw new IllegalArgumentException(String.valueOf(size));
        }
        int nWords = (size + Long.SIZE - 1) / Long.SIZE;
        this.words = new long[nWords];
        this.size = size;
    }

    /**
     * Constructs a {@code BitArray} instance from the specified values.
     * @param values a sequence of bits
     * @param size the number of bits
     * @throws IllegalArgumentException if {@code size < 0}
     * @throws IllegalArgumentException if
     * {@code values.length != (size + Long.SIZE - 1) / Long.SIZE}
     * @throws NullPointerException if {@code values == null}
     */
    public BitArray(long[] values, int size) {
        if (size<0) {
            throw new IllegalArgumentException(String.valueOf(size));
        }
        int nWords = (size + Long.SIZE - 1) / Long.SIZE;
        if (values.length!=nWords) {
            throw new IllegalArgumentException(String.valueOf(values.length));
        }
        this.words = Arrays.copyOf(values, nWords);
        this.size = size;
    }

    /**
     * Constructs a new {@code BitArray} instance with the same sequence of
     * bits and the same size as the specified {@code BitArray}.
     * @param bitList a sequence of bits to be copied
     * @throws NullPointerException if {@code bitList == null}
     */
    public BitArray(BitArray bitList) {
        this.words = bitList.words.clone();
        this.size = bitList.size;
    }

    /**
     * Returns the number of bits in this {@code BitArray}.
     * @return the number of bits in this {@code BitArray}
     */
    public int size() {
        return size;
    }

    /**
     * Returns the specified bit as a {@code boolean} value. A 1 bit returns
     * {@code true} and a 0 bit returns {@code false}.
     * @param index a bit index
     * @return the specified bit as a {@code boolean} value.
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public boolean get(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int wordIndex =  index >> LOG2_BITS_PER_WORD;
        return (words[wordIndex] & (1L << index)) != 0L;
    }

    /**
     * Returns the specified bit as a {@code int} value.
     * @param index a bit index
     * @return the specified bit as a {@code int} value.
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public int getAsInt(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int wordIndex =  index >> LOG2_BITS_PER_WORD;
        return  (int) ((words[wordIndex] & (1L << index)) >> index);
    }

    /**
     * Sets the specified bit.
     * @param index a bit index
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public void set(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int wordIndex =  index >> LOG2_BITS_PER_WORD;
        words[wordIndex] |= (1L << index);
    }

    /**
     * Clears the specified bit.
     * @param index a bit index
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public void clear(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int wordIndex =  index >> LOG2_BITS_PER_WORD;
        words[wordIndex] &= ~(1L << index);
    }

    /**
     * Clears all bits.
     */
    public void clear() {
        Arrays.fill(words, 0L);
    }

    /**
     * Returns a new {@code BitArray} of size {@code (from - to)} that is a
     * copy of the specified bit indices of this {@code BitArray}.
     * @param from the first bit to be copied (inclusive)
     * @param to the last bit to be copied (exclusive)
     * @return a new {@code BitArray} of size {@code (from - to)} that is a
     * copy of the specified bit indices of this {@code BitArray}
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || to > this.size}
     */
    public BitArray restrict(int from, int to) {
        if (from<0 || from>to || to>size) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return new BitArray(0);
        }
        BitArray result = new BitArray(to-from);
        int nResultWordsM1 = result.words.length - 1;
        final boolean isWordAligned = ((from & BIT_INDEX_MASK) == 0);
        int srcWord = from >> LOG2_BITS_PER_WORD;
        for (int w=0; w<nResultWordsM1; w++, srcWord++) {
            result.words[w] = isWordAligned ? words[srcWord]
                    : (words[srcWord] >>> from) | (words[srcWord+1] << -from);
        }
        long endWordMask = WORD_MASK >>> -to;
        result.words[nResultWordsM1] =
                ((to-1) & BIT_INDEX_MASK) < (from & BIT_INDEX_MASK)
                ? ((words[srcWord] >>> from) | (words[srcWord+1] & endWordMask) << -from)
                : (words[srcWord] & endWordMask) >>> from;

        return result;
    }

    /**
     * Replaced the specified bits in this {@code Bitlist} with the corresponding
     * bits in the specified {@code BitArray}.
     * @param src the {@code BitArray} to be copied from
     * @param from the first bit to be copied (inclusive)
     * @param to the last bit to be copied (exclusive)
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || to > this.size || to > src.size()}
     * @throws NullPointerException if {@code src == null}
     */
    public void copyFrom(BitArray src, int from, int to) {
        if (from<0 || from>to || to>size || to>src.size()) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return;
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            words[startWord] ^= ((words[startWord] ^ src.words[startWord]) & mask);
        }
        else {
            words[startWord] ^= ((words[startWord] ^ src.words[startWord]) & startWordMask);
            for (int j=startWord+1; j<endWord; ++j) {
                words[j] = src.words[j];
            }
            words[endWord] ^= ((words[endWord] ^ src.words[endWord]) & endWordMask);
        }
    }

    /**
     * Returns a hash code for the specified bits in this {@code Bitlist}
     * @param from the first bit (inclusive)
     * @param to the last bit (exclusive)
     * @return a hash code for the specified bits in this {@code Bitlist}
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || to > this.size}
     */
    public int hash(int from, int to) {
        if (from<0 || from>to || to>size) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return 0;
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            return longHashCode(words[startWord] & mask);
        }
        else {
            long longHash = (words[startWord] & startWordMask);
            for (int j=startWord+1; j<endWord; ++j) {
                longHash ^= words[j];
            }
            longHash ^= (words[endWord] & endWordMask);
            return longHashCode(longHash);
        }
    }

    public static int longHashCode(long value) {
        return (int)(value ^ (value >>> 32));
    }

    /**
     * Swaps the specified bits of the two specified {@code Bitlist} objects.
     * @param a the first {@code BitArray}
     * @param b the second {@code BitArray}
     * @param from the first bit to be copied (inclusive)
     * @param to the last bit to be copied (exclusive)
     * @throws IllegalArgumentException if
     * {@code s.size() != b.size()}
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || to > a.size()}
     * @throws NullPointerException if {@code a == null || b == null}
     */
    public static void swapBits(BitArray a, BitArray b, int from, int to) {
        if (a.size()!=b.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (from<0 || from>to || to>a.size()) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return;
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            a.words[startWord] ^= (b.words[startWord] & mask);
            b.words[startWord] ^= (a.words[startWord] & mask);
            a.words[startWord] ^= (b.words[startWord] & mask);
        }
        else {
            a.words[startWord] ^= (b.words[startWord] & startWordMask);
            b.words[startWord] ^= (a.words[startWord] & startWordMask);
            a.words[startWord] ^= (b.words[startWord] & startWordMask);
            for (int j=startWord+1; j<endWord; ++j) {
                a.words[j] ^= b.words[j];
                b.words[j] ^= a.words[j];
                a.words[j] ^= b.words[j];
            }
            a.words[endWord] ^= (b.words[endWord] & endWordMask);
            b.words[endWord] ^= (a.words[endWord] & endWordMask);
            a.words[endWord] ^= (b.words[endWord] & endWordMask);
        }
    }

    /**
     * Returns {@code true} if this {@code Bitlist} and the specified
     * {@code BitArray} have identical sequences of bits for the specified
     * indices, and returns {@code false} otherwise. Returns {@code true}
     * if {@code (from == to) && (0 <= from) && (from < other.size())}.
     * @param other the {@code BitArray} to be compared with {@code this} for
     * equality.
     * @param from the first bit to be compared (inclusive)
     * @param to the last bit to be compared (exclusive)
     * @return {@code true} if this {@code Bitlist} and the specified
     * {@code BitArray} have identical sequences of bits for the specified
     * indices.
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || to > this.size || to > other.size()}
     * @throws NullPointerException if {@code other == null}
     */
    public boolean equal(BitArray other, int from, int to) {
        if (from < 0 || from>to || to>size || to>other.size()) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return true;
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            return ((words[startWord] ^ other.words[startWord]) & mask) == 0L;
        }
        else {
            boolean areEqual = true;
            areEqual &= (((words[startWord] ^ other.words[startWord]) & startWordMask) == 0L);
            for (int j=startWord+1; j<endWord; ++j) {
                areEqual &= (words[j]==other.words[j]);
            }
            areEqual &= (((words[endWord] ^ other.words[endWord]) & endWordMask) == 0L);
            return areEqual;
        }
    }

    /**
     * Returns this {@code BitArray} as a {@code long} array.
     * @return this {@code BitArray} as a {@code long} array
     */
    public long[] toLongArray() {
        return words.clone();
    }

    /**
     * Returns a string representation of this {@code BitArray}.
     * The exact details of the representation are unspecified and subject
     * to change.
     *
     * @return a string representation of this {@code BitArray}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(size);
        for (int j=0; j<size; ++j) {
            sb.append(get(j) ? '1' : '0');
        }
        return sb.toString();
    }

    /**
     * Returns {@code true} if the specified {@code BitArray} objects
     * represent identical bit sequences having the same size,
     * and returns {@code false} otherwise.
     * @param a a sequence of long values
     * @param b a sequence of long values
     * @return {@code true} if the specified {@code BitArray} objects
     * represent identical bit sequences having the same size
     * @throws NullPointerException if {@code a == null || b == null}
     */
    public static boolean equals(BitArray a, BitArray b) {
        return a.size()==b.size() && Arrays.equals(a.words, b.words);
    }
}

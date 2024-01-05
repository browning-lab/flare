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

import java.util.Arrays;

/**
 * <p>Class {@code IntIntMap} represents a map with integer keys and integer
 * values.
 * </p>
 * <p>Class {@code IntIntMap} is not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IntIntMap {

    private static final int NIL = -1;
    private static final float LOAD_FACTOR = 0.75f;

    private int size;
    private int nBuckets;

    private int[] next;
    private int[] data; // stores list index of keys and values
    private int[] keys;
    private int[] values;
    private int firstFreeIndex;

    /**
     * Creates a new {@code IntMap} instance.
     *
     * @param capacity the initial capacity of this map
     * @throws IllegalArgumentException if
     * {@code capacity < 0 || (capacity > (1 << 30))}
     */
    public IntIntMap(int capacity) {
        if (capacity < 1 || capacity > (1<<30)) {
            throw new IllegalArgumentException(String.valueOf(capacity));
        }
        int numBuckets = (int) Math.ceil(capacity/LOAD_FACTOR) + 1;
        allocateArrays(capacity, numBuckets);
        initializeFields(numBuckets);
    }

    private void allocateArrays(int capacity, int numBuckets) {
        this.next = new int[numBuckets + capacity];
        this.data = new int[numBuckets + capacity];
        this.keys = new int[capacity];
        this.values = new int[capacity];
    }

    private void initializeFields(int numBuckets) {
        size = 0;
        nBuckets = numBuckets;
        firstFreeIndex = nBuckets;
        Arrays.fill(next, 0, nBuckets, NIL);
        for (int j=nBuckets; j<next.length; ++j) {
            next[j] = j+1;
        }
    }

    /*
     * Increases the capacity of the internal hash table.
     */
    private void rehash(int newCapacity) {
        if (newCapacity > size) {
            int oldSize = size;
            int[] oldKeys = keys.clone();
            int[] oldValues = values.clone();
            int newNumBuckets = (int) Math.ceil(newCapacity/LOAD_FACTOR);
            allocateArrays(newCapacity, newNumBuckets);
            initializeFields(newNumBuckets);
            for (int j=0; j<oldSize; ++j) {
                put(oldKeys[j], oldValues[j]);
            }
        }
    }

    /**
     * Removes all keys from this map.
     */
    public void clear() {
        initializeFields(nBuckets);
    }

    /**
     * Returns {@code true} if the map contains the specified key,
     * and returns {@code false} otherwise.
     * @param key a key
     * @return {@code true} if the map contains the specified key
     */
    public boolean contains(int key) {
        return indexOf(key)>=0;
    }

    /**
     * Returns the current index of the specified key, or returns {@code -1}
     * if the key is not found in this map.
     * @param key a key index
     * @return the current index of the specified key
     */
    private int indexOf(int key) {
        int index = next[bucket(key)];
        while (index!=NIL && keys[data[index]]<key) {
            index = next[index];
        }
        return (index!=NIL && keys[data[index]]==key) ? index : -1;
    }

    /**
     * Adds the specified key and value to this map.  Returns {@code true}
     * if this map was changed by the operation, and returns {@code false}
     * otherwise. The indexing of keys immediately before and after this
     * method is invoked may differ if this map was changed by the operation.
     * @param key they key
     * @param value the value
     * @return {@code true} if this map was changed by the operation
     */
    public boolean put(int key, int value) {
        int prevIndex = prevIndex(key);
        int nextIndex = next[prevIndex];
        if (nextIndex==NIL || keys[data[nextIndex]]!=key) {
            int index = firstFreeIndex;
            firstFreeIndex = next[firstFreeIndex];
            next[prevIndex] = index;
            data[index] = size;
            next[index] = nextIndex;
            keys[size] = key;
            values[size++] = value;
            if (size == keys.length) {
                int newCapacity = 3*keys.length/2 + 1;
                rehash(newCapacity);
            }
            return true;
        }
        else if (values[data[nextIndex]]!=value) {
            values[data[nextIndex]] = value;
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Removes the specified key from this map. Returns {@code true}
     * if this map was changed by the operation, and returns {@code false}
     * otherwise. The indexing of keys immediately before and after this
     * method is invoked may differ if this map was changed by the operation.
     *
     * @param key a key index
     * @return {@code true} if this map was changed by the operation
     */
    public boolean remove(int key) {
        int prevIndex = prevIndex(key);
        int index = next[prevIndex];
        if (index==NIL || keys[data[index]]!=key) {
            return false;
        }
        else {
            int oldListIndex = data[index];
            next[prevIndex] = next[index];
            next[index] = firstFreeIndex;
            firstFreeIndex = index;

            --size;
            if (oldListIndex!=size) {
                // overwrite removed key
                index = indexOf(keys[size]);
                data[index] = oldListIndex;
                keys[oldListIndex] = keys[size];
                values[oldListIndex] = values[size];
            }
            return true;
        }
    }

    private int bucket(int key) {
        return Math.abs((71*key) % nBuckets);
    }

    private int prevIndex(int key) {
        int prevIndex = bucket(key);
        int index = next[prevIndex];
        while (index!=NIL && keys[data[index]]<key) {
            prevIndex = index;
            index = next[index];
        }
        return prevIndex;
    }

    /**
     * Returns the specified key.
     * @param index an index of a key in this map
     * @return the specified key
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int key(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return keys[index];
    }

    /**
     * Returns the value for the specified key or the specified sentinel
     * value if the specified key is not present in this map.
     * @param key the key
     * @param sentinel the value to be returned if the specified key is not
     * present in this map
     * @return the specified value
     */
    public int get(int key, int sentinel) {
        int index = indexOf(key);
        if (index == -1) {
            return sentinel;
        }
        return values[data[index]];
    }

    /**
     * Returns the number of keys in this map.
     *
     * @return the number of keys in this map
     */
    public int size() {
        return size;
    }

    /**
     * Returns an array containing the keys in this map. The returned
     * array will satisfy:
     * {@code this.toArray()[j]==this.key(j)} for each
     * {@code j} satisfying {@code (0 <= j && j < this.size())}.
     * @return an array containing the keys in this map
     */
    public int[] keys() {
        return Arrays.copyOf(keys, size);
    }

    /**
     * Returns an array containing the values in this map. The returned
     * array will satisfy:
     * {@code this.toArray()[j]==this.value(this.key(j))} for each
     * {@code j} satisfying {@code (0 <= j && j < this.size())}.
     * @return an array containing the values in this map
     */
    public int[] values() {
        return Arrays.copyOf(values, size);
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.keys())}.
     *
     * @return {@code java.util.Arrays.toString(this.keys())}
     */
    @Override
    public String toString() {
        return Arrays.toString(keys());
    }

    // xxx code for testing class
    private static void main(String[] args) {
        IntIntMap map1 = new IntIntMap(4);
        java.util.Map<Integer, Integer> map2 = new java.util.HashMap<>(100);
        java.util.Random rand = new java.util.Random(0);
        int nTests = 5000;
        int maxKey = 100;
        for (int j=0; j<nTests; ++j) {
            int key = j % maxKey;
            int value = rand.nextInt();
            double d = rand.nextDouble();
            if (d < 0.005) {
                map1.clear();
                map2.clear();
            }
            else if (d < 0.4) {
                int i1 = map1.get(key, -1);
                Integer i2 = map2.get(key);
                assert (i1==-1 && i2==null) || (i1==i2);
            }
            else if (d < 0.7) {
                int i1 = map1.get(key, -1);
                map1.put(key, value);
                Integer i2 = map2.put(key, value);
                assert (i1==-1 && i2==null) || (i1==i2);
            }
            else {
                int i1 = map1.get(key, -1);
                map1.remove(key);
                Integer i2 = map2.remove(key);
                assert (i1==-1 && i2==null) || (i1==i2);
            }
            assert map1.size()==map2.size();
        }
    }
    // xx code for testing class
}

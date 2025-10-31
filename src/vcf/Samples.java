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
package vcf;

import java.util.Arrays;

/**
 * <p>Class {@code Samples} stores a list of samples.
 * </p>
 * Instances of class {@code Samples} are immutable.
 *
 * @author Brian L. Browning
 */
public final class Samples {

    private final String[] ids;
    private final boolean[] isDiploid;

   /**
     * Constructs a new {@code Samples} instance corresponding to the
     * specified list of sample identifiers.  A warning is printed to standard
     * error if any string occurs more than once in the {@code ids} array.
     * @param ids an array of sample identifiers
     * @param isDiploid a boolean array whose {@code k}-th value is {@code true}
     * if the {@code k}-th sample is diploid, and is {@code false} if the
     * {@code k}-th sample is haploid
     *
     * @throws IllegalArgumentException if {@code ids.length != isDiploid.length}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code ((0 <= j) && j < ids.length) && (ids[j].length()==0)}
     * @throws NullPointerException if {@code ids == null || isDiploid == null}
     * @throws NullPointerException if there exists {@code j} such that
     * {@code ((0 <= j) && j < ids.length) && (ids[j]==null))}
     */
    public Samples(String[] ids, boolean[] isDiploid) {
        if (ids.length!=isDiploid.length) {
            throw new IllegalArgumentException(String.valueOf(isDiploid));
        }
        checkForNullsAndDuplicates(ids);
        this.ids = ids.clone();
        this.isDiploid = isDiploid.clone();
    }

    private static void checkForNullsAndDuplicates(String[] ids) {
        String[] sortedCopy = Arrays.stream(ids)
                .parallel()
                .sorted()
                .toArray(String[]::new);
        if (sortedCopy.length>0 && sortedCopy[0].length()==0) {
            throw new IllegalArgumentException("Empty string identifier");
        }
        for (int j=1; j<sortedCopy.length; ++j) {
            if (sortedCopy[j].length()==0) {
                throw new IllegalArgumentException("Empty string identifier");
            }
            if (sortedCopy[j].equals(sortedCopy[j-1])) {
                System.err.println("Warning: duplicate sample identifier: "
                        + sortedCopy[j]);
            }
        }
    }

    /**
     * Returns a new samples instance by combining the two list of samples
     * in the specified order
     * @param first the first list of samples
     * @param second the second list of samples
     * @return the combined samples
     * @throws IllegalArgumentException if the two lists of samples are not
     * disjoint
     * @throws NullPointerException if
     * {@code first == null || second == null}
     */
    public static Samples combine(Samples first, Samples second) {
        int n1 = first.size();
        int n2 = second.size();
        int n = n1 + n2;
        String[] ids = new String[n];
        boolean[] isDiploid = new boolean[n];
        System.arraycopy(first.ids, 0, ids, 0, n1);
        System.arraycopy(second.ids, 0, ids, n1, n2);
        System.arraycopy(first.isDiploid, 0, isDiploid, 0, n1);
        System.arraycopy(second.isDiploid, 0, isDiploid, n1, n2);
        return new Samples(ids, isDiploid);
    }

    /**
     * Returns a hash code value for the object.
     * @return a hash code value for the object.
     */
    @Override
    public int hashCode() {
        int hash = 59;
        hash += 29*Arrays.hashCode(this.isDiploid);
        hash += 29*Arrays.hashCode(this.ids);
        return hash;
    }

    /**
     * Returns {@code true} if the specified object is a
     * {@code Samples} object which represents the same ordered
     * list of samples as {@code this}, and returns {@code false}
     * otherwise.
     * @param obj the object to be tested for equality with {@code this}
     * @return {@code true} if the specified object is a
     * {@code Samples} object which represents the same ordered
     * list of samples as {@code this}
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null || this.getClass() != obj.getClass()) {
            return false;
        }
        final Samples other = (Samples) obj;
        if (Arrays.equals(this.isDiploid, other.isDiploid)==false) {
            return false;
        }
        return Arrays.equals(this.ids, other.ids);
    }

    /**
     * Returns the number of samples in this list.
     * @return the number of samples in this list
     */
    public int size() {
        return ids.length;
    }

    /**
     * Returns the identifier for the sample with the specified
     * index in this list of samples.
     * @param index a sample index
     * @return the identifier for the sample with the specified
     * index in this list of samples
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public String id(int index) {
        return ids[index];
    }

    /**
     * Returns this list of samples as an array of sample identifiers.
     * The returned array has length {@code this.size()}, and it
     * satisfies {@code this.ids()[j].equals(this.id(j))} for
     * {@code 0 <= j && j < this.size()}
     * @return this list of samples as an array of sample identifiers
     */
    public String[] ids() {
        return ids.clone();
    }

     /**
      * Returns {@code true} if the specified sample has two alleles per
      * genotype, and returns {@code false} if the sample has one allele
      * per genotype.
      * @param sample a sample index
      * @return {@code true} if the specified sample is diploid
      * @throws IndexOutOfBoundsException if
      * {@code sample < 0 || sample >= this.size()}
      */
    public boolean isDiploid(int sample) {
        return isDiploid[sample];
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.ids())}.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return Arrays.toString(ids());
    }
}

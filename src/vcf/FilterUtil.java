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

import blbutil.Utilities;
import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Predicate;

/**
 * Class {@code FilterUtil} contains static methods for constructing
 * marker filters.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class FilterUtil {

    private FilterUtil() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns a predicate that accepts all non-null objects.
     * @param <E> the type of object that is filtered
     * @return a filter that accepts all non-null objects
     */
    public static <E> Predicate<E> acceptAllPredicate() {
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return true;
        };
    }

   /**
     * <p>Returns a predicate that accepts all strings that are not contained
     * in the specified file, and returns a predicate that accepts all
     * strings if {@code (file == null)}. Blank lines in the file are ignored,
     * and each line in the file will be trimmed to eliminate beginning
     * and ending white-space. The returned predicate will throw
     * a {@code NullPointerException} if its {@code test()} method is
     * invoked on a null reference</p>
     *
     * <p>If an {@code IOException} is thrown while accessing the specified
     * file, an error message will be printed to standard error and the
     * Java virtual machine will terminate.</p>
     *
     * @param file a text file with one identifier per line
     * @return a predicate that rejects all non-null objects that are not
     * contained in the specified file
     * @throws IllegalArgumentException if the specified file does not exist
     * @throws IllegalArgumentException if the specified file is a directory
     * @throws IllegalArgumentException if any line of the specified
     * file contains two non-white-space characters separated by one or
     * more white-space characters
     */
    public static Predicate<String> includePredicate(File file) {
        if (file==null) {
            return acceptAllPredicate();
        }
        final HashSet<String> includeSet = Utilities.idSet(file);
        return (String str) -> {
            if (str == null) {
                throw new NullPointerException("e==null");
            }
            return includeSet.contains(str);
        };
    }

    /**
     * <p>Returns a predicate that rejects all strings that are not contained
     * in the specified file, and returns a predicate that accepts all
     * strings if {@code (file == null)}. Blank lines in the file are ignored,
     * and each line in the file will be trimmed to eliminate beginning
     * and ending white-space. The returned predicate will throw
     * a {@code NullPointerException} if its {@code test()} method is
     * invoked on a null reference</p>
     *
     * <p>If an {@code IOException} is thrown while accessing the specified
     * file, an error message will be printed to standard error and the
     * Java virtual machine will terminate.</p>
     *
     * @param file a text file with one identifier per line
     * @return a predicate that rejects all non-null objects that are not
     * contained in the specified file
     * @throws IllegalArgumentException if the specified file does not exist
     * @throws IllegalArgumentException if the specified file is a directory
     * @throws IllegalArgumentException if any line of the specified
     * file contains two non-white-space characters separated by one or
     * more white-space characters
     */
    public static Predicate<String> excludePredicate(File file) {
        if (file==null) {
            return acceptAllPredicate();
        }
        final HashSet<String> excludeSet = Utilities.idSet(file);
        return (String str) -> {
            if (str == null) {
                throw new NullPointerException("e==null");
            }
            return excludeSet.contains(str)==false;
        };
    }

    /**
     * Returns a predicate that accepts all non-null objects that are
     * contained in the specified collection.
     * @param <E> the type of object that is filtered
     * @param include the collection of objects that will be accepted by
     * the filter
     * @return a predicate that accepts all non-null objects that are
     * contained in the specified collection
     * @throws NullPointerException if {@code include == null}
     */
    public static <E> Predicate<E> includePredicate(Collection<E> include) {
        final HashSet<E> includeSet = new HashSet<>(include);
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return includeSet.contains(e);
        };
    }

    /**
     * Returns a predicate that accepts all non-null objects that are not
     * contained in the specified collection.
     * @param <E> the type of object that is filtered
     * @param exclude the collection of objects that will be rejected
     * by the filter
     * @return a predicate that accepts all non-null objects that are not
     * contained in the specified collection
     * @throws NullPointerException if {@code exclude == null}
     */
    public static <E> Predicate<E> excludePredicate(Collection<E> exclude) {
        final HashSet<E> excludeSet = new HashSet<>(exclude);
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return excludeSet.contains(e)==false;
        };
    }

    /**
     * Returns a filter that accepts all non-null objects that are contained
     * in the specified {@code include} collection and that are not in the
     * contained in the specified {@code exclude} collection.
     * @param <E> the type of object that is filtered.
     * @param include the collection of elements that may be accepted
     * by the filter.
     * @param exclude the collection of elements that will be rejected
     * by the filter.
     * @return a filter that accepts all non-null objects that are contained
     * in the specified {@code include} collection and that are not in the
     * contained in the specified {@code exclude} collection.
     * @throws NullPointerException if {@code ((include==null) || (exclude==null))}
     */
    public static <E> Predicate<E> includeExcludePredicate(
            final Collection<E> include, final Collection<E> exclude) {
        final Set<E> includeSet = new HashSet<>(include);
        final Set<E> excludeSet = new HashSet<>(exclude);
        return new Predicate<E>() {
            @Override
            public boolean test(E e) {
                if (e==null) {
                    throw new NullPointerException("e==null");
                }
                return includeSet.contains(e) && !excludeSet.contains(e);
            }
        };
    }

    /**
     * Returns a filter that excludes markers that have an identifier
     * or genome coordinates that matches a line of the specified file,
     * or returns a filter that accepts all markers if the
     * {@code excludeMarkersFile} parameter is {@code null}.
     * Genome coordinates must be in "CHROM:POS" format.
     * @param excludeMarkersFile a file that contains an identifier
     * or genome coordinate of one excluded marker on each line
     * @return a filter that excludes markers that have an identifier
     * or genome coordinates that matches a line of the specified file,
     * or {@code null} if the {@code excludeMarkersFile} parameter is
     * {@code null}
     *
     * @throws IllegalArgumentException if the specified file does not exist
     * @throws IllegalArgumentException if the specified file is a directory
     * @throws IllegalArgumentException if any line of the specified
     * file contains two non-white-space characters separated by one or
     * more white-space characters
     */
    public static Predicate<Marker> markerFilter(File excludeMarkersFile) {
        Set<String> excludeIds;
        if (excludeMarkersFile==null) {
            return FilterUtil.acceptAllPredicate();
        }
        else {
            excludeIds = Utilities.idSet(excludeMarkersFile);
            return excludeMarkerFilter(excludeIds);
        }
    }

    /**
     * Returns {@code true} if the specified marker has an identifier
     * is in the specified set, or if ("marker.chromID()" + ":" + "marker.pos()")
     * is in the specified set, and returns {@code false} otherwise.
     * @param marker a marker
     * @param set a set of marker identifiers and chromosome positions in
     * "CHROM:POS" format
     * @return {@code true} if the specified marker has an identifier
     * is in the specified set or if ("marker.chromID()" + ":" + "marker.pos()")
     * is in the specified set
     * @throws NullPointerException if {@code marker == null || set == null}
     */
    public static boolean markerIsInSet(Marker marker, Set<String> set) {
        String[] ids = MarkerUtils.ids(marker);
        for (int j=0; j<ids.length; ++j) {
            if (set.contains(ids[j])) {
                return true;
            }
        }
        String posId = marker.chromID() + ':' + marker.pos();
        return set.contains(posId);
    }

    /**
     * Returns a filter that accepts all markers which do not have an
     * identifier or chromomsome position present in the specified
     * collection.
     * A marker is excluded if {@code exclude.contains(marker.id(j)) == true}
     * for any {@code 0 <= j < marker.nIds()} or if
     * {@code exclude.contains(marker.chromID() + ":" + marker.pos()) == true}.
     * @param exclude a collection of marker identifiers and chromosome
     * positions in "CHROM:POS" format
     * @return a filter that accepts all markers which do not have an
     * identifier or chromomsome position present in the specified
     * collection
     * @throws NullPointerException if {@code exclude == null}
     */
    public static Predicate<Marker> excludeMarkerFilter(Collection<String> exclude) {
        final Set<String> excludeSet = new HashSet<>(exclude);
        if (excludeSet.isEmpty()) {
            return marker -> true;
        }
        else {
            return (Marker marker) -> !markerIsInSet(marker, excludeSet);
        }
    }
}

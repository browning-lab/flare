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

import java.io.File;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;

/**
 * Class {@code Utilities} contains miscellaneous static utility methods.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Utilities {

    private Utilities() {
        // private constructor to prevent instantiation
    }

   /**
    * Returns a string representation of the command line arguments.
    * The exact details of the representation are unspecified and
    * subject to change.
    *
    * @param program the name of the program's jar file.
    * @param args command line arguments.
    * @return a string representation of the command line arguments.
    */
    public static String commandLine(String program, String[] args) {
        StringBuilder sb = new StringBuilder(args.length*20);
        long maxMemory = Runtime.getRuntime().maxMemory();
        sb.append(Const.nl);
        sb.append("Command line: java");
        if (maxMemory!=Long.MAX_VALUE) {
            long maxMb = maxMemory / (1024*1024);
            sb.append(" -Xmx");
            sb.append(maxMb);
            sb.append("m");
        }
        sb.append(" -jar ");
        sb.append(program);
        sb.append(Const.nl);
        for (int j = 0; j < args.length; ++j) {
            sb.append("  ");
            sb.append(args[j]);
            sb.append(Const.nl);
        }
        return sb.toString();
    }

    /**
     * Prints a summary of memory use at the time of method invocation
     * to standard output.
     * @param msg a string a message to be printed with the summary
     * of memory use
     */
    public static void printMemoryUse(String msg) {
        long Mb = 1024*1024;
        Runtime rt = Runtime.getRuntime();
        System.out.println(Const.nl + msg
                + Const.tab + "maxMb=" + (rt.maxMemory()/Mb)
                + Const.tab + "totalMb=" + (rt.totalMemory()/Mb)
                + Const.tab + "freeMb=" + (rt.freeMemory()/Mb)
                + Const.tab + "usedMb=" + ((rt.totalMemory() - rt.freeMemory())/Mb));
    }

    /**
     * Returns the current local time as a string.  The
     * exact details of the string representation
     * are unspecified and subject to change.
     *
     * @return the current local time as a string.
     */
    public static String timeStamp() {
        Date now = new Date();
        SimpleDateFormat sdf =
                new SimpleDateFormat("hh:mm a z 'on' dd MMM yyyy");
        return sdf.format(now);
    }

    /**
     * Returns the current minutes and seconds as a string.  The
     * exact details of the string representation
     * are unspecified and subject to change.
     *
     * @return the current local time as a string.
     */
    public static String minutesAndSeconds() {
        Date now = new Date();
        SimpleDateFormat sdf =
                new SimpleDateFormat("mm:ss");
        return sdf.format(now);
    }

    /**
     * <p>Returns a set of identifiers found in a text file that has
     * one identifier per line.  The empty set is returned if
     * {@code file == null}. Blank lines are ignored, and white-space that
     * begins or ends a line is ignored.
     * </p>
     * If an {@code IOException} is thrown, an error message is printed
     * to standard error and the Java virtual machine is forced to terminate.
     *
     * @param file a text file with one identifier per line
     * @return a set of identifiers
     *
     * @throws IllegalArgumentException if the specified file does not exist
     * @throws IllegalArgumentException if the specified file is a directory
     * @throws IllegalArgumentException if any line of the specified
     * file contains two non-white-space characters separated by one or
     * more white-space characters
     */
    public static HashSet<String> idSet(File file) {
        if (file==null) {
            return new HashSet<>();
        }
        else {
            if (file.exists()==false) {
                String s = "file does not exist: " + file;
                throw new IllegalArgumentException(s);
            }
            if (file.isDirectory()) {
                String s = "file is a directory: " + file;
                throw new IllegalArgumentException(s);
            }
            HashSet<String> idSet = new HashSet<>();
            try (FileIt<String> it = InputIt.fromGzipFile(file)) {
                while (it.hasNext()) {
                    String line = it.next().trim();
                    if (line.length() > 0) {
                        if (StringUtil.countFields(line) > 1) {
                            String s = "Line has more than one white-space delimited field (file: '"
                                    + file + "'; line: '" + line + "')";
                            throw new IllegalArgumentException(s);
                        }
                        idSet.add(line);
                    }
                }
            }
            return idSet;
        }
    }

    /**
     * Prints the specified string to the specified {@code PrintWriter} and
     * to standard out.  The line separator string is not appended to the
     * specified string before printing.
     *
     * @param out a print writer
     * @param s a string to be printed
     *
     * @throws NullPointerException if {@code out == null}
     */
    public static void duoPrint(PrintWriter out, CharSequence s) {
        System.out.print(s);
        out.print(s);
    }

   /**
     * Prints the specified string to the specified {@code PrintWriter} and
     * to standard out.  The line separator string is appended to the
     * specified string before printing.
     *
     * @param out a print writer
     * @param s a string to be printed
     *
     * @throws NullPointerException if {@code out == null}
     */
    public static void duoPrintln(PrintWriter out, CharSequence s) {
        System.out.println(s);
        out.println(s);
    }

     /**
     * Returns a string representation of the specified elapsed time
     * in the format "H hours M minutes S seconds".
     *
     * @param nanoseconds the elapsed time in nanoseconds
     *
     * @return a string representation of the specified elapsed time
     */
    public static String elapsedNanos(long nanoseconds) {
        long seconds = Math.round(nanoseconds /1000000000.0);
        StringBuilder sb = new StringBuilder(80);
        if (seconds >= 3600) {
            long hours = seconds / 3600;
            sb.append(hours);
            sb.append(hours==1 ? " hour " : " hours ");
            seconds %= 3600;

        }
        if (seconds >= 60) {
            long minutes = seconds / 60;
            sb.append(minutes);
            sb.append(minutes==1 ? " minute " : " minutes ");
            seconds %= 60;
        }
        sb.append(seconds);
        sb.append(seconds==1 ? " second" : " seconds");
        return sb.toString();
    }

    /**
     * Prints the specified exception, its stack trace, and
     * the specified string to standard error and then terminates the
     * Java virtual machine.
     *
     * @param str a string to be printed to standard error
     * @param thr an exception or error to be printed to standard error
     *
     * @throws NullPointerException if {@code t == null}
     */
    public static void exit(Throwable thr, String str) {
        System.err.println();
        System.err.println(str);
        System.err.println();
        thr.printStackTrace(System.err);
        exit();
    }

    /**
     * Prints the specified exception, its stack trace and then terminates the
     * Java virtual machine.
     *
     * @param t an exception or error to be printed to standard error
     *
     * @throws NullPointerException if {@code e == null}
     */
    public static void exit(Throwable t) {
        t.printStackTrace(System.err);
        exit();
    }

    /**
     * Prints the specified string to standard error and then terminates the
     * Java virtual machine.
     *
     * @param s a string to be written to standard error
     */
    public static void exit(String s) {
        System.err.println();
        System.err.print(s);
        exit();
    }

    /**
     * Prints "Terminating program." to standard error and terminates
     * the Java virtual machine.
     */
    public static void exit() {
        System.err.println();
        System.err.println("Terminating program.");
        System.err.flush();
        System.exit(1);
    }

    /**
     * Randomly shuffles the elements of the specified array.
     * @param ia an array to be shuffled
     * @param random a random number generator
     * @throws NullPointerException if {@code ia == null || random == null}
     */
    public static void shuffle(int[] ia, Random random) {
        for (int j=0; j<ia.length; ++j) {
            int x = random.nextInt(ia.length-j);
            int tmp = ia[j];
            ia[j] = ia[j+x];
            ia[j+x] = tmp;
        }
    }

    /**
     * Shuffles the specified array so that a random set of ${code nElements}
     * elements from the array are the first ${code nElements} elements.
     * The specified array is unchanged if {@code nElements <= 0}.
     * @param ia an array to be shuffled
     * @param nElements the size of the random set of elements which
     * are shuffled to the beginning of the array
     * @param random a random number generator
     * @throws IndexOutOfBoundsException if {@code nElements > ia.length}
     * @throws NullPointerException if {@code ia == null || random == null}
     */
    public static void shuffle(int[] ia, int nElements, Random random) {
        for (int j=0; j<nElements; ++j) {
            int x = random.nextInt(ia.length-j);
            int tmp = ia[j];
            ia[j] = ia[j+x];
            ia[j+x] = tmp;
        }
    }

    /**
     * Returns a two-dimensional array with indices of common array
     * elements. If there are {@code n} elements in common between
     * the {@code a1} and {@code a2} arrays, then for {@code j} satisfying
     * {@code (0 <= j && j < n)}, the {@code [0][j]} and {@code P[1][j]}
     * elements of the returned array are the indices of the {@code j}-th
     * common element in the {@code a1} and {@code a2} arrays respectively.
     * Common elements are indexed according to their order in the
     * {@code a1} array.
     * @param <E> the type of array element
     * @param a1 an array of elements
     * @param a2 an array of elements
     * @return a two-dimensional array with indicies of common array
     * elements
     * @throws IllegalArgumentException if the {@code a1} or {@code a2} arrays
     * have a duplicate element
     * @throws NullPointerException if {@code a1 == null || a2 == null}
     */
    public static <E> int[][] commonIndices(E[] a1, E[] a2) {
        Map<E, Integer> map1 = arrayToMap(a1);
        Map<E, Integer> map2 = arrayToMap(a2);
        map1.keySet().retainAll(map2.keySet());
        int[][] indices = new int[2][map1.size()];
        int index = 0;
        for (int j=0; j<a1.length; ++j) {
            E id = a1[j];
            if (map1.containsKey(id)) {
                assert map2.keySet().contains(id);
                indices[0][index] = map1.get(id);
                indices[1][index++] = map2.get(id);
            }
        }
        return indices;
    }

    /**
     * Returns a map from array element to array index
     * @param <E> the type of array element
     * @param array an array
     * @return a map from array element to array index
     * @throws IllegalArgumentException if the {@code array}
     * has a duplicate element
     * @throws NullPointerException if {@code array == null}
     */
    public static <E> HashMap<E, Integer> arrayToMap(E[] array) {
        HashMap<E, Integer> map = new HashMap<>();
        for (int j=0; j<array.length; ++j) {
            E e = array[j];
            if (map.containsKey(e)) {
                String s = "duplicate element in array: " + e;
                throw new IllegalArgumentException(s);
            }
            map.put(e, j);
        }
        return map;
    }
}


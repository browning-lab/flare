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
package admix;

import blbutil.Const;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.FloatArray;
import blbutil.InputIt;
import blbutil.StringUtil;
import blbutil.Utilities;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.IntStream;
import vcf.Samples;
import vcf.VcfHeader;

/**
 * <p>Class {@code AdmixUtils} contains static utility methods for the
 * admix package.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixUtils {

    private AdmixUtils() {
        // private constructor to prevent instantiation
    }

    /**
     * Deletes any existing file with the specified filename and returns the
     * a {@code File} corresponding to the specified filename.
     * @param filename a filename
     * @return a {@code File} corresponding to the specified filename
     * @throws NullPointerException if {@code filename == null}
     */
    public static File clobberAndReturnFile(String filename) {
        File file = new File(filename);
        if (file.exists()) {
            file.delete();
        }
        return file;
    }

    /**
     * Returns the list of unique elements in the order they appear in the
     * specified array.
     *
     * @param sa a list of strings
     * @return the list of unique elements in the order they appear in the
     * specified array
     * @throws NullPointerException if {@code (sa == null)} or if
     * {@code (0 <= j && j < sa.length && sa[j] == null)}
    */
    public static String[] uniqueElements(String[] sa) {
        ArrayList<String> idList = new ArrayList<>();
        HashSet<String> idSet = new HashSet<>();
        for (String panel : sa) {
            if (idSet.add(panel)) {
                idList.add(panel);
            }
        }
        return idList.toArray(new String[0]);
    }

    /**
     * Returns the sum of the elements in the specified array.
     * @param da the elements to be summed
     * @return the sum of the elements in the specified array
     * @throws NullPointerException if {@code da == null}
     */
    public static double sum(double[] da) {
        double sum = 0.0;
        for (double d : da) {
            sum += d;
        }
        return sum;
    }

    /**
     * Sets all elements of the specified two-dimensional array to the
     * specified value.
     * @param da a two-dimensional array
     * @param value a value
     * @throws NullPointerException if {@code da == null}
     * @throws NullPointerException if {@code (da[j] == null)}
     * for any {@code j} such that {@code (0 <= j && j < fa.length)}
     */
    public static void fill(double[][] da, double value) {
        for (double[] row : da) {
            Arrays.fill(row, value);
        }
    }

    /**
     * Copies the specified initial columns of the two-dimensional source array
     * to the two-dimensional destination array.  The first dimension of each
     * array is the row index and the second dimension is the column index.
     * @param src a two-dimensional source array
     * @param dst a two-dimensional destination array
     * @param nCols the number of initial columns to be copied
     * @throws IndexOutOfBoundsException if {@code dst.length < src.length}
     * @throws IndexOutOfBoundsException if
     * {@code (src[j].length < nCols || dst[j].length < nCols}
     * for some {@code j} such that {@code (0 <= j && j < src.length)}
     * @throws IndexOutOfBoundsException if {@code nCols < 0}
     * @throws NullPointerException if {@code src == null || dst == null}
     * @throws NullPointerException if
     * {@code (src[j] == null || dst[j] == null)} for any {@code j} such
     * that {@code (0 <= j && j < src.length)}
     */
    public static void copy(double[][] src, double[][] dst, int nCols) {
        for (int j=0; j<src.length; ++j) {
            System.arraycopy(src[j], 0, dst[j], 0, nCols);
        }
    }

    /**
     * Replaces each element of the specified two-dimensional array with the
     * product of the element and the specified value.
     * @param da a two-dimensional array
     * @param scaleFactor the scale factor
     * @throws NullPointerException if {@code da == null}
     * @throws NullPointerException if {@code (da[j] == null)}
     * for any {@code j} such that {@code (0 <= j && j < da.length)}
     */
    public static void scale(double[][] da, double scaleFactor) {
        for (double[] row : da) {
            scale(row, scaleFactor);
        }
    }

    /**
     * Replaces each element of the specified array with the
     * product of the element and the specified value.
     * @param da an array
     * @param scaleFactor the scale factor
     * @throws NullPointerException if {@code da == null}
     */
    public static void scale(double[] da, double scaleFactor) {
        for (int j=0; j<da.length; ++j) {
            da[j] *= scaleFactor;
        }
    }

    /**
     * Replaces each of the specified elements of the specified array with the
     * product of the element and the specified value.
     * @param da an array
     * @param start the start index (inclusive)
     * @param end the end index (exclusive)
     * @param scaleFactor the scale factor
     * @throws IllegalArgumentException if {@code end < start}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end > fa.length}
     * @throws NullPointerException if {@code da == null}
     */
    public static void scale(double[] da, int start, int end, double scaleFactor) {
        if (end<start) {
            throw new IllegalArgumentException(String.valueOf(end));
        }
        for (int j=start; j<end; ++j) {
            da[j] *= scaleFactor;
        }
    }
    /**
     * Returns an array whose {@code k}-th element is the probability of
     * transitioning to a random HMM state in the genomic interval between
     * the {@code (k-1)}-st marker and the {@code k}-th marker. The 0-th
     * entry of the returned array is {@code 0.0}.
     * @param recombIntensity the intensity of the exponential distribution
     * that determines the probability of transitioning to a random HMM state
     * @param genDist a list whose {@code k}-th element is the cM
     * distance between the {@code (k-1)}-st marker and the {@code k}-th
     * marker, or {@code 0.0} if {@code (k == 0)}
     * @return the probability of transitioning to a random HMM state in each
     * intermarker interval.
     * @throws IllegalArgumentException if
     * {@code recombIntensity <= 0.0 || Float.isFinite(recombIntensity) == false}
     * @throws NullPointerException if {@code genDist == null}
     */
    public static double[] pRec(double recombIntensity, FloatArray genDist) {
        if (recombIntensity <= 0.0 || Double.isFinite(recombIntensity)==false) {
            throw new IllegalArgumentException(String.valueOf(recombIntensity));
        }
        return IntStream.range(0, genDist.size())
                .parallel()
                .mapToDouble(m -> -Math.expm1(-0.01*genDist.get(m)*recombIntensity))
                .toArray();
    }

    /**
     * Returns an array of length {@code (2*ia.length)} obtained by
     * repeating each element in the specified array once.  The returned array
     * is equal to
     * <pre>
            IntStream.range(0, ia.length&lt;&lt;1)
                    .parallel()
                    .map(j -&gt; ia[j&gt;&gt;1])
                    .toArray();
     * </pre>
     * @param ia an integer array
     * @return an array obtained by repeating each element in the
     * specified array once
     * @throws NullPointerException if {@code ia == null}
     */
    public static int[] repeatEachElement(int[] ia) {
        return IntStream.range(0, ia.length<<1)
                .parallel()
                .map(j -> ia[j>>1])
                .toArray();
    }

    /**
     * Returns an array mapping each sample index to the sample's string map
     * value. The map file is a white-space delimited text file whose first
     * column contains sample identifiers and whose second column contains
     * the value associated with each sample identifier.  Multiple consecutive
     * white-space characters are treated as a single delimiter, and any
     * white-space at the beginning or end of a line is ignored.  The java
     * virtual machine will exit with an error message if any sample identifier
     * in the specified {@code sampleIds} array does not have a value in the
     * specified {@code mapFile}.
     * @param samples the list of sample identifiers
     * @param mapFile a white-space delimited text file that maps each sample
     * identifier to the sample's map value
     * @return an array mapping each sample index to the sample's string map
     * value
     * @throws IllegalArgumentException if any line in the specified
     * {@code mapFile} that contains a non-white-space character does not
     * contain exactly two white-space delimited fields
     * @throws NullPointerException if
     * {@code sampleIds == null || mapFile == null}
     */
    public static String[] sampleMap(Samples samples, File mapFile) {
        HashMap<String, String> map = readMap(mapFile);
        return IntStream.range(0, samples.size())
                .parallel()
                .mapToObj(j -> {
                        String panel = map.get(samples.id(j));
                        if (panel==null) {
                            String err = "No value associated with sample identifier";
                            String info = Const.nl + "Error     :  " + err
                                    + Const.nl     + "Filename  :  " + mapFile
                                    + Const.nl     + "Sample    :  " + samples.id(j);
                            Utilities.exit(new Throwable(err), info);
                        }
                        return panel;
                    }
                )
                .toArray(String[]::new);
    }

    /**
     * Returns the a hash map whose keys are the elements of the first column
     * and whose values are the elements of the second column in the specified
     * white-space delimited text file.  Multiple consecutive white-space
     * characters are treated as a single delimiter, and any white-space
     * at the beginning or end of a line is ignored.
     * The Java virtual machine will exit with an error message if any
     * non-blank line does not have exactly two white-space delimited fields
     * or if any entry in the first column is duplicated.
     * @param mapFile a text file
     * @return the hash map
     * @throws NullPointerException if {@code mapFile == null}
     */
    public static HashMap<String, String> readMap(File mapFile) {
        String[] lines = readLines(mapFile);
        HashMap<String, String> map = new HashMap<>(lines.length);
        for (int j=0; j<lines.length; ++j) {
            String line = lines[j].trim();
            if (line.length()>0) {
                String[] fields = StringUtil.getFields(line, 3);
                if (fields.length != 2) {
                    String err = "Line does not have exactly two white-space delimited fields";
                    String info = Const.nl + "Error     :  " + err
                            + Const.nl     + "Filename  :  " + mapFile
                            + Const.nl     + "Line      :  " + lines[j]
                            + Const.nl     + "Fields    :  " + Arrays.toString(fields);
                    Utilities.exit(new Throwable(err), info);
                }
                String prevValue = map.put(fields[0], fields[1]);
                if (prevValue!=null) {
                    String err = "Duplicate identifier in the first column";
                    String info = Const.nl + "Error      :  " + err
                            + Const.nl     + "Filename   :  " + mapFile
                            + Const.nl     + "Identifier :  " + fields[0]
                            + Const.nl     + "Value_1    :  " + prevValue
                            + Const.nl     + "Value_2    :  " + fields[1];
                    Utilities.exit(new Throwable(err), info);
                }
            }
        }
        return map;
    }

    /**
     * Returns a {@code HashMap} that maps each string array element to its
     * index.
     * @param sa a string array
     * @return a {@code HashMap} that maps each string array element to its
     * index
     * @throws NullPointerException if {@code sa == null} or if any array
     * element is {@code null}
     * @throws IllegalArgumentException if any two array elements are equal
     */
    public static HashMap<String, Integer> inverseMap(String[] sa) {
        HashMap<String, Integer> inverseMap = new HashMap<>();
        for (int j=0; j<sa.length; ++j) {
            if (sa[j]==null) {
                throw new NullPointerException("sa[" + j + "]==null");
            }
            Integer prevValue = inverseMap.put(sa[j], j);
            if (prevValue!=null) {
                throw new IllegalArgumentException("Duplicate entry: " + sa[j]);
            }
        }
        return inverseMap;
    }

    /**
     * Returns the lines of the specified file.  The java virtual machine will
     * exit with a stack trance and an error message if the specified file
     * does not exist or is a directory.
     * @param file a text file
     * @return the non-blank lines of the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static String[] readLines(File file) {
        String err = null;
        if (file.exists()==false) {
            err = "File does not exist";
        }
        else if (file.isDirectory()) {
            err = "File cannot be a directory";
        }
        if (err!=null) {
            String info = Const.nl + "Error     :  " + err
                    + Const.nl     + "Filename  :  " + file;
            Utilities.exit(new Throwable(err), info);
        }
        ArrayList<String> lines = new ArrayList<>();
        try (FileIt<String> it = InputIt.fromGzipFile(file)) {
            while (it.hasNext()) {
                String trimmed = it.next().trim();
                if (trimmed.length()>0) {
                    lines.add(trimmed);
                }
            }
        }
        return lines.toArray(new String[0]);
    }

    /**
     * Returns the list of filtered samples in the specified VCF file.
     * The Java Virtual Machine will exit with a stack trace and an error
     * message if the VCF meta-information lines are not immediately followed
     * by a properly formatted VCF header line.
     * @param vcfFile a VCF file
     * @param sampFilter a sample filter
     * @return the list of filtered samples in the specified VCF file
     * @throws NullPointerException if
     * {@code vcfFile == null || sampFilter == null}}
     */
    public static String[] readSampleIds(File vcfFile, Filter<String> sampFilter) {
        String headerPrefix = VcfHeader.HEADER_PREFIX + Const.tab;
        try (FileIt<String> it = InputIt.fromGzipFile(vcfFile)) {
            String line = "##";
            while (it.hasNext() && line.startsWith("##")) {
                line = it.next();
            }
            if (line.startsWith(headerPrefix)) {
                int firstSampleField = 9;
                String[] fields = StringUtil.getFields(line, Const.tab);
                return Arrays.stream(fields)
                        .parallel()
                        .skip(firstSampleField)
                        .filter(sampFilter::accept)
                        .toArray(String[]::new);
            }
            else {
                String msg = "Incorrectly formatted beginning of VCF header line";
                Exception ex = new Exception(msg);
                String err = Const.nl + "Error  :  " + msg
                        + Const.nl    + "File   :  " + vcfFile
                        + Const.nl    + "Line   :  "
                        + line.substring(0, headerPrefix.length()+10) + "[...]";
                Utilities.exit(ex, err);
            }
        }
        return null;
    }

    /**
     * Appends the specified vector to the specified {@code StringBuilder}.
     * The elements of the matrix will be tab-delimited, and the vector
     * will be terminated with a new line character.
     * @param sb the string builder to which the matrix will be appended
     * @param da an array of {@code double} values
     * @throws NullPointerException if {@code (sb == null) || (da == null)}
     */
    public static void appendLineWithVector(StringBuilder sb, double[] da) {
        for (int j=0; j<da.length; ++j) {
            if (j>0) {
                sb.append(Const.tab);
            }
            sb.append((float) da[j]);
        }
        sb.append(Const.nl);
    }

    /**
     * Appends the specified matrix to the specified {@code StringBuilder}.
     * The elements of the matrix will be tab-delimited, and each row of
     * the matrix will be terminated with a new line character.
     * @param sb the string builder to which the matrix will be appended
     * @param da a two-dimensional array of {@code double} values
     * @throws NullPointerException if {@code (sb == null) || (da == null)}
     * or if any row of the matrix is {@code null}
     */
    public static void appendLinesWithMatrix(StringBuilder sb, double[][] da) {
        for (int i=0; i<da.length; ++i) {
            appendLineWithVector(sb, da[i]);
        }
    }

    /**
     * Returns a copy of the specified two-dimensional array
     * @param array a two-dimensional array
     * @return a copy of the specified two-dimensional array
     * @throws NullPointerException if {@code arrayt == null}
     */
    public static double[][] cloneArray(double[][] array) {
        return IntStream.range(0, array.length)
                .mapToObj(i -> array[i].clone())
                .toArray(double[][]::new);
    }
}

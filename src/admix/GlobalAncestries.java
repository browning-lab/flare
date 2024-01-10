/*
 * Copyright 2021-2023 Brian L. Browning
 * Copyright 2024 Genomics plc
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
import blbutil.InputIt;
import blbutil.StringUtil;
import blbutil.Utilities;
import ints.IntList;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.IntStream;

/**
 * <p>Class {@code GlobalAncestries} stores global ancestry proportions for
 * some or all of the target samples.</p>
 *
 * <p>Instances of class {@code GlobalAncestries} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class GlobalAncestries {

    private final int nAnc;
    private final int nTargSamples;
    private final int[] sampleIndices;
    private final double[][] ancestryProbs;

    /**
     * Constructs a new {@code GlobalAncestries} instance from the specified
     * data.
     * @param nAnc the number of ancestries
     * @param nTargSamples the number of target samples
     * @throws IllegalArgumentException if {@code nAnc < 0}
     * @throws IllegalArgumentException if {@code nTargSamples < 0}
     */
    public GlobalAncestries(int nAnc, int nTargSamples) {
        if (nAnc < 1) {
            throw new IllegalArgumentException(String.valueOf(nAnc));
        }
        if (nTargSamples < 1) {
            throw new IllegalArgumentException(String.valueOf(nTargSamples));
        }
        this.nAnc = nAnc;
        this.nTargSamples = nTargSamples;
        this.sampleIndices = new int[0];
        this.ancestryProbs = new double[0][];
    }

    /**
     * Constructs a new {@code GlobalAncestries} instance from the specified
     * data. The Java Virtual Machine will exit with an error message if
     * an I/O error is encountered or if the list of ancestry identifiers
     * in {@code gtAncFile} is not equal to {@code expectedAncIds}.
     * @param gtAncFile a file with sample ancestry proportions
     * @param expectedAncIds the expected list of ancestry identifiers
     * @param targSampleIds the target sample identifiers
     * @throws IllegalArgumentException if {@code nAnc < 0}
     * @throws NullPointerException if {@code (gtAncFile == null)}
     * @throws NullPointerException if {@code ((expAncIds == null)} or if
     * any element of {@code expAncIds} is null
     * @throws NullPointerException if {@code ((targSampleIds == null)} or if
     * any element of {@code targSampleIds} is null
     */
    public GlobalAncestries(File gtAncFile, String[] expectedAncIds,
            String[] targSampleIds) {
        if (expectedAncIds.length==0) {
            throw new IllegalArgumentException(String.valueOf(expectedAncIds.length));
        }
        HashMap<String, Integer> sampleToIndexMap = AdmixUtils.inverseMap(targSampleIds);
        int maxLines = 5000;
        ArrayList<String> lines = new ArrayList<>(maxLines);
        IntList tempSampleIndices = new IntList();
        ArrayList<double[]> tempAncestryProbs = new ArrayList<>();
        this.nAnc = expectedAncIds.length;
        String[] foundAncIds = readAncIds(gtAncFile);
        if (!Arrays.equals(foundAncIds, expectedAncIds)) {
            String err = "Error: expected ancestries in gt-ancestries file to be: "
                    + Arrays.toString(expectedAncIds) + ", found: " + Arrays.toString(expectedAncIds);
            printErrAndExit(err, gtAncFile, null);
        }

        try (FileIt<String> it = InputIt.fromGzipFile(gtAncFile)) {
            readLines(it, lines, 1); // read header line
            if (lines.size()==1) {
                checkAncestryIds(gtAncFile, lines.get(0), expectedAncIds);
            }
            while (it.hasNext()) {
                readLines(it, lines, maxLines);
                GtAncLine[] gtAncLines = parseLines(gtAncFile, nAnc, lines,
                        sampleToIndexMap);
                for (GtAncLine gtAncLine : gtAncLines) {
                    tempSampleIndices.add(gtAncLine.sampleIndex);
                    tempAncestryProbs.add(gtAncLine.ancProbs);
                }
                lines.clear();
            }
        }
        checkForDuplicateSamples(gtAncFile, targSampleIds, tempSampleIndices);
        this.nTargSamples = targSampleIds.length;
        this.sampleIndices = tempSampleIndices.toArray();
        this.ancestryProbs = tempAncestryProbs.toArray(new double[0][]);
    }

    private static void readLines(FileIt<String> it, ArrayList<String> lines,
            int maxLines) {
        while (it.hasNext() && lines.size()<maxLines) {
            String candidate = it.next().trim();
            if (candidate.length()>0) {
                lines.add(candidate);
            }
        }
    }

    private static GtAncLine[] parseLines(File gtAncFile, int nAnc,
            ArrayList<String> lines,
            HashMap<String, Integer> sampleToIndexMap) {
        return lines.stream()
                .parallel()
                .map(line -> StringUtil.getFields(line))
                .filter(fields -> sampleToIndexMap.get(fields[0])!=null)
                .map(fields -> new GtAncLine(gtAncFile, nAnc, fields, sampleToIndexMap))
                .toArray(GtAncLine[]::new);
    }

    private static class GtAncLine {
        private final int sampleIndex;
        private final double[] ancProbs;

        private GtAncLine(File file, int nAnc, String[] fields,
                HashMap<String, Integer> sampleToIndexMap) {
            if (fields.length!=(nAnc+1)) {
                String err = "Error: line does not have " + (nAnc + 1) + " fields";
                printErrAndExit(err, file, Arrays.toString(fields));
            }
            double[] da = new double[nAnc];
            for (int j=0; j<da.length; ++j) {
                da[j] = parseDouble(file, fields, fields[j+1]);
            }
            this.sampleIndex = sampleToIndexMap.get(fields[0]);
            this.ancProbs = da;
        }
    }

    private static double parseDouble(File file, String[] fields, String value) {
        double d = Double.NaN;
        try {
            d = Double.parseDouble(value);
        }
        catch (NumberFormatException e) {
            String err = "String cannot be parsed as a number: " + value;
            printErrAndExit(err, file, Arrays.toString(fields));
        }
        return d;
    }

    private static void checkForDuplicateSamples(File gtAncFile,
            String[] sampleIds, IntList sampleIndexList) {
        int[] indices = sampleIndexList.toArray();
        Arrays.sort(indices);
        for (int j=1; j<indices.length; ++j) {
            if (indices[j]==indices[j-1]) {
                String err = "Duplicate sample identifier in first column of "
                        + "gt-ancestries file: " + sampleIds[indices[j]];
                printErrAndExit(err, gtAncFile, null);
            }
        }
    }

    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return nAnc;
    }

    /**
     * Returns an array with length
     * {@code this.fixedParams().targSamples().size()} whose {@code j}-th
     * entry is an array with global ancestry proportions for the {@code j}-th
     * target sample. The contract for this method is undefined if the
     * the elements of the specified {@code defaultAncProbs} array are not
     * finite, non-negative numbers.
     * @param defaultAncProbs the global ancestry proportions for samples
     * that do not have global ancestry proportions stored in {@code this}
     * @return the ancestry proportions for target samples
     * @throws NullPointerException if {@code defaultAncProbs == null}
     * @throws IllegalArgumentException if
     * {@code defaultAncProbs.length != this.fixedParams().nAnc()}
     */
    public double[][] globalAncestries(double[] defaultAncProbs) {
        if (defaultAncProbs.length != this.nAnc) {
            throw new IllegalArgumentException(String.valueOf(defaultAncProbs.length));
        }
        double[] da = normalize(defaultAncProbs);
        double[][] ancProbs = IntStream.range(0, nTargSamples)
                .parallel()
                .mapToObj(j -> da)
                .toArray(double[][]::new);
        for (int j=0; j<sampleIndices.length; ++j) {
            ancProbs[sampleIndices[j]] = this.ancestryProbs[j];
        }
        return ancProbs;
    }

    private static double[] normalize(double[] defaultAncProbs) {
        double sum = AdmixUtils.sum(defaultAncProbs);
        double[] copy = defaultAncProbs.clone();
        AdmixUtils.scale(copy, sum);
        return copy;
    }

    /**
     * Returns the list of ancestries in a sample ancestries file header line.
     * The Java virtual machine will exit with an error message if an
     * I/O error is encountered, or if the first line of the sample ancestries
     * file is incorrectly formatted.
     * @param gtAncFile a sample ancestries file
     * @return the list of ancestries in a sample ancestries file header line
     * @throws NullPointerException if {@code gtAncFile == null}
     */
    public static String[] readAncIds(File gtAncFile) {
        String headerLine = "";
        try (FileIt<String> it = InputIt.fromGzipFile(gtAncFile)) {
            while (it.hasNext() && headerLine.length()==0) {
                String candidate = it.next().trim();
                if (candidate.length()>0) {
                    headerLine = candidate;
                }
            }
        }
        return extractAncIds(gtAncFile, headerLine);
    }

    private static void checkAncestryIds(File gtAncFile, String headerLine,
            String[] ancIds) {
        String[] fields = StringUtil.getFields(headerLine);
        for (int j=0; j<ancIds.length; ++j) {
            if (ancIds[j].equals(fields[j+1])==false) {
                String err = "Error: expected " + j + "-th ancestry to be"
                        + ancIds[j] + ", found " + fields[j+1];
                printErrAndExit(err, gtAncFile, headerLine);
            }
        }
    }

    private static String[] extractAncIds(File gtAncFile, String headerLine) {
        String[] fields = StringUtil.getFields(headerLine);
        if (fields.length < 2) {
            String s = "Missing ancestries in gt-ancestry file header line";
            printErrAndExit(s, gtAncFile, headerLine);
        }
        if (fields[0].equalsIgnoreCase("SAMPLE")==false) {
            String s = "First field in gt-ancestry file header line is not \"SAMPLE\"";
            printErrAndExit(s, gtAncFile, headerLine);
        }
        return Arrays.copyOfRange(fields, 1, fields.length);
    }

    private static void printErrAndExit(String err, File gtAncFile, String line) {
        StringBuilder sb = new StringBuilder(1<<9);
        sb.append("Error:  ");
        sb.append(err);
        sb.append(Const.nl);
        sb.append("File:   ");
        sb.append(gtAncFile);
        sb.append(Const.nl);
        if (line != null) {
            sb.append("Line:   ");
            sb.append(line);
        }
        Utilities.exit(sb.toString());
    }
}

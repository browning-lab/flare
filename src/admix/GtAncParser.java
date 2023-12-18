/*
 * Copyright 2023 Genomics plc
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

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import blbutil.Const;
import blbutil.StringUtil;
import blbutil.Utilities;

/**
 * <p>Class {@code GtAncParser} parses the file passed in using the `gt-ancestries`
 * option.</p>
 *
 * @author Thomas Hickman {@code <thomas.hickman@genomicsplc.com>}
 */
public class GtAncParser {
    private double[][] sampleToMu;
    private Set<String> samplesParsed;

    private final Map<String, Integer> sampleIdsToSampleIndex;
    private final int nAncestries;
    private final File gtAncestriesFile;

    /**
     * Constructs a new {@code GtAncParser} instance for the specified file.
     * analysis
     * @param gtAncestriesFile the file to parse
     * @param sampleIdsInGtFile the sample IDs in the target genotype file 
     * (used for cross-checking the gt-ancestries file)
     * @param ancIds the ancestry IDs used in ancestry estimation (used 
     * for cross-checking the gt-ancestries file)
     */
    public GtAncParser(File gtAncestriesFile, String[] sampleIdsInGtFile, String[] ancIds) {
        this.sampleToMu = new double[sampleIdsInGtFile.length][ancIds.length];

        this.sampleIdsToSampleIndex = new HashMap<String, Integer>(sampleIdsInGtFile.length);
        for (int i = 0; i < sampleIdsInGtFile.length; i++) {
            sampleIdsToSampleIndex.put(sampleIdsInGtFile[i], i);
        }

        this.gtAncestriesFile = gtAncestriesFile;
        this.nAncestries = ancIds.length;

        this.samplesParsed = new HashSet<String>();
        String[] lines = AdmixUtils.readLines(gtAncestriesFile);

        checkNonEmpty(lines);
        checkHeader(lines[0], ancIds);

        for (String line : Arrays.copyOfRange(lines, 1, lines.length)) {
            parseLine(line);
        }

        checkHaveParsedRequiredSamples(sampleIdsInGtFile);
    }
    
    private void parseLine(String line) {
        String[] fields = StringUtil.getFields(line);
        // Note: this is never going to fail, as `readLines` (called in the caller) does not pass us blank lines.
        String sampleId = fields[0];
        checkNonDuplicateSample(sampleId);

        double[] mu = Arrays.stream(fields).skip(1).mapToDouble(Double::parseDouble).toArray();

        checkMuIsAProportion(sampleId, mu);
        checkMuLength(sampleId, mu);

        Integer sampleIndex = sampleIdsToSampleIndex.get(sampleId);
        if (sampleIndex == null) {
            // Ignore samples in the gt-ancestries file that aren't in the target file.
            return;
        }

        sampleToMu[sampleIdsToSampleIndex.get(sampleId)] = mu;
        samplesParsed.add(sampleId);
    }

    private void checkMuLength(String sampleId, double[] mu) {
        if (mu.length != nAncestries) {
            String err = "The gt-ancestries file contains a line with an incorrect number of ancestry proportions";
            String info = Const.nl + "Error                                   :  " + err
                        + Const.nl + "Filename                                :  " + gtAncestriesFile
                        + Const.nl + "Sample ID                               :  " + sampleId
                        + Const.nl + "Expected number of ancestry proportions :  " + nAncestries
                        + Const.nl + "Found number of ancestry proportions    :  " + mu.length;

            Utilities.exit(new Throwable(err), info);
        }
    }

    private void checkMuIsAProportion(String sampleId, double[] mu) {
        for (double prob : mu) {
            if (prob < 0f || prob > 1f) {
                String err = "The gt-ancestries file contains a probability that is less than 0.0 or greater than 1.0";
                String info = Const.nl + "Error       :  " + err
                            + Const.nl + "Filename    :  " + gtAncestriesFile
                            + Const.nl + "Sample ID   :  " + sampleId
                            + Const.nl + "Probability :  " + prob;

                Utilities.exit(new Throwable(err), info);
            }
        }

        double tolerance = 0.01f;
        double muSum = AdmixUtils.sum(mu);
        if (Math.abs(1f - muSum) > tolerance) {
            String err = "The gt-ancestries file contains a sample which has probabilities that do not sum to 1";
            String info = Const.nl + "Error           :  " + err
                        + Const.nl + "Filename        :  " + gtAncestriesFile
                        + Const.nl + "Sample ID       :  " + sampleId
                        + Const.nl + "Probability sum :  " + muSum;

            Utilities.exit(new Throwable(err), info);
        }
    }

    private void checkNonDuplicateSample(String sampleId) {
        if (samplesParsed.contains(sampleId)) {
            String err = "The gt-ancestries file contains a duplicate sample ID";
            String info = Const.nl + "Error     :  " + err
                        + Const.nl + "Filename  :  " + gtAncestriesFile
                        + Const.nl + "Sample ID :  " + sampleId;

            Utilities.exit(new Throwable(err), info);
        }
    }

    private void checkHaveParsedRequiredSamples(String[] sampleIdsInGenotypeFile) {
        if (!samplesParsed.containsAll(Arrays.asList(sampleIdsInGenotypeFile))) {
            String err = "The gt-ancestries file does not contain all the samples in the genotype file";
            String info = Const.nl + "Error          :  " + err
                        + Const.nl + "Filename       :  " + gtAncestriesFile;

            Utilities.exit(new Throwable(err), info);
        }
    }

    private void checkNonEmpty(String[] lines) {
        if (lines.length == 0) {
            String err = "The gt-ancestries file does not contain any content";
            String info = Const.nl + "Error          :  " + err
                        + Const.nl + "Filename       :  " + gtAncestriesFile;

            Utilities.exit(new Throwable(err), info);
        }
    }

    private void checkHeader(String headerLine, String[] ancIds) {
        String[] headerFields = StringUtil.getFields(headerLine);
        String[] headerAncestries = Arrays.copyOfRange(headerFields, 1, headerFields.length);

        if (!Arrays.equals(headerAncestries, ancIds)) {
            String err = "The gt-ancestries file does not contain a header with the same ancestries as are defined in the model file";
            String info = Const.nl + "Error             :  " + err
                        + Const.nl + "Filename          :  " + gtAncestriesFile
                        + Const.nl + "File Ancestries   :  " + Arrays.toString(headerAncestries)
                        + Const.nl + "Model Ancestries  :  " + Arrays.toString(ancIds);

            Utilities.exit(new Throwable(err), info);
        }
    }

    /**
     * Gets the ancestry proportions for a given sample.
     * @param sampleI the index of the sample to query
     * @returns an array of ancestry proportions (which can be used in the 
     * `mu` parameter).
     */
    public double[] getMu(int sampleI) {
        return sampleToMu[sampleI];
    }
}

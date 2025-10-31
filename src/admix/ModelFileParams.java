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
import blbutil.InputIt;
import blbutil.StringUtil;
import blbutil.Utilities;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import vcf.Samples;

/**
 * <p>Class {@code ModelFileParams} represents the analysis parameters for a
 * local ancestry inference analysis.</p>
 *
 * <p>Instances of class {@code ModelFileParams} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ModelFileParams implements ParamsInterface {

    private final SampleData sampleData;
    private final double t;                    // gen since admixture
    private final double[] mu;                 // ancestry proportions
    private final double[][] theta;            // miscopy probability matrix
    private final double[][] p;                // copying probability matrix
    private final double[] rho;                // population switch probabilities

    /**
     * Constructs a new {@code ModelFileParams} instance for the specified
     * data. The Java virtual machine will exit with an error if an error in
     * the format of the {@code sampleData.par().model()} file is detected.
     * @param sampleData reference and target sample metadata
     * @throws NullPointerException if
     * {@code (sampleData == null) || (sampleData.par().modelFile() == null)}
     */
    public ModelFileParams(SampleData sampleData) {
        AdmixPar par = sampleData.par();
        Samples refSamples = sampleData.refSamples();
        String[] sampleToRefPanel = AdmixUtils.sampleMap(refSamples, par.ref_panel());
        ArrayList<String> lines = readParamsData(par.model());

        String[] ancIds = StringUtil.getFields(lines.get(0));
        String[] refPanelIds = StringUtil.getFields(lines.get(1));
        SampleData.confirmMoreThanOneAncestry(ancIds);
        AdmixRefPanels.checkRefPanelIds(refSamples, sampleToRefPanel, refPanelIds);

        int nAnc = ancIds.length;
        double[][] dblValues = extractDoubles(par.model(), lines, nAnc, refPanelIds);

        this.sampleData = sampleData;
        this.t = dblValues[0][0];
        this.mu = dblValues[1];
        this.p = Arrays.copyOfRange(dblValues, 2, nAnc+2);
        this.theta = Arrays.copyOfRange(dblValues, nAnc+2, 2*nAnc+2);
        this.rho = dblValues[2*nAnc+2];
    }

    private static ArrayList<String> readParamsData(File paramsFile) {
        ArrayList<String> lines = new ArrayList<>();
        int maxLines = Integer.MAX_VALUE;
        int nAnc = -1;
        try (FileIt<String> it = InputIt.fromGzipFile(paramsFile)) {
            while (it.hasNext() && lines.size()<maxLines) {
                String candidate = it.next().trim();
                if (candidate.length()>0 && candidate.charAt(0)!='#') {
                    lines.add(candidate);
                    if (lines.size()==1) {
                        nAnc = StringUtil.countFields(candidate);
                        maxLines = 10 + 5 + 2*nAnc;
                    }
                }
            }
        }
        checkDataLineCount(paramsFile, nAnc, lines);
        return lines;
    }

    private static void checkDataLineCount(File paramsFile, int nAnc,
            ArrayList<String> lines) {
        int nExpectedLines = 5 + 2*nAnc;
        if (lines.size()!=nExpectedLines) {
            String nFound = (lines.size()<nExpectedLines)
                    ? String.valueOf(lines.size())
                    : "more than " + nExpectedLines;
            String msg = "Parameters file does not have " + nExpectedLines
                    + " data lines";
            Exception ex = new Exception(msg);
            StringBuilder sb = new StringBuilder(1<<10);
            sb.append("Error:       ");
            sb.append(msg);
            sb.append(Const.nl);
            sb.append("File:        ");
            sb.append(paramsFile);
            sb.append(Const.nl);
            sb.append("Ancestries:  ");
            sb.append(nAnc);
            sb.append(Const.nl);
            sb.append("Data lines:  ");
            sb.append(nFound);
            sb.append(Const.nl);
            sb.append(Const.nl);
            sb.append("First ");
            sb.append(lines.size());
            sb.append(" data lines:");
            sb.append(Const.nl);
            for (int j=0, n=lines.size(); j<n; ++j) {
                sb.append(Const.nl);
                sb.append(lines.get(j));
            }
            Utilities.exit(ex, sb.toString());
        }
    }

    private static double[][] extractDoubles(File paramsFile,
            List<String> lines, int nAnc, String[] refPanelIds) {
        List<String> doubleLines = lines.subList(2, lines.size());
        double[][] params = doubleLines.stream()
                .parallel()
                .map(line -> parseDoubleArray(paramsFile, line))
                .toArray(double[][]::new);

        checkParams(paramsFile, nAnc, refPanelIds.length, doubleLines, params);
        return params;
    }

    private static double[] parseDoubleArray(File file, String line) {
        String[] fields = StringUtil.getFields(line);
        double[] da = new double[fields.length];
        for (int j=0; j<fields.length; ++j) {
            da[j] = parseDouble(file, line, fields[j]);
        }
        return da;
    }

    private static double parseDouble(File file, String line, String value) {
        double d = Double.NaN;
        try {
            d = Double.parseDouble(value);
        }
        catch (NumberFormatException e) {
            String msg = "Non-numerical value: " + value;
            Exception ex = new Exception(msg);
            StringBuilder sb = new StringBuilder(1<<9);
            sb.append("Error:  ");
            sb.append(msg);
            sb.append(Const.nl);
            sb.append("File:   ");
            sb.append(file);
            sb.append(Const.nl);
            sb.append("Line:   ");
            sb.append(line);
            Utilities.exit(ex, sb.toString());
        }
        return d;
    }

    private static void checkParams(File paramsFile, int nAnc, int nRefPanels,
            List<String> lines, double[][] params) {
        // T
        checkLength(paramsFile, nAnc, nRefPanels, lines.get(0), params[0], 1);
        // mu
        checkLength(paramsFile, nAnc, nRefPanels, lines.get(1), params[1], nAnc);
        checkProbs(paramsFile, lines.get(1), params[1]);
        checkProbSum(paramsFile, lines.get(1), params[1]);
        // p
        int start = 2;
        int end = start + nAnc;
        for (int j=start; j<end; ++j) {
            checkLength(paramsFile, nAnc, nRefPanels, lines.get(j), params[j], nRefPanels);
            checkProbs(paramsFile, lines.get(j), params[j]);
            checkProbSum(paramsFile, lines.get(j), params[j]);
        }
        // theta
        start = end;
        end = start + nAnc;
        for (int j=start; j<end; ++j) {
            checkLength(paramsFile, nAnc, nRefPanels, lines.get(j), params[j], nRefPanels);
            checkProbs(paramsFile, lines.get(j), params[j]);
        }
        // rho
        checkLength(paramsFile, nAnc, nRefPanels, lines.get(end), params[end], nAnc);
        checkNonNegative(paramsFile, lines.get(end), params[end]);
    }

    private static void checkLength(File paramsFile, int nAnc, int nRefPanels,
            String line, double[] params, int expLength) {
        if (params.length != expLength) {
            String msg = "Line does not have have " + expLength + " fields";
            Exception ex = new Exception(msg);
            StringBuilder sb = new StringBuilder(1<<9);
            sb.append("Error:       ");
            sb.append(msg);
            sb.append(Const.nl);
            sb.append("File:        ");
            sb.append(paramsFile);
            sb.append(Const.nl);
            sb.append("Ancestries:  ");
            sb.append(nAnc);
            sb.append(Const.nl);
            sb.append("Ref panels:  ");
            sb.append(nRefPanels);
            sb.append(Const.nl);
            sb.append("Line:        ");
            sb.append(line);
            Utilities.exit(ex, sb.toString());
        }
    }

    private static void checkProbs(File paramsFile, String line, double[] probs) {
        for (double prob : probs) {
            if (prob<0f || prob>1f) {
                String msg = "Probability is less than 0.0 or greater than 1.0 ["
                        + prob + "]";
                Exception ex = new Exception(msg);
                StringBuilder sb = new StringBuilder(1<<9);
                sb.append("Error:  ");
                sb.append(msg);
                sb.append(Const.nl);
                sb.append("File:   ");
                sb.append(paramsFile);
                sb.append(Const.nl);
                sb.append("Line:   ");
                sb.append(line);
                Utilities.exit(ex, sb.toString());
            }
        }
    }

    private static void checkProbSum(File paramsFile, String line, double[] params) {
        double tolerance = 0.01f;
        double sum = AdmixUtils.sum(params);
        if (Math.abs(1f - sum)>tolerance) {
            String msg = "Probabilities do not sum to 1 [sum: " + sum + "]";
            Exception ex = new Exception(msg);
            StringBuilder sb = new StringBuilder(1<<9);
            sb.append("Error:  ");
            sb.append(msg);
            sb.append(Const.nl);
            sb.append("File:   ");
            sb.append(paramsFile);
            sb.append(Const.nl);
            sb.append("Line:   ");
            sb.append(line);
            Utilities.exit(ex, sb.toString());
        }
    }

    private static void checkNonNegative(File paramsFile, String line,
            double[] params) {
        for (double param : params) {
            if (param<0) {
                String msg = "Parameter value is negative [" + param + "]";
                Exception ex = new Exception(msg);
                StringBuilder sb = new StringBuilder(1<<9);
                sb.append("Error:  ");
                sb.append(msg);
                sb.append(Const.nl);
                sb.append("File:   ");
                sb.append(paramsFile);
                sb.append(Const.nl);
                sb.append("Line:   ");
                sb.append(line);
                Utilities.exit(ex, sb.toString());
            }
        }
    }

    @Override
    public SampleData sampleData() {
        return sampleData;
    }

    @Override
    public double T() {
        return t;
    }

    @Override
    public double[] studyMu() {
        return mu.clone();
    }

    @Override
    public double rho(int i) {
        return rho[i];
    }

    @Override
    public double[] rho() {
        return rho.clone();
    }

    @Override
    public double p(int i, int j) {
        return p[i][j];
    }

    @Override
    public double[][] p() {
        return AdmixUtils.cloneArray(p);
    }

    @Override
    public double[][] theta() {
        return AdmixUtils.cloneArray(theta);
    }

    @Override
    public String toString() {
        return ParamUtils.toString(this);
    }
}

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
import ints.IntArray;
import java.io.File;
import java.util.Arrays;
import vcf.Samples;

/**
 * <p>Class {@code SampleData} stores reference and target
 * sample metadata.</p>
 *
 * <p>Instances of class {@code SampleData} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SampleData {

    private final AdmixPar par;
    private final Samples targSamples;
    private final Samples refSamples;
    private final String[] refPanelIds;
    private final IntArray refHapToPanel;
    private final IntArray nPanelHaps;
    private final String[] ancIds;
    private final GlobalAncestries globalAncestries;

    /**
     * Returns a {@code SampleData} instance that is constructed from the
     * specified data. The Java virtual machine will exit with an error message
     * if an inconsistency in the input data is detected.
     * @param par the command line parameters
     * @param refSamples the list of reference samples identifiers
     * @param targSamples the list of target samples identifiers
     * @return a {@code SampleData} instance
     * @throws NullPointerException if
     * {@code (par == null) || (refSamples == null) || (targSamples == null}}
     */
    public static SampleData create(AdmixPar par, Samples refSamples,
            Samples targSamples) {
        if (par.model()==null) {
            return new SampleData(par, refSamples, targSamples);
        }
        else {
            String[] lines = readAncAndPanelLines(par.model());
            String[] ancIds = StringUtil.getFields(lines[0]);
            String[] refPanelIds = StringUtil.getFields(lines[1]);
            return new SampleData(par, refSamples, targSamples, ancIds,
                    refPanelIds);
        }
    }

    /**
     * Constructs and returns an {@code SampleData} instance for the
     * specified data. The Java virtual machine will exit with an error
     * message if an inconsistency in the input data is detected.
     * @param par the command line parameters
     * @param refSamples the list of reference samples
     * @param targSamples the list of target samples
     * @param refSampleIds the list of reference sample identifiers
     * @param ancIds the list of ancestry identifiers
     * @param refPanelIds the list of reference panel identifiers
     * @throws NullPointerException if any parameter is {@code null}
     */
    private SampleData(AdmixPar par, Samples refSamples, Samples targSamples,
            String[] ancIds, String[] refPanelIds) {
        this.par = par;
        this.refSamples = refSamples;
        this.targSamples = targSamples;
        String[] sampleToRefPanel = AdmixUtils.sampleMap(refSamples, par.ref_panel());
        this.refPanelIds = refPanelIds.clone();
        confirmMoreThanOneAncestry(ancIds);
        AdmixRefPanels.checkRefPanelIds(refSamples, sampleToRefPanel, refPanelIds);
        this.refHapToPanel = AdmixRefPanels.hapToPanel(sampleToRefPanel, refPanelIds);
        this.nPanelHaps = AdmixRefPanels.nPanelHaps(refHapToPanel, refPanelIds.length);
        this.ancIds = ancIds;
        this.globalAncestries = par.gt_ancestries()==null
                ? new GlobalAncestries(ancIds.length, targSamples.size())
                : new GlobalAncestries(par.gt_ancestries(), ancIds, targSamples.ids());
    }

    private SampleData(AdmixPar par, Samples refSamples, Samples targSamples,
            String[] ancIds) {
        if (targSamples == null) {
            throw new NullPointerException(Samples.class.toString());
        }
        this.par = par;
        this.refSamples = refSamples;
        this.targSamples = targSamples;
        String[] refSampToPanelId = AdmixUtils.sampleMap(refSamples, par.ref_panel());
        this.refPanelIds = AdmixRefPanels.indexPanels(refSampToPanelId);
        this.refHapToPanel = AdmixRefPanels.hapToPanel(refSampToPanelId,
                refPanelIds);
        this.ancIds = ancIds.clone();
        this.nPanelHaps = AdmixRefPanels.nPanelHaps(refHapToPanel,
                refPanelIds.length);
        this.globalAncestries = par.gt_ancestries()==null
                ? new GlobalAncestries(ancIds.length, targSamples.size())
                : new GlobalAncestries(par.gt_ancestries(), ancIds, targSamples.ids());
    }

    private SampleData(AdmixPar par, Samples refSamples, Samples targSamples) {
        if (targSamples == null) {
            throw new NullPointerException(Samples.class.toString());
        }
        this.par = par;
        this.refSamples = refSamples;
        this.targSamples = targSamples;
        String[] refSampToPanelId = AdmixUtils.sampleMap(refSamples, par.ref_panel());
        this.refPanelIds = AdmixRefPanels.indexPanels(refSampToPanelId);
        this.refHapToPanel = AdmixRefPanels.hapToPanel(refSampToPanelId,
                refPanelIds);
        this.nPanelHaps = AdmixRefPanels.nPanelHaps(refHapToPanel,
                refPanelIds.length);
        this.ancIds = refPanelIds;
        this.globalAncestries = par.gt_ancestries()==null
                ? new GlobalAncestries(ancIds.length, targSamples.size())
                : new GlobalAncestries(par.gt_ancestries(), ancIds, targSamples.ids());
    }

    private static String[] readAncAndPanelLines(File modelFile) {
        String[] lines = new String[2];
        int index = 0;
        try (FileIt<String> it = InputIt.fromGzipFile(modelFile)) {
            while (it.hasNext() && index<2) {
                String candidate = it.next().trim();
                if (candidate.length()>0 && candidate.charAt(0)!='#') {
                    lines[index++] = candidate;
                }
            }
        }
        if (index<2) {
            Utilities.exit("ERROR: Missing lines in model file");
        }
        return lines;
    }

    /**
     * Exits the Java virtual machine with an error message if
     * {@code ancIds.length < 2}
     * @param ancIds the list of ancestry identifiers
     * @throws NullPointerException if {@code ancIds == null}
     */
    public static void confirmMoreThanOneAncestry(String[] ancIds) {
        if (ancIds.length<2) {
            String err = "Fewer than two ancestries are specified";
            String info = Const.nl + "Error       :  " + err
                    + Const.nl     + "Ancestries  :  " + Arrays.toString(ancIds);
            Utilities.exit(new Exception(err), info);
        }
    }

    /**
     * Returns the command line parameters.
     * @return the command line parameters
     */
    public AdmixPar par() {
        return par;
    }

    /**
     * Return the list of reference sample identifiers.
     * @return the list of reference sample identifiers
     */
    public Samples refSamples() {
        return refSamples;
    }

    /**
     * Return the list of target sample identifiers.
     * @return the list of target sample identifiers
     */
    public Samples targSamples() {
        return targSamples;
    }

    /**
     * Returns the number of reference haplotypes.
     * @return the number of reference haplotypes
     */
    public int nRefHaps() {
        return refHapToPanel.size();
    }

    /**
     * Returns the number of reference panels.
     * @return the number of reference panels
     */
    public int nRefPanels() {
        return refPanelIds.length;
    }


    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return ancIds.length;
    }

    /**
     * Returns the specified reference panel identifier.
     * @param refPanel a reference panel index
     * @return the specified reference panel identifier
     * @throws IndexOutOfBoundsException if
     * {@code refPanel < 0 || refPanel >= this.nRefPanels()}
     */
    public String refPanelId(int refPanel) {
        return refPanelIds[refPanel];
    }

    /**
     * Returns the list of reference panel identifiers.  The returned list
     * has length {@code this.nRefPanels()}.
     * @return the list of reference panel identifiers
     */
    public String[] refPanelIds() {
        return refPanelIds.clone();
    }

    /**
     * Returns a map of reference haplotype index to reference panel index.
     * @return a map of reference haplotype index to reference panel index
     */
    public IntArray refHapToPanel() {
        return refHapToPanel;
    }

    /**
     * Returns an {@code IntArray} of length {@code this.nRefPanels()} whose
     * {@code j}-th element is the number of haplotypes in the {@code j}-th
     * reference panel.
     * @return an {@code IntArray} of length {@code this.nRefPanels()} whose
     * {@code j}-th element is the number of haplotypes in the {@code j}-th
     * reference panel
     */
    public IntArray nPanelHaps() {
        return nPanelHaps;
    }

    /**
     * Returns the specified ancestry identifier.
     * @param anc an ancestry index
     * @return the specified ancestry identifier
     * @throws IndexOutOfBoundsException if
     * {@code anc < 0 || anc >= this.nAnc()}
     */
    public String ancId(int anc) {
        return ancIds[anc];
    }

    /**
     * Returns the list of ancestry identifiers.
     * @return the list of ancestry identifiers
     */
    public String[] ancIds() {
        return ancIds.clone();
    }

    /**
     * Returns an array with length {@code this.targSamples().size()}
     * whose {@code j}-th entry is an array with global ancestry proportions
     * for the {@code j}-th target sample. The contract for this method
     * is undefined if the the elements of the specified
     * {@code defaultAncProbs} array are not finite, non-negative numbers.
     * @param defaultAncProbs the global ancestry proportions for samples
     * that do not have user-specified global ancestry proportions.
     * @return the ancestry proportions for target samples
     * @throws NullPointerException if {@code defaultAncProbs == null}
     * @throws IllegalArgumentException if
     * {@code defaultAncProbs.length != this.fixedParams().nAnc()}
     */
    public double[][] globalAncestries(double[] defaultAncProbs) {
        return globalAncestries.globalAncestries(defaultAncProbs);
    }
}

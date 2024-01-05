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
import blbutil.Utilities;
import ints.IntArray;
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import vcf.Samples;

/**
 * <p>Class {@code AdmixRefPanels} has static methods for constructing a map
 * from reference haplotype index to reference panel index.</p>
 *
 * <p>Instances of class {@code AdmixRefPanels} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
final class AdmixRefPanels {

    private AdmixRefPanels() {
        // private constructor to prevent instantiation
    }

    /**
     * Checks the reference panel identifiers.
     * The Java virtual machine will exit with an error message if there
     * are more than {@code Short.MAX_VALUE} distinct reference panel
     * identifiers, if the set of distinct reference panel identifiers
     * in {@code sampToRefPanelId} is not equal to the set of distinct
     * reference panel identifiers in {@code refPanelIds}, or if any
     * two elements of {@code refPanelIds} are identical.
     *
     * @param refSamples the list of distinct reference samples identifiers
     * @param sampToRefPanelId a map from reference sample index to reference
     * panel identifier
     * @param refPanelIds the list of reference panel identifiers
     * @throws IllegalArgumentException if
     * {@code refSampId.length != sampToRefPanelId.length}
     * @throws NullPointerException if any argument is {@code null} or if
     * any array element is {@code null}
     */
    public static void checkRefPanelIds(Samples refSamples,
            String[] sampToRefPanelId, String[] refPanelIds) {
        if (refSamples.size()!=sampToRefPanelId.length) {
            throw new IllegalArgumentException("inconsistent data");
        }
        checkRefPanelCnt(refPanelIds.length);
        HashSet<String> sampPanelIdSet = new HashSet<>();
        HashSet<String> refPanelIdSet = checkForDuplicates(refPanelIds);
        for (int j=0; j<sampToRefPanelId.length; ++j) {
            String id = sampToRefPanelId[j];
            if (sampPanelIdSet.add(id)==true && refPanelIdSet.remove(id)==false) {
                String err = "Reference sample is not assigned to a known reference panel";
                String info = Const.nl + "Error              :  " + err
                        + Const.nl     + "Reference sample   :  " + refSamples.id(j)
                        + Const.nl     + "Assigned ref-panel :  " + id
                        + Const.nl     + "Known ref-panels   :  " + Arrays.toString(refPanelIds);
                Utilities.exit(new Throwable(err), info);
            }
        }
        if (refPanelIdSet.isEmpty()==false) {
            String err = "No reference samples in reference panel";
            String info = Const.nl + "Error            :  " + err
                    + Const.nl     + "Empty ref-panels :  " + refPanelIdSet.toString();
            Utilities.exit(new Throwable(err), info);
        }
    }

    private static HashSet<String> checkForDuplicates(String[] refPanelIds) {
        HashSet<String> panelIdSet = new HashSet<>(refPanelIds.length);
        for (int j=0; j<refPanelIds.length; ++j) {
            if (panelIdSet.add(refPanelIds[j])==false) {
                String err = "Duplicate reference panel in list of reference panels";
                String info = Const.nl + "Error               :  " + err
                        + Const.nl     + "Duplicate ref-panel :  " + refPanelIds[j]
                        + Const.nl     + "Ref-panel list      :  " + Arrays.toString(refPanelIds);
                Utilities.exit(new Throwable(err), info);
            }
        }
        return panelIdSet;
    }

    /**
     * Checks that there are not more than {@code Short.MAX_VALUE} reference
     * panels.  The Java virtual machine will exit with an error message if
     * {@code (nRefPanels > Short.MAX_VALUE)}.
     * @param nRefPanels the number of reference panels
     */
    private static void checkRefPanelCnt(int nRefPanels) {
        if (nRefPanels>Short.MAX_VALUE) {
            String err = "More than " + Short.MAX_VALUE + " reference panels";
            String info = Const.nl + "Error                 :  " + err
                    + Const.nl     + "Number of ref panels  :  " + nRefPanels;
            Utilities.exit(new Throwable(err), info);
        }
    }

    /**
     * Returns a map from reference haplotype index to reference panel index.
     * @param sampToPanelId a map from reference sample index to reference panel
     * identifier
     * @param panelIds the list of distinct reference panel identifiers
     * @return a map from reference haplotype index to reference panel index
     * @throws IllegalArgumentException if any two elements of the
     * {@code panelIds} array are equal
     * @throws NullPointerException if the set of identifiers in the
     * {@code sampToPanelId} array is not a subset of the set of identifiers in
     * the {@code panelIds} array
     * @throws NullPointerException if any argument is {@code null} or if
     * any element of the specified arrays is {@code null}
     */
    public static IntArray hapToPanel(String[] sampToPanelId, String[] panelIds) {
        HashMap<String, Integer> panelToIndexMap = AdmixUtils.inverseMap(panelIds);
        int[] sampToPanelIndex = Arrays.stream(sampToPanelId)
                .parallel()
                .mapToInt(panel -> panelToIndexMap.get(panel))
                .toArray();
        int[] hapToPanelIndex = AdmixUtils.repeatEachElement(sampToPanelIndex);
        return IntArray.create(hapToPanelIndex, panelToIndexMap.size());
    }

    /**
     * Constructs a list of unique reference panel identifiers in the order
     * that they appear.  The constructor exits with an error message if
     * there are more than {@code SHORT.MAX_VALUE} distinct reference panels.
     *
     * @param refSampToPanelId a map from sample index to reference panel
     * identifier
     * @return a list of the unique reference panel identifiers in the order
     * they appear
     * @throws NullPointerException if {@code (refSampToPanelId == null)} or
     * if {@code (0 <= j) && (j < refSampToPanelId.length)
     * && (refSampToPanelId[j] == null)}
    */
    public static String[] indexPanels(String[] refSampToPanelId) {
        String[] refPanelIds = AdmixUtils.uniqueElements(refSampToPanelId);
        checkRefPanelCnt(refPanelIds.length);
        return refPanelIds;
    }

    /**
     * Returns a list which maps reference panel index to the number of
     * haplotypes in the reference panel.
     * @param refHapToPanel a map from reference haplotype index to reference
     * panel index
     * @param nPanels the number of reference panels
     * @return a list which maps reference panel index to the number of
     * haplotypes in the reference panel
     * @throws IndexOutOfBoundsException if
     * {@code (0 <= j && j < refHapToPanel.size() && refHapToPanel.get(j) >= nPanels)}
     * @throws NullPointerException if {@code refHapToPanel == null}
     */
    public static IntArray nPanelHaps(IntArray refHapToPanel, int nPanels) {
        int[] cnts = new int[nPanels];
        for (int j=0, n=refHapToPanel.size(); j<n; ++j) {
            ++cnts[refHapToPanel.get(j)];
        }
        return new WrappedIntArray(cnts);
    }
}

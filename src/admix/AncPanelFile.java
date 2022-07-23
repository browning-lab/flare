/*
 * Copyright 2021 Brian L. Browning
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
import blbutil.StringUtil;
import blbutil.Utilities;
import ints.IntArray;
import ints.WrappedIntArray;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * <p>Class {@code AncPanelFile} reads a stores  a list of ancestries and
 * relevant reference panels for each ancestry.</p>
 *
 * <p>Instances of class {@code AncPanelFile} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AncPanelFile {

    private final String[] ancIds;
    private final IntArray[] ancToRefPanels;

    /**
     * Constructs a new {@code AncPanelFile} instance for the specified
     * data. The Java virtual machine will exit with an error message if any
     * two reference panel identifiers are equal, if any non-blank line
     * in the anc-panel file does not have at least two white-space
     * delimited fields, if there are less than two non-blank lines,
     * if any entry in the first column of the anc-panel file is duplicated,
     * if any field after the first field on a line is not equal to a
     * reference panel identifier, or if any two fields after the first field
     * on a line are identical.
     *
     * @param refPanelIds the list of reference panel identifiers
     * @param ancPanelFile an anc-panel file
     * @throws NullPointerException if
     * {@code (refPanelIdst == null || ancPanelFile == null)} or if
     * {@code (0 <= j && j < refPanelIds.length && refPanelIds[j] == null)}
     */
    public AncPanelFile(String[] refPanelIds, File ancPanelFile) {
        checkForDuplicateRefPanel(refPanelIds);
        HashMap<String, Integer> panelToIndexMap = AdmixUtils.inverseMap(refPanelIds);
        ArrayList<String> ancList = new ArrayList<>();
        HashSet<String> ancSet = new HashSet<>();
        ArrayList<IntArray> panelLists = new ArrayList<>();
        String[] lines = AdmixUtils.readLines(ancPanelFile);
        for (String line : lines) {
            extractAncPanelData(ancPanelFile, line, ancList, ancSet,
                    panelLists, panelToIndexMap);
        }
        this.ancIds = ancList.toArray(new String[0]);
        this.ancToRefPanels = panelLists.toArray(new IntArray[0]);
        confirmMoreThanOneAncestry(ancIds);
    }

    private static void extractAncPanelData(File ancPanelFile, String line,
            ArrayList<String> ancList, HashSet<String> ancSet,
            ArrayList<IntArray> panelLists,
            HashMap<String, Integer> panelToIndexMap) {
        String[] fields = splitAndCheckFieldCnt(ancPanelFile, line);
        if (fields.length > 0) {
            boolean nonDuplicateAncestry = ancSet.add(fields[0]);
            if (nonDuplicateAncestry==false) {
                dupAncestryError(ancPanelFile, fields[0]);
            } else {
                ancList.add(fields[0]);
                int[] indices = new int[fields.length - 1];
                for (int j=1; j<fields.length; ++j) {
                    Integer index = panelToIndexMap.get(fields[j]);
                    if (index==null) {
                        unknownPanelError(ancPanelFile, line, fields[j]);
                    } else {
                        indices[j-1] = index;
                    }
                }
                Arrays.sort(indices);
                checkForRedundantPanel(ancPanelFile, line, fields, indices);
                panelLists.add(new WrappedIntArray(indices));
            }
        }
    }

    private static void checkForDuplicateRefPanel(String[] refPanels) {
        HashSet<String> idSet = new HashSet<>(refPanels.length);
        for (String id : refPanels) {
            if (idSet.add(id)==false) {
                String err = "Duplicate reference panel identifier";
                String info = Const.nl + "Error            :  " + err
                        + Const.nl     + "Reference panel  :  " + id;
                Utilities.exit(new Exception(err), info);
            }
        }
    }

    private static String[] splitAndCheckFieldCnt(File ancPanelFile, String line) {
        String[] fields = StringUtil.getFields(line);
        if (fields.length == 1) {
            String err = "Line in anc-panel file has only one field";
            String info = Const.nl + "Error           :  " + err
                    + Const.nl     + "Anc-panel file  :  " + ancPanelFile
                    + Const.nl     + "Line            :  " + line;
            Utilities.exit(new Exception(err), info);
        }
        return fields;
    }

    private static void dupAncestryError(File ancPanelFile, String ancestry) {
        String err = "Duplicate ancestry identifiers in first column of anc-panel file";
        String info = Const.nl + "Error           :  " + err
                + Const.nl     + "Anc-panel file  :  " + ancPanelFile
                + Const.nl     + "Ancestry        :  " + ancestry;
        Utilities.exit(new Exception(err), info);
    }

    private static void unknownPanelError(File ancPanelFile, String line,
            String panel) {
        String err = "Unrecognized reference panel in anc-panel file";
        String info = Const.nl + "Error           :  " + err
                + Const.nl     + "Anc-panel file  :  " + ancPanelFile
                + Const.nl     + "Reference panel :  " + panel
                + Const.nl     + "Line            :  " + line;
        Utilities.exit(new Exception(err), info);
    }

    private static void checkForRedundantPanel(File file, String line,
            String[] fields, int[] sortedIndices) {
        for (int j=1; j<sortedIndices.length; ++j) {
            if (sortedIndices[j] == sortedIndices[j-1]) {
                String err = "Duplicate reference panel associated with ancestry";
                String info = Const.nl + "Error:          :  " + err
                        + Const.nl     + "Anc-panel file  :  " + file
                        + Const.nl     + "Ancestry        :  " + fields[0]
                        + Const.nl     + "Reference panel :  " + fields[j]
                        + Const.nl     + "Line            :  " + line;
                Utilities.exit(new Exception(err), info);
            }
        }
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
     * Returns the list of ancestry identifiers.t
     * @return the list of ancestry identifiers
     */
    public String[] ancIds() {
       return ancIds.clone();
    }

    /**
     * Returns an array of length {@code this.nAnc()} containing
     * a sorted, increasing list of reference panel indices associated
     * with each ancestry.
     * @return the reference panels associated with each ancestry
     */
    public IntArray[] ancToRefPanels() {
        return ancToRefPanels.clone();
    }
}

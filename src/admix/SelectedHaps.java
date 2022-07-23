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

import blbutil.Utilities;
import ints.IntArray;
import ints.IntList;
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * <p>Class {@code SelectedHaps} contains selected target and reference
 * haplotypes.</p>
 *
 * <p>Instances of {@code SelectedHaps} are immutable</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SelectedHaps {

    private final FixedParams fixedParams;
    private final boolean includeRefHaps;
    private final WrappedIntArray hapList;
    private final WrappedIntArray[] panelToHapListIndices;

    /**
     * Constructs a {@code SelectedHaps} object from the specified data.
     * All target haplotypes are selected.  If {@code includeRefHaps == true},
     * then {@code fixedParams.par().em_haps()} haplotypes are randomly
     * selected from each reference panel.  If {@code includeRefHaps == true}
     * and the reference panel contains fewer than
     * {@code fixedParams.par().em_haps()} haplotypes, then all haplotypes
     * in the reference panel are selected.
     *
     * @param fixedParams the fixed parameters.
     * @param includeRefHaps {@code true} if the selected haplotypes should
     * include reference haplotypes
     * @throws NullPointerException if {@code fixedParams == null}
     */
    public SelectedHaps(FixedParams fixedParams, boolean includeRefHaps) {
        int nTargSamples = fixedParams.targSamples().size();
        int nTargHaps =  nTargSamples << 1;
        this.fixedParams = fixedParams;
        this.includeRefHaps = includeRefHaps;
        WrappedIntArray panelToRefHaps[] = panelToShiftedRefHaps(fixedParams,
                includeRefHaps);
        this.hapList = hapList(nTargHaps, panelToRefHaps);
        this.panelToHapListIndices = refPanelToHapListIndices(nTargHaps,
                panelToRefHaps);
    }

    /**
     * Returns a {@code SelectedHaps} instance that is obtained by removing all
     * reference haplotypes from {@code this}.
     * @return a {@code SelectedHaps} instance that is obtained by removing all
     * reference haplotypes from {@code this}
     */
    public SelectedHaps removeRefHaps() {
        if (includeRefHaps) {
            boolean includeRefHaplotypes = false;
            return new SelectedHaps(fixedParams, includeRefHaplotypes);
        }
        else {
            return this;
        }
    }

    private static WrappedIntArray[] panelToShiftedRefHaps(FixedParams fixedParams,
            boolean includeRefHaps) {
        int nRefPanels = fixedParams.nRefPanels();
        if (includeRefHaps) {
            // ref haplotype indices are shifted by nTargHaps;
            long seed = fixedParams.par().seed();
            int maxHaps = fixedParams.par().em_haps();
            IntList[] panelToHaps = panelToShiftedRefHaps(fixedParams);
            return IntStream.range(0, nRefPanels)
                    .parallel()
                    .mapToObj(j -> randomSubset(panelToHaps[j], maxHaps, (seed+j)))
                    .toArray(WrappedIntArray[]::new);
        } else {
            IntArray ia = new WrappedIntArray(new int[0]);
            return IntStream.range(0, nRefPanels)
                    .mapToObj(j -> ia)
                    .toArray(WrappedIntArray[]::new);
        }
    }

    private static IntList[] panelToShiftedRefHaps(FixedParams fixedParams) {
        // ref haplotype indices are shifted by nTargHaps;
        int nTargSamples = fixedParams.targSamples().size();
        int nTargHaps = nTargSamples << 1;
        IntList[] panelToShiftedRefHaps = IntStream.range(0, fixedParams.nRefPanels())
                .parallel()
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        IntArray refHapToPanel = fixedParams.refHapToPanel();
        for (int h=0, n=refHapToPanel.size(); h<n; ++h) {
            int panel = refHapToPanel.get(h);
            int shiftedRefHap = nTargHaps + h;
            panelToShiftedRefHaps[panel].add(shiftedRefHap);
        }
        return panelToShiftedRefHaps;
    }

    private static IntArray randomSubset(IntList list, int size, long seed) {
        Random rand = new Random(seed);
        int[] ia = list.toArray();
        if (size<ia.length) {
            Utilities.shuffle(ia, size, rand);
            ia = Arrays.copyOf(ia, size);
        }
        return new WrappedIntArray(ia);
    }

    private static WrappedIntArray hapList(int nTargHaps,
            IntArray[] panelToRefHaps) {
        IntList hapList = new IntList(nTargHaps);
        for (int h=0; h<nTargHaps; ++h) {
            hapList.add(h);
        }
        for (int j=0; j<panelToRefHaps.length; ++j) {
            IntArray refHaps = panelToRefHaps[j];
            for (int h=0, n=refHaps.size(); h<n; ++h) {
                hapList.add(refHaps.get(h));
            }
        }
        return new WrappedIntArray(hapList);
    }

    private static WrappedIntArray[] refPanelToHapListIndices(int nTargHaps,
            WrappedIntArray[] panelToRefHaps) {
        int nPanels = panelToRefHaps.length;
        WrappedIntArray[] panelToHapListIndices = new WrappedIntArray[nPanels];
        int offset = nTargHaps;
        for (int j=0; j<nPanels; ++j) {
            int[] indices = new int[panelToRefHaps[j].size()];
            for (int h=0; h<indices.length; ++h) {
                indices[h] = offset + h;
            }
            panelToHapListIndices[j] = new WrappedIntArray(indices);
            offset += indices.length;
        }
        return panelToHapListIndices;
    }

    /**
     * Returns the fixed parameters.
     *
     * @return the fixed parameters
     */
    public FixedParams fixedParams() {
        return fixedParams;
    }

    /**
     * Returns {@code true} if the selected haplotypes include
     * reference haplotypes, and returns {@code false} otherwise.
     *
     * @return {@code true} if the selected haplotypes include  reference
     * haplotypes
     */
    public boolean includeRefHaps() {
        return includeRefHaps;
    }

    /**
     * Returns the list of selected haplotypes.
     * @return the list of selected haplotypes
     */
    public WrappedIntArray hapList() {
        return hapList;
    }

    /**
     * Returns a list of indices of the {@code this.hapList()} elements
     * that are reference haplotypes from the specified reference panel.
     * @param refPanel a reference panel index
     * @return a list of indices of the {@code this.hapList()} elements
     * that are reference haplotypes from the specified reference panel
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || (refPanel >= this.fixedParams().nRefPanels())}
     */
    public WrappedIntArray hapListIndices(int refPanel) {
        return panelToHapListIndices[refPanel];
    }
}

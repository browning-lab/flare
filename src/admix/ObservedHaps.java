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

import blbutil.Utilities;
import ints.IntArray;
import ints.IntList;
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * <p>Class {@code ObservedHaps} contains a set of target and reference
 * haplotypes that are used to estimate model parameters.</p>
 *
 * <p>Instances of {@code ObservedHaps} are immutable</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ObservedHaps {

    private final SampleData sampleData;
    private final boolean includeRefHaps;
    private final WrappedIntArray observedHaps;
    private final WrappedIntArray[] panelToObservedHapsIndices;

    /**
     * Constructs a {@code ObservedHaps} object from the specified data.
     * All target haplotypes are selected.  If {@code includeRefHaps == true},
     * then {@code sampleData.par().em_haps()} haplotypes are randomly
     * selected from each reference panel.  If {@code includeRefHaps == true}
     * and the reference panel contains fewer than
     * {@code sampleData.par().em_haps()} haplotypes, then all haplotypes
     * in the reference panel are selected.
     *
     * @param sampleData reference and target sample metadata
     * @param includeRefHaps {@code true} if the selected haplotypes should
     * include reference haplotypes
     * @throws NullPointerException if {@code (sampleData == null)}
     */
    public ObservedHaps(SampleData sampleData, boolean includeRefHaps) {
        int nTargSamples = sampleData.targSamples().size();
        int nTargHaps =  nTargSamples << 1;
        this.sampleData = sampleData;
        this.includeRefHaps = includeRefHaps;
        WrappedIntArray[] panelToSelectedRefHaps = panelToSelectedRefHaps(
                sampleData, includeRefHaps);
        this.observedHaps = observedHaps(nTargHaps, panelToSelectedRefHaps);
        this.panelToObservedHapsIndices = panelToObservedHapsIndices(nTargHaps,
                panelToSelectedRefHaps);
    }

    /**
     * Returns a {@code ObservedHaps} instance that is obtained by removing all
     * reference haplotypes from {@code this}.
     * @return a {@code ObservedHaps} instance that is obtained by removing all
     * reference haplotypes from {@code this}
     */
    public ObservedHaps removeRefHaps() {
        if (includeRefHaps) {
            boolean includeRefHaplotypes = false;
            return new ObservedHaps(sampleData, includeRefHaplotypes);
        }
        else {
            return this;
        }
    }

    private static WrappedIntArray[] panelToSelectedRefHaps(SampleData sampleData,
            boolean includeRefHaps) {
        int nRefPanels = sampleData.nRefPanels();
        if (includeRefHaps) {
            // ref haplotype indices are shifted by nTargHaps;
            long seed = sampleData.par().seed();
            int maxHaps = sampleData.par().em_haps();
            IntList[] panelToHaps = panelToShiftedRefHaps(sampleData);
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

    private static IntList[] panelToShiftedRefHaps(SampleData sampleData) {
        // ref haplotype indices are shifted by nTargHaps;
        int nTargSamples = sampleData.targSamples().size();
        int nTargHaps = nTargSamples << 1;
        IntList[] panelToShiftedRefHaps = IntStream.range(0, sampleData.nRefPanels())
                .parallel()
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        IntArray refHapToPanel = sampleData.refHapToPanel();
        for (int h=0, n=refHapToPanel.size(); h<n; ++h) {
            int panel = refHapToPanel.get(h);
            int shiftedRefHap = nTargHaps + h;
            panelToShiftedRefHaps[panel].add(shiftedRefHap);
        }
        return panelToShiftedRefHaps;
    }

    /**
     * Returns a sorted list of {@code size} elements randomly selected
     * without repetition from the specified list. The original list is
     * returned if the specified list has fewer than {@code size} elements.
     * @param list a list of integers
     * @param size the maximum number of elements in the returned list
     * @param seed a seed for random number generation
     * @return a list of elements randomly selected from the specified list
     * @throws NullPointerException if {@code list == null}
     */
    public static IntArray randomSubset(IntList list, int size, long seed) {
        if (size<list.size()) {
            Random rand = new Random(seed);
            int[] ia = list.toArray();
            Utilities.shuffle(ia, size, rand);
            Arrays.sort(ia, 0, size);
            return new WrappedIntArray(Arrays.copyOf(ia, size));
        }
        else {
            return new WrappedIntArray(list);
        }
    }

    private static WrappedIntArray observedHaps(int nTargHaps,
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

    private static WrappedIntArray[] panelToObservedHapsIndices(int nTargHaps,
            WrappedIntArray[] panelToSelectedRefHaps) {
        int nPanels = panelToSelectedRefHaps.length;
        WrappedIntArray[] panelToHapListIndices = new WrappedIntArray[nPanels];
        int offset = nTargHaps;
        for (int j=0; j<nPanels; ++j) {
            int[] indices = new int[panelToSelectedRefHaps[j].size()];
            for (int h=0; h<indices.length; ++h) {
                indices[h] = offset + h;
            }
            panelToHapListIndices[j] = new WrappedIntArray(indices);
            offset += indices.length;
        }
        return panelToHapListIndices;
    }

    /**
     * Returns the reference and target sample metadata.
     *
     * @return the reference and target sample metadata
     */
    public SampleData sampleData() {
        return sampleData;
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
        return observedHaps;
    }

    /**
     * Returns a list of indices of {@code this.observedHaps()} elements
     * that are reference haplotypes from the specified reference panel.
     * @param refPanel a reference panel index
     * @return a list of indices of {@code this.observedHaps()} elements
     * that are reference haplotypes from the specified reference panel
     * @throws IndexOutOfBoundsException if
     * {@code (refPanel < 0) || (refPanel >= this.sampleData().nRefPanels())}
     */
    public WrappedIntArray panelToObservedHapsIndices(int refPanel) {
        return panelToObservedHapsIndices[refPanel];
    }
}

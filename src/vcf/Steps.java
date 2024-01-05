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
package vcf;

import blbutil.DoubleArray;
import ints.IntList;

/**
 * <p>Class {@code Steps} represents a partition of a list of markers into
 * a sequence of sets of consecutive markers (the steps).</p>
 *
 * <p>Instances of class {@code Steps} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Steps {

    private final MarkerMap map;
    private final int[] stepEnds;

    /**
     * Returns the estimated number of bytes consumed by this object,
     * excluding the overhead bytes required by {@code this}.
     * @return the estimated number of bytes required to store this object
     */    
    public long estBytes() {
        // assume this.map is not owned
        int overhead = 12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 2*8;  // assume 8 bytes per reference
        estBytes += 4 + 4*stepEnds.length; // includes 4 bytes for array length
        return estBytes;
    }

    /**
     * Constructs a new {@code Steps} instance from the specified data.
     * @param map the marker map
     * @param minStep the minimum distance between the first markers in
     * consecutive steps
     *
     * @throws IllegalArgumentException if
     * {@code minStep <= 0f || Float.isFinite(minStep) == false}
     * @throws NullPointerException if {@code map == null}
     */
    public Steps(MarkerMap map, float minStep) {
        if (minStep <= 0f || Float.isFinite(minStep)==false) {
            throw new IllegalArgumentException(String.valueOf(minStep));
        }
        this.map = map;
        this.stepEnds = stepEnds(map, minStep);
    }

     private static int[] stepEnds(MarkerMap map, double minStep) {
        DoubleArray genPos = map.genPos();
        int nMarkers = genPos.size();
        IntList indices = new IntList(nMarkers>>1);
        int end = 0;
        while (end<nMarkers) {
            double minGenPos = genPos.get(end) + minStep;
            ++end;
            while (end<genPos.size() && genPos.get(end)<minGenPos) {
                ++end;
            }
            indices.add(end);
        }
        return indices.toArray();
    }

    /**
     * Returns the number of steps.
     * @return the number of steps
     */
    public int size() {
        return stepEnds.length;
    }

    /**
     * Returns the index of the first marker in the specified step.
     * @param step a step index
     * @return the index of the first marker in the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int start(int step) {
        return step==0 ? 0 : stepEnds[step-1];
    }

    /**
     * Returns the index of the last marker (exclusive) in the specified step.
     * @param step a step index
     * @return the index of the last marker (exclusive) in the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int end(int step) {
        return stepEnds[step];
    }

    /**
     * Return the marker map.
     * @return the marker map
     */
    public MarkerMap map() {
        return map;
    }
}

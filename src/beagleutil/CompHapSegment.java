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
package beagleutil;

/**
 * <p>Class {@code CompHapSegment} represents a copied haplotype segment
 * in a composite reference haplotype.</p>
 *
 * <p>Instances of class {@code CompHapSegment} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CompHapSegment implements Comparable<CompHapSegment> {

    private int hap;
    private int startMarker;
    private int lastIbsStep;
    private final int compHapIndex;

    /**
     * Constructs a new {@code CompHapSegment} from the specified data.
     * @param hap the haplotype
     * @param startMarker the index of the first marker in the haplotype segment
     * @param ibsStep the last recorded IBS step
     * @param compHapIndex the composite haplotype index
     */
    public CompHapSegment(int hap, int startMarker, int ibsStep, int compHapIndex) {
        this.hap = hap;
        this.startMarker = startMarker;
        this.lastIbsStep = ibsStep;
        this.compHapIndex = compHapIndex;
    }

    /**
     * Update the haplotype, the first marker in the haplotype segment,
     * and the last recorded IBS step.
     * @param hap the haplotype
     * @param startMarker the first marker in the haplotype segment
     * @param lastIbsStep the last recorded IBS step
     */
    public void updateSegment(int hap, int startMarker, int lastIbsStep) {
        this.hap = hap;
        this.startMarker = startMarker;
        this.lastIbsStep = lastIbsStep;
    }

    /**
     * Updates the last recorded IBS step to the specified value
     * @param ibsStep the last recorded IBS Step
     */
    public void setLastIbsStep(int ibsStep) {
        this.lastIbsStep = ibsStep;
    }

    /**
     * Returns the haplotype.
     * @return the haplotype
     */
    public int hap() {
        return hap;
    }

    /**
     * Returns the first marker in the haplotype segment
     * @return the first marker in the haplotype segment
     */
    public int startMarker() {
        return startMarker;
    }

    /**
     * Returns the last recorded IBS step for {@code this.hap()}.
     * @return the last recorded IBS step for {@code this.hap()}
     */
    public int lastIbsStep() {
        return lastIbsStep;
    }

    /**
     * Returns the composite haplotype index.
     * @return the composite haplotype index
     */
    public int compHapIndex() {
        return compHapIndex;
    }

    /**
     * Compares the specified segment to {@code this} for order.  Returns
     * -1, 0, or 1 according to whether {@code this.end()} is less than,
     * equal, or greater than {@code seg.end()}.
     * @param seg the object to be compared
     * @return -1, 0, or 1 according to whether {@code this.end()} is less
     * than, equal, or greater than {@code seg.end()}
     */
    @Override
    public int compareTo(CompHapSegment seg) {
        if (this.lastIbsStep!=seg.lastIbsStep) {
            return this.lastIbsStep<seg.lastIbsStep ? -1 : 1;
        } else {
            return 0;
        }
    }

}

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

import beagleutil.ChromIds;

/**
 * <p>Class {@code PositionMap} represents a genetic map obtained by
 * multiplying chromosome position by a scale factor.
 * </p>
 * <p>Instances of class {@code PositionMap} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PositionMap implements GeneticMap {

    private final double scaleFactor;
    private final double invScaleFactor;

    /**
     * Returns the scale factor that is multiplied by the chromosome position
     * to obtain the corresponding genetic map position
     * @return the scale factor.
     */
    public double scaleFactor() {
        return scaleFactor;
    }

    /**
     * Constructs a new {@code PositionMap} instance.
     * @param scaleFactor the factor that is multiplied by
     * a base position to obtain the corresponding genetic map position
     * @throws IllegalArgumentException if
     * {@code scaleFactor <= 0d || Double.isFinite(scaleFactor) == false}
     */
    public PositionMap(double scaleFactor) {
        if (Double.isFinite(scaleFactor) == false || scaleFactor <= 0d) {
            throw new IllegalArgumentException(String.valueOf(scaleFactor));
        }
        this.scaleFactor = scaleFactor;
        this.invScaleFactor = 1.0/scaleFactor;
    }

    @Override
    public int basePos(int chrom, double geneticPosition) {
        if (chrom < 0 || chrom >= ChromIds.instance().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(chrom));
        }
        long pos = Math.round(geneticPosition*invScaleFactor);
        if (pos>Integer.MAX_VALUE) {
            String s = "An estimated base position exceeds the maximum integer value"
                    + blbutil.Const.nl
                    + "Is the window parameter in cM units?";
            throw new IllegalArgumentException(s);
        }
        return (int) pos;
    }

    @Override
    public double genPos(Marker marker) {
        return scaleFactor*marker.pos();
    }

    @Override
    public double genPos(int chrom, int basePosition) {
        if (chrom < 0 || chrom >= ChromIds.instance().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(chrom));
        }
        return scaleFactor*basePosition;
    }
}

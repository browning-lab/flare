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
import blbutil.FloatArray;
import ints.IntArray;
import java.util.stream.IntStream;

/**
 * <p>Class {@code MarkerRecombMap} represents genetic map positions and
 * inter-marker genetic distances for a sequence of genomic loci.
 * </p>
 * <p>Instances of class {@code MarkerRecombMap} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MarkerMap {

    private final DoubleArray genPos;
    private final FloatArray genDist;

    /**
     * Returns the estimated number of bytes consumed by this object,
     * excluding the overhead bytes required by {@code this}.
     * @return the estimated number of bytes required to store this object
     */
    public long estBytes() {
        int overhead = 2*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 2*8;  // assume 8 bytes per reference
        estBytes += 4 + 8*genPos.size(); // includes 4 bytes for array length
        estBytes += 4 + 4*genDist.size(); // includes 4 bytes for array length
        return estBytes;
    }

    /**
     * Returns a new {@code MarkerMap} instance constructed
     * from the specified data.
     * @param genMap the genetic map
     * @param markers a list of markers
     * @return a returns new {@code MarkerMap} instance
     * @throws IllegalArgumentException if
     * {@code markers.marker(0).chromIndex() != markers.marker(markers.nMarkers() - 1).chromIndex()}
     * @throws IllegalArgumentException if the specified genetic map has no
     * map positions for the specified chromosome
     * @throws NullPointerException if
     * {@code genMap == null || markers == null}
     */
    public static MarkerMap create(GeneticMap genMap, Markers markers) {
        double meanGenDiff = MarkerMap.meanSingleBaseGenDist(genMap, markers);
        return new MarkerMap(GeneticMap.genPos(genMap, meanGenDiff, markers));
    }

   /**
     * Returns a new {@code MarkerMap} instance constructed
     * from the specified data.
     * @param genMap the genetic map
     * @param minGenDist the required minimum cM distance between successive
     * markers
     * @param markers a list of markers
     * @return a returns new {@code MarkerMap} instance
     * @throws IllegalArgumentException if
     * {@code markers.marker(0).chromIndex() != markers.marker(markers.nMarkers() - 1).chromIndex()}
     * @throws IllegalArgumentException if {@code Double.isFinite(minDist) == false}
     * @throws IllegalArgumentException if the specified genetic map has no
     * map positions for the specified chromosome
     * @throws NullPointerException if
     * {@code genMap == null || markers == null}
     */
    public static MarkerMap create(GeneticMap genMap, double minGenDist,
            Markers markers) {
        return new MarkerMap(GeneticMap.genPos(genMap, minGenDist, markers));
    }

    /**
     * Returns the mean genetic distance between two consecutive base positions.
     * @param genMap the genetic map
     * @param markers a list of markers
     * @return the mean genetic distance between two consecutive base positions
     * @throws IllegalArgumentException if
     * {@code markers.marker(0).chromIndex() != markers.marker(markers.nMarkers() - 1).chromIndex()}
     * @throws IllegalArgumentException if
     * {@code markers.marker(0).pos() == markers.marker(markers.nMarkers() - 1).pos()}
     * @throws IllegalArgumentException if the specified genetic map has no
     * map positions for the specified chromosome
     * @throws NullPointerException if {@code genMap == null || markers == null}
     */
    public static double meanSingleBaseGenDist(GeneticMap genMap, Markers markers) {
        Marker a = markers.marker(0);
        Marker b = markers.marker(markers.size()-1);
        if (a.chromIndex()!=b.chromIndex()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (a.pos()==b.pos()) {
            String s = "Window has only one position: CHROM=" + a.chromID() + " POS=" + a.pos();
            throw new IllegalArgumentException(s);
        }
        double meanSingleBaseDist = Math.abs(genMap.genPos(b) - genMap.genPos(a))
                / Math.abs(b.pos() - a.pos());
        // require meanSingleBaseDist to be >= 0.01 * mean human single base genetic distance
        return Math.max(meanSingleBaseDist, 1e-8);
    }

    private MarkerMap(double[] gPos) {
        this.genPos = new DoubleArray(gPos);
        this.genDist = genDist(gPos);
    }

    /**
     * Return a marker map restricted to the specified markers
     * @param indices a list of distinct marker indices in increasing order
     * @return a marker map restricted to the specified markers
     * @throws IndexOutOfBoundsException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (indices[j] < 0 || indices[j] >= this.nMarkers())}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (1 <= j && j < indices.length)} such that
     * {@code (indices[j] <= indice[j - 1])}
     * @throws NullPointerException if {@code indices == null}
     */
    public MarkerMap restrict(int[] indices) {
        double[] gPos = new double[indices.length];
        gPos[0] = genPos.get(indices[0]);
        for (int j=1; j<indices.length; ++j) {
            if (indices[j] <= indices[j-1]) {
                throw new IllegalArgumentException(String.valueOf(indices[j]));
            }
            gPos[j] = genPos.get(indices[j]);
        }
        return new MarkerMap(gPos);
    }

    /**
     * Return a marker map restricted to the specified markers
     * @param indices a list of distinct marker indices in increasing order
     * @return a marker map restricted to the specified markers
     * @throws IndexOutOfBoundsException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (indices.get(j) < 0 || indices.get(j) >= this.nMarkers())}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (1 <= j && j < indices.length)} such that
     * {@code (indices.get(j) <= indice.get(j - 1))}
     * @throws NullPointerException if {@code indices == null}
     */
    public MarkerMap restrict(IntArray indices) {
        double[] gPos = new double[indices.size()];
        gPos[0] = genPos.get(indices.get(0));
        for (int j=1, n=indices.size(); j<n; ++j) {
            if (indices.get(j) <= indices.get(j-1)) {
                throw new IllegalArgumentException(String.valueOf(indices.get(j)));
            }
            gPos[j] = genPos.get(indices.get(j));
        }
        return new MarkerMap(gPos);
    }

    /**
     * Returns an array whose {@code (k+1}-st element is the genetic
     * distance between the {@code (k+1)}-st genetic position and
     * the {@code k}-th genetic position.  The first element of the returned
     * array is {@code 0f}.
     * @param genPos an array of strictly increasing genetic positions
     * @return an array whose {@code (k+1}-st element is the genetic
     * distance between the {@code (k+1)}-st genetic position and
     * the {@code k}-th genetic position.
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (0 < j) && (j < genPos.length) && (genPos[j-1] >= genPos[j])}
     * @throws NullPointerException if {@code (genPos == null)}
     */
    private static FloatArray genDist(double[] genPos) {
        float[] da = new float[genPos.length];
        for (int j=1; j<da.length; ++j) {
            da[j] = (float) (genPos[j] - genPos[j-1]);
            if (da[j]<=0) {
                String s = "Nonpositive genetic distance: dist[" + j + "]="
                        + da[j];
                throw new IllegalArgumentException(s);
            }
        }
        return new FloatArray(da);
    }

    /**
     * Returns a {@code DoubleArray} of size {@code this.markers().nMarkers()}
     * whose {@code k}-th element is the genetic map position of the
     * {@code k}-th marker.
     * @return the array of genetic map positions
     */
    public DoubleArray genPos() {
        return genPos;
    }

    /**
     * Return a {@code FloatArray} of size {@code this.markers().nMarkers()}
     * whose {@code k}-th element is the genetic distance between the
     * {@code k}-th target marker and the previous marker, or {@code 0.0}
     * if {@code (k == 0)}.
     * @return a {@code FloatArray} of size {@code this.nTargMarkers()}
     * whose {@code k}-th element is the genetic distance between the
     * {@code k}-th target marker and the previous marker,
     */
    public FloatArray genDist() {
        return genDist;
    }

    /**
     * Returns a map of marker index to the probability of recombination
     * in the interval between the marker and the preceding marker.
     * @param recombIntensity the intensity of the exponential distribution
     * that gives the probability of transitioning to a random HMM state
     * in a specified cM distance
     * @return a map of marker index to the probability of recombination
     * in the interval between the marker and the preceding marker
     * @throws IllegalArgumentException if
     * {@code intensity <= 0.0 || Float.isFinite(intensity)==false}
     */
    public FloatArray pRecomb(float recombIntensity) {
        if (recombIntensity <= 0.0 || Float.isFinite(recombIntensity)==false) {
            throw new IllegalArgumentException(String.valueOf(recombIntensity));
        }
        double c = -recombIntensity;
        double[] pRecomb = IntStream.range(0, genDist.size())
                .parallel()
                .mapToDouble(m -> -Math.expm1(c*genDist.get(m)))
                .toArray();
        return new FloatArray(pRecomb);
    }
}

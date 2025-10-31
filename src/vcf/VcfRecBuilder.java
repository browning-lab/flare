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

import blbutil.Const;
import java.io.PrintWriter;

/**
 * <p>Class {@code VcfRecBuilder} contains methods for constructing
 * and printing a VCF record in VCF 4.2 format.  The FORMAT field data
 * for each sample is added sequentially to the record via the
 * {@code addSampleData()} method.
 *
 * </p>
 * <p>Instances of class {@code VcfRecBuilder} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRecBuilder {

    /**
     * The default initial size for the string buffer, which is 50
     * characters.
     */
    public static final int DEFAULT_INIT_SIZE = 50;

    private final StringBuilder sb;
    private final Marker marker;
    private final int nAlleles;

    /**
     * Constructs a new {@code VcfRecBuilder} instance with an initial
     * capacity for the specified number of samples.
     *
     * @param marker the marker
     * @param nSamples the number of samples
     * @throws IllegalArgumentException if {@code nSamples < 0}
     */
    public VcfRecBuilder(Marker marker, int nSamples) {
        if (nSamples < 0) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        StringBuilder buffer = new StringBuilder(128 + (nSamples<<2));
        MarkerUtils.appendFirst8Fields(marker, buffer);
        buffer.append("\tGT");              // FORMAT
        this.marker = marker;
        this.nAlleles = marker.nAlleles();
        this.sb = buffer;
    }

    /**
     * Returns the current marker.  Returns {@code null} if
     * {@code this.reset()} has not been previously invoked.
     * @return the current marker.
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Adds the specified diploid phased genotype to the VCF record for
     * the current marker.
     * @param a1 the first allele
     * @param a2 the second allele
     * @throws IllegalStateException if {@code this.marker() == null}
     * @throws IndexOutOfBoundsException if
     * {@code a1 < 0 || a1 >= this.marker().nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code a2 < 0 || a2 >= this.marker().nAlleles()}
     */
    public void addSampleData(int a1, int a2) {
        if (a1 < 0 || a1 >= nAlleles) {
            throw new IndexOutOfBoundsException(String.valueOf(a1));
        }
        if (a2 < 0 || a2 >= nAlleles) {
            throw new IndexOutOfBoundsException(String.valueOf(a2));
        }
        sb.append(Const.tab);
        sb.append(a1);
        sb.append(Const.phasedSep);
        sb.append(a2);
    }

    /**
     * Adds the specified haploid phased genotype to the VCF record for
     * the current marker.
     * @param allele the allele
     * @throws IllegalStateException if {@code this.marker() == null}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker().nAlleles()}
     */
    public void addSampleData(int allele) {
        if (allele < 0 || allele >= nAlleles) {
            throw new IndexOutOfBoundsException(String.valueOf(allele));
        }
        sb.append(Const.tab);
        sb.append(allele);
    }

    /**
     * Prints the current VCF record for the current marker to the specified
     * {@code PrintWriter}.
     * @param out the {@code PrintWriter} to which the VCF record will be
     * printed
     * @throws NullPointerException if {@code out == null}
     */
    public void writeRec(PrintWriter out) {
        out.println(sb);
    }
}

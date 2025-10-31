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

import ints.IntArray;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code RefGT} stores a list of samples and a
 * haplotype pair for each sample.
 * </p>
 * <p>Instances of class {@code BasicRefGT} are immutable.<p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefGT implements GT {

    private final Markers markers;
    private final Samples samples;
    private final RefGTRec[] recs;

    @Override
    public long estBytes() {
        int overhead = (1 + recs.length)*12;             // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + (3 + recs.length)*8;  // assume 8 bytes per reference
        estBytes += 4; // include 4 bytes to store array length
        long alleleCodedCnt = 0;
        long seqCodedCnt0 = 0;
        long seqCodedCnt = 0;
        long alleleCodedBytes = 0;
        long seqCodedBytes0 = 0;
        long seqCodedBytes1 = 0;
        long sumLength0 = 0;
        long sumLength1 = 0;
        IntArray[] lastMaps = new IntArray[2];
        for (RefGTRec rec : recs) {
            if (rec.isAlleleRecord()) {
                ++alleleCodedCnt;
                alleleCodedBytes += rec.estBytes();
            }
            else {
                // Cannot delegate to rec.estBytes() because there can be
                // multiple references to the same maps[0] array
                overhead = 2*12;            // assume 12 bytes overhead per "owned" object
                estBytes = overhead + 4*8;  // assume 8 bytes per reference
                IntArray[] maps = rec.maps();
                assert maps.length==2;
                if (maps[0]!=lastMaps[0]) {
                    ++seqCodedCnt0;
                    sumLength0 += maps[0].size();
                    seqCodedBytes0 += (1 + maps[0].size())*4;   // includes 4 bytes to store array length
                }
                ++seqCodedCnt;
                sumLength1 += maps[1].size();
                seqCodedBytes1 += (1 + (maps[1].size()>>3))*4;  // assume 1 bit per allele plus 4 bytes for array length
                lastMaps = maps;
            }
        }
        long otherBytes = estBytes;
        estBytes += (seqCodedBytes0 + seqCodedBytes1 + alleleCodedBytes);
        boolean printSummary = false;
        if (printSummary) {
            System.out.println("RefGT:57 markers:        " + recs.length);
            System.out.println("RefGT:58 alleleCodedCnt: " + alleleCodedCnt);
            System.out.println("RefGT:59 seqCodedCnt:    " + seqCodedCnt);
            System.out.println("RefGT:60 seqCodedCnt0:   " + seqCodedCnt0);
            System.out.println("RefGT:61 meanLength0:      " + (sumLength0 / seqCodedCnt0));
            System.out.println("RefGT:62 meanLength1:      " + (sumLength1 / seqCodedCnt));

            System.out.println();
            System.out.println("RefGT:71 overhead/fields KiB: " + (otherBytes >> 10));
            System.out.println("RefGT:72 alleleCoded     KiB: " + (alleleCodedBytes >> 10));
            System.out.println("RefGT:73 seqCoded0       KiB: " + (seqCodedBytes0 >> 10));
            System.out.println("RefGT:74 seqCoded1       KiB: " + (seqCodedBytes1 >> 10));
            System.out.println("RefGT:75 TOTAL           KiB: " + (estBytes       >> 10));
        }
        return estBytes;
    }

    /**
     * Constructs a new {@code RefGT} instance.
     * @param markers the sequence of markers
     * @param samples the sequence of samples
     * @param refVcfRecs the sequence of per-marker genotype data
     *
     * @throws IllegalArgumentException if
     * {@code markers.nMarkers() != refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].samples().equals(samples) == false} for any
     * {@code k} satisfying {@code 0 <= k  && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].marker().equals(markers.marker(k)) == false}
     * for any {@code k} satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].isPhased() == false} for any {@code k}
     * satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws NullPointerException if
     * {@code markers == null || samples == null || refVcfRecs == null
     * || refVcfRecs[k] == null} for any {@code k} satisfying
     * {@code 0 <= k && k <= refVcfRecs.length}
     */
    public RefGT(Markers markers, Samples samples, RefGTRec[] refVcfRecs) {
        checkData(markers, samples, refVcfRecs);
        this.markers = markers;
        this.samples = samples;
        this.recs = refVcfRecs.clone();
    }

    /**
     * Constructs a new {@code RefHapPairs} instance.
     * @param refVcfRecs the sequence of per-marker genotype data
     *
     * @throws IllegalArgumentException if {@code refVcfRecs.length == 0}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].samples().equals(samples) == false} for any
     * {@code k} satisfying {@code 0 <= k  && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].isPhased() == false} for any {@code k}
     * satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws NullPointerException if
     * {@code samples == null || refVcfRecs == null}
     * @throws NullPointerException if
     * {@code (refVcfRecs[k] == null)} for any {@code k} satisfying
     * {@code (0 <= k && k <= refVcfRecs.length)}
     */
    public RefGT(RefGTRec[] refVcfRecs) {
        this.samples = checkData(refVcfRecs);
        Marker[] ma = Arrays.stream(refVcfRecs)
                .parallel()
                .map(rec -> rec.marker())
                .toArray(Marker[]::new);
        this.markers = Markers.create(ma);
        this.recs = refVcfRecs.clone();
    }

    private static Samples checkData(GTRec[] refVcfRecs) {
        if (refVcfRecs.length==0) {
            String s = "Missing data in VCF file";
            throw new IllegalArgumentException(s);
        }
        Samples samples = refVcfRecs[0].samples();
        for (int j=0; j<refVcfRecs.length; ++j) {
            if (refVcfRecs[j].samples().equals(samples)==false) {
                String s = "sample inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].isPhased()==false) {
                String s = "non-reference data at marker index " + j;
                throw new IllegalArgumentException(s);
            }
        }
        return samples;
    }

    private static void checkData(Markers markers, Samples samples,
            GTRec[] refVcfRecs) {
        if (markers.size()!=refVcfRecs.length) {
            String s = "markers.nMarkers()=" + markers.size()
                    + " refVcfRecs.length=" + refVcfRecs.length;
            throw new IllegalArgumentException(s);
        }
        for (int j=0; j<refVcfRecs.length; ++j) {
            if (refVcfRecs[j].samples().equals(samples)==false) {
                String s = "sample inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].marker().equals(markers.marker(j))==false) {
                String s = "marker inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].isPhased()==false) {
                String s = "non-reference data at marker index " + j;
                throw new IllegalArgumentException(s);
            }
        }
    }

    @Override
    public boolean isReversed() {
        return false;
    }

    @Override
    public int nMarkers() {
       return markers.size();
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int nHaps() {
        return 2*samples.size();
    }

    @Override
    public int nSamples() {
        return samples.size();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public boolean isPhased() {
        return true;
    }

    @Override
    public int allele(int marker, int haplotype) {
        return recs[marker].get(haplotype);
    }

    /**
     * Returns a {@code RefGT} instance restricted to genotype data for
     * the specified markers.
     * @param refGT the {@code RefGT} instance to be restricted
     * @param indices a list of distinct marker indices (from
     * {@code this.markers())} in increasing order
     * @return a {@code RefGT} instance restricted to genotype data for
     * the specified markers
     *
     * @throws IndexOutOfBoundsException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (indices[j] < 0 || indices[j] >= gt.nMarkers())}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (1 <= j && j < indices.length)} such that
     * {@code (indices[j] <= indice[j - 1])}
     * @throws NullPointerException if
     * {@code gt == null || indices == null}
     */
    public static RefGT restrict(RefGT refGT, int[] indices) {
        RefGTRec[] rra = new RefGTRec[indices.length];
        for (int j=0; j<rra.length; ++j) {
            if (j>0 && indices[j] <= indices[j-1]) {
                throw new IllegalArgumentException(String.valueOf(indices[j]));
            }
            rra[j] = refGT.recs[indices[j]];
        }
        return new RefGT(refGT.markers, refGT.samples, rra);
    }

    @Override
    public RefGT restrict(Markers markers, int[] indices) {
        RefGTRec[] rra = new RefGTRec[indices.length];
        for (int j=0; j<rra.length; ++j) {
            if (j>0 && indices[j] <= indices[j-1]) {
                throw new IllegalArgumentException(String.valueOf(indices[j]));
            }
            rra[j] = this.recs[indices[j]];
        }
        return new RefGT(markers, samples, rra);
    }

    @Override
    public RefGT restrict(int start, int end) {
        Markers restrictMarkers = markers.restrict(start, end);
        RefGTRec[] restrictRecs = IntStream.range(start, end)
                .mapToObj(j -> recs[j])
                .toArray(RefGTRec[]::new);
        return new RefGT(restrictMarkers, samples, restrictRecs);
    }

    /**
     * Returns the {@code RefGTRec} for the specified marker.
     * @param marker the marker index
     * @return the {@code RefGTRec} for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public RefGTRec get(int marker) {
        return recs[marker];
    }
}

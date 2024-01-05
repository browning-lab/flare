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

import beagleutil.ChromInterval;
import java.io.Closeable;
import java.util.Optional;
import vcf.GeneticMap;
import vcf.MarkerMap;
import vcf.Markers;
import vcf.RefGT;
import vcf.Samples;

/**
 * <p>Class {@code AdmixChromData} represents the immutable input data for
 * a local ancestry inference on a chromosome.</p>
 *
 * <p>Class {@code AdmixChromData} is immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixChromData {

    private final AdmixPar par;
    private final RefGT targRefGT;
    private final MarkerMap map;
    private final int nRefHaps;
    private final int nTargHaps;

    private AdmixChromData(AdmixPar par, RefGT targRefGT, int nRefHaps,
            MarkerMap map) {
        this.par = par;
        this.targRefGT = targRefGT;
        this.map = map;
        this.nRefHaps = nRefHaps;
        this.nTargHaps = targRefGT.nHaps() - nRefHaps;
    }

    /**
     * Returns the command line parameters
     * @return the command line parameters
     */
    public AdmixPar par() {
        return par;
    }

    /**
     * Returns the phased reference and target genotypes with
     * target samples preceding reference samples.
     * @return the phased reference and target genotypes with
     * target samples preceding reference samples
     */
    public RefGT targRefGT() {
        return targRefGT;
    }

    /**
     * Returns the genetic position of each marker.
     * @return the genetic position of each marker
     */
    public MarkerMap map() {
        return map;
    }

    /**
     * Returns the number of reference haplotypes.
     * @return the number of reference haplotypes
     */
    public int nRefHaps() {
        return nRefHaps;
    }

    /**
     * Returns the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return nTargHaps;
    }

    /**
     * Returns the total number of reference and target haplotypes.
     * @return the total number of reference and target haplotypes
     */
    public int nHaps() {
        return targRefGT.nHaps();
    }

    /**
     * Class {@code AdmixChromData.It} constructs and returns
     * {@code AdmixChromData} objects.
     */
    public static class It implements Closeable {

        private final AdmixPar par;
        private final AdmixReader reader;
        private final GeneticMap genMap;

        /**
         * Constructs an {@code AdmixChromData.It} instance from the
         * specified data.
         * @param par the command line parameters
         * @throws NullPointerException if {@code par == null}
         */
        public It(AdmixPar par) {
            ChromInterval chromInt = null;
            this.par = par;
            this.reader = AdmixReader.instance(par);
            this.genMap = GeneticMap.geneticMap(par.map(), chromInt);
        }

        /**
         * Reads and returns the immutable input data for the next chromosome.
         * Returns {@code Optional.empty()} if {@code this.close()} has
         * previously been invoked or if there is no more phased genotype data.
         * @return the input data for the next chromosome
         * @throws IllegalArgumentException if any chromosome is missing from
         * the genetic map specified on the command line
         */
        public Optional<AdmixChromData> nextChrom() {
            Optional<RefGT> optGT = reader.nextChrom();
            int nRefHaps = reader.refSamples().size()<<1;
            if (optGT.isPresent()) {
                RefGT targRefGT = optGT.get();
                Markers markers = targRefGT.markers();
                double meanGenDiff = MarkerMap.meanSingleBaseGenDist(genMap, markers);
                MarkerMap map = MarkerMap.create(genMap, meanGenDiff, targRefGT.markers());
                AdmixChromData ad = new AdmixChromData(par, targRefGT, nRefHaps, map);
                return Optional.of(ad);
            }
            else {
                return Optional.empty();
            }
        }

        /**
         * Returns the list of reference samples.
         * @return the list of reference samples
         */
        public Samples refSamples() {
            return reader.refSamples();
        }

        /**
         * Returns the list of target samples.
         * @return the list of target samples
         */
        public Samples targSamples() {
            return reader.targSamples();
        }

        /**
         * Returns the list of input reference and target samples.
         * The reference samples precede the target samples.
         * @return the list of input reference and target samples
         */
        public Samples allSamples() {
            return reader.allSamples();
        }

        /**
         * Returns the cumulative number of markers returned by previous
         * invocations of {@code this.nextChrom()}.
         * @return the cumulative number of markers
         */
        public int nMarkers() {
            return reader.nMarkers();
        }

        /**
         * Releases any I/O resources held by this object.
         */
        @Override
        public void close() {
            reader.close();
        }
    }
}

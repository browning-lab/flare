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

/**
 * <p>Class {@code AdmixData} represents the input data for a chromosome
 * and analysis parameters for local ancestry inference.</p>
 *
 * <p>Instances of class {@code AdmixData} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixData {

    private final AdmixChromData chromData;
    private final ParamsInterface params;
    private final AdmixHmmProbs hmmProbs;

    /**
     * Constructs a new {@code AdmixData} for the specified data.
     * @param chromData the immutable input data for a chromosome
     * @param params the analysis parameters
     * @throws IllegalArgumentException if
     * {@code chromData.par() != params.sampleData().par()}
     * @throws NullPointerException if
     * {@code chromData == null || params == null}
     */
    public AdmixData(AdmixChromData chromData, ParamsInterface params) {
        if (chromData.par()!=params.sampleData().par()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        this.chromData = chromData;
        this.params = params;
        this.hmmProbs = new AdmixHmmProbs(params, chromData.map());
    }

    /**
     * Returns the immutable input data for the chromosome.
     * @return the immutable input data for the chromosome
     */
    public AdmixChromData chromData() {
        return chromData;
    }

    /**
     * Returns the analysis parameters.
     * @return the analysis parameters
     */
    public ParamsInterface params() {
        return params;
    }

    /**
     * Returns HMM emission and the components of HMM transition probabilities.
     * @return HMM emission and the components of HMM transition probabilities
     */
    public AdmixHmmProbs hmmProbs() {
        return hmmProbs;
    }
}

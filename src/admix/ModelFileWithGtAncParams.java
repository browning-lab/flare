/*
 * Copyright 2023 Genomics plc
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
 * <p>Class {@code ModelFileWithGtAncParams} represents the analysis
 * parameters for a local ancestry inference analysis when the `gt-ancestries`
 * parameter is used.</p>
 *
 * @author Thomas Hickman {@code <thomas.hickman@genomicsplc.com>}
 */
public class ModelFileWithGtAncParams extends ModelFileParams {
    private final GtAncParser gtAncestries;
    private final Integer sampleIndex;

    public ModelFileWithGtAncParams(FixedParams fixedParams, GtAncParser gtAncestries, int sampleI) {
        super(fixedParams);
        this.gtAncestries = gtAncestries;
        this.sampleIndex = sampleI;
    }

    @Override
    public double[] mu() {
        return this.gtAncestries.getMu(sampleIndex).clone();
    }
}

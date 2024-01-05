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
package bref;

import java.io.Closeable;
import vcf.RefGTRec;
import vcf.Samples;

/**
 * <p>Interface {@code BrefWrites} writes phased, non-missing genotypes to a
 * binary reference format (bref) file.  The {@code close()} method must
 * be called after the last invocation of the {@code write()} method
 * in order to ensure that any buffered data are written to the output
 * binary reference file.
 * </p>
 * <p>Instances of class {@code BrefWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface BrefWriter extends Closeable {

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();

    /**
     * Writes the specified phased genotype data in binary reference format.
     * The Java virtual machine will exit with an error message if an I/O
     * error occurs during method execution, if {@code this.close()}
     * has previously been invoked, or if
     * {@code rec.samples().equals(this.samples()) == false}.
     *
     * @param rec phased genotype data
     *
     * @throws NullPointerException if {@code rec == null}
     */
    void write(RefGTRec rec);

    /**
     * Flushes any buffered output and releases any system resources that are
     * held by this {@code BrefWriter}.  The Java virtual machine will exit
     * with an error message if an I/O error occurs during method execution.
     */
    @Override
    void close();
}

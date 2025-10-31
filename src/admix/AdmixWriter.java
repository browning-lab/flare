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

import blbutil.BGZIPOutputStream;
import blbutil.FileUtil;
import blbutil.Utilities;
import ints.UnsignedByteArray;
import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import vcf.Samples;
import vcf.VcfHeader;
import vcf.VcfWriter;

/**
 * <p>Class {@code AdmixWriter} contains methods for writing
 * inferred local ancestry to an output VCF file and for writing
 * per-global global ancestry probabilities.</p>
 *
 * <p>Instances of class {@code AdmixWriter} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AdmixWriter implements Closeable {

    private final AdmixPar par;
    private final int nTargHaps;
    private final int nAnc;
    private final File outFile;
    private final OutputStream os;

    /**
     * Constructs a new {@code AdmixWriter} instance from the specified data.
     * The Java virtual machine will exit with an error message if there
     * is an error writing the VCF header lines.
     * @param sampleData the fixed parameters for a local ancestry inference
     * analysis
     * @param targVcfHeader the target VCF meta-information lines and header
     * line
     * @throws NullPointerException if
     * {@code ((sampleData == null) || (targVcfHeader == null))}
     */
    public AdmixWriter(SampleData sampleData, VcfHeader targVcfHeader) {
        Samples targSamples = sampleData.targSamples();
        this.par = sampleData.par();
        this.outFile = new File(par.out() + ".anc.vcf.gz");
        this.nTargHaps = targSamples.size() << 1;
        this.nAnc = sampleData.nAnc();
        this.os = FileUtil.bufferedOutputStream(outFile, false);
        printVcfHeader(sampleData, targVcfHeader, outFile, os);
    }

    private static void printVcfHeader(SampleData sampleData,
            VcfHeader targVcfHeader, File outFile, OutputStream os) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter out = new PrintWriter(
                new BGZIPOutputStream(baos, false))) {
            printVcfHeader(sampleData, targVcfHeader, out);
        }
        UnsignedByteArray byteArray = new UnsignedByteArray(baos);
        print(byteArray, os, outFile);
    }

    private static void printVcfHeader(SampleData sampleData,
            VcfHeader targVcfHeader, PrintWriter out) {
        int nAnc = sampleData.nAnc();
        String[] ancIds = sampleData.ancIds();
        Samples targSamples = sampleData.targSamples();
        boolean printProbs = sampleData.par().probs();
        VcfWriter.copyMetaInfoLines(targVcfHeader, out);
        out.println(VcfHeader.metaInfoLine("flareCommand",
                AdmixPar.flareCommand(), true));
        out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        out.println("##FORMAT=<ID=AN1,Number=1,Type=Integer,Description=\"Ancestry of first haplotype\">");
        out.println("##FORMAT=<ID=AN2,Number=1,Type=Integer,Description=\"Ancestry of second haplotype\">");
        if (printProbs) {
            out.print("##FORMAT=<ID=ANP1,Number=");
            out.print(nAnc);
            out.println(",Type=Float,Description=\"Ancestry probabilities for first haplotype\">");
            out.print("##FORMAT=<ID=ANP2,Number=");
            out.print(nAnc);
            out.println(",Type=Float,Description=\"Ancestry probabilities for second haplotype\">");
        }
        for (int j=0; j<nAnc; ++j) {
            out.print(j==0 ? "##ANCESTRY=<" : ",");
            out.print(ancIds[j]);
            out.print("=");
            out.print(j);
        }
        out.println(">");
        VcfWriter.writeHeaderLine(targSamples.ids(), out);
    }

    /**
     * Returns the command line parameters.
     * @return the command line parameters
     */
    public AdmixPar par() {
        return par;
    }

    /**
     * Returns the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return nTargHaps;
    }

    /**
     * Returns the number of ancestries.
     * @return the number of ancestries
     */
    public int nAnc() {
        return nAnc;
    }

    /**
     * Write VCF Records with estimated local ancestry.
     * The Java virtual machine will exit with an error message if there
     * is an error writing the output VCF records.
     * @param estAnc the estimated local ancestry
     * @throws IllegalArgumentException if
     * {@code (estAnc.nTargHaps() != this.nTargHaps()}
     * @throws IllegalArgumentException if
     * {@code (estAnc.nAnc() != this.nAnc()}
     * @throws NullPointerException {@code estAnc == null)}
     */
    public void writeAncestry(EstimatedAncestry estAnc) {
        if (estAnc.nTargHaps()!=nTargHaps) {
            throw new IllegalArgumentException(String.valueOf(estAnc.nTargHaps()));
        }
        if (estAnc.nAnc() != nAnc) {
            throw new IllegalArgumentException(String.valueOf(estAnc.nAnc()));
        }
        int nMarkers = estAnc.nMarkers();
        int nMarkersPerBatch = nMarkersPerBatch();
        for (int start=0; start<nMarkers; start+=nMarkersPerBatch) {
            int end = Math.min(start + nMarkersPerBatch, nMarkers);
            UnsignedByteArray[] byteArrays = estAnc.writeAncestry(start, end);
            for (UnsignedByteArray byteArray : byteArrays) {
                print(byteArray, os, outFile);
            }
        }
    }

    private int nMarkersPerBatch() {
        int BUFFER_MB_PER_THREAD = 100;
        long totalMb = BUFFER_MB_PER_THREAD*par.nthreads();
        int minOutBytesPerHap = par.probs() ? (8 + 4*nAnc) : 8;
        long minOutBytesPerMarker = nTargHaps*minOutBytesPerHap;
        return (int) ( (totalMb *(1<<20)) / minOutBytesPerMarker );
    }

    private static void print(UnsignedByteArray outBytes, OutputStream os,
            File outFile)  {
        try {
            outBytes.write(os);
        } catch (IOException ex) {
            fileOutputError(ex, outFile);
        }
    }

    private static void fileOutputError(Exception ex, File file) {
        Utilities.exit(ex, "I/O error for file: " + file);
    }

    /**
     * Closes this {@code AdmixWriter} and releases system resources held
     * by it. If this {@code AdmixWriter} is already closed, then invoking
     * this method has no effect.
     */
    @Override
    public void close() {
        try {
            BGZIPOutputStream.writeEmptyBlock(os);
            os.close();
        }
        catch (IOException ex) {
            fileOutputError(ex, outFile);
        }
    }
}

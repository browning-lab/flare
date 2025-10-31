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
import blbutil.Const;
import blbutil.FloatArray;
import ints.IntArray;
import ints.UnsignedByteArray;
import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.stream.IntStream;
import vcf.MarkerUtils;
import vcf.RefGT;

/**
 * <p>Class {@code AdmixRecBuilder} contains methods for constructing
 * and printing VCF records in VCF 4.3 format.</p>
 *
 * <p>Instances of class {@code AdmixRecBuilder} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AdmixRecBuilder {

    private final SampleData sampleData;
    private final RefGT targRefGT;
    private final int mStart;
    private final int mEnd;
    private final int nTargHaps;
    private final int nAnc;
    private final int probStart;
    private final boolean printProbs;
    private final StringBuilder[] recData;

    private int hapCnt;

    /**
     * Constructs a new {@code AdmixdRecBuilder} instance for the specified
     * data.
     * @param sampleData reference and target sample metadata
     * @param targRefGT the phased reference and target haplotypes with
     * reference samples preceding target samples
     * @param start the first marker index (inclusive)
     * @param end the last marker index (exclusive)
     * @throws IllegalArgumentException if
     * {@code (sampleData.refSamples().size() + sampleData.targSamples().size())<<1 != targRefGT.nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || end > targRefGT.nMarkers() || end < start)}
     * @throws NullPointerException if {@code sampleData == null || targRefGT == null}
     */
    public AdmixRecBuilder(SampleData sampleData, RefGT targRefGT, int start,
            int end) {
        if (start<0) {
            throw new IndexOutOfBoundsException(String.valueOf(start));
        }
        if (end<start || end>targRefGT.nMarkers()) {
            throw new IndexOutOfBoundsException(String.valueOf(end));
        }
        int nHaps = (sampleData.refSamples().size() + sampleData.targSamples().size())<<1;
        if (targRefGT.nHaps() != nHaps) {
            throw new IllegalArgumentException(String.valueOf(targRefGT.nHaps()));
        }
        int mSize = end - start;
        this.sampleData = sampleData;
        this.targRefGT = targRefGT;
        this.mStart = start;
        this.mEnd = end;
        this.nTargHaps = (sampleData.targSamples().size() << 1);
        this.nAnc = sampleData.nAnc();
        this.printProbs = sampleData.par().probs();
        this.probStart = nAnc*start;
        this.recData = new StringBuilder[mSize];
        int charPerHap = printProbs ? nAnc<<3 : nAnc<<1;
        for (int j=0; j<mSize; ++j) {
            recData[j] = new StringBuilder(80 + charPerHap*nTargHaps);
        }
        this.hapCnt = 0;
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
     * Returns the fixed parameters for a local ancestry inference
     * analysis
     * @return the fixed parameters for a local ancestry inference
     * analysis
     */
    public SampleData sampleData() {
        return sampleData;
    }

    /**
     * Returns the start marker index (inclusive).
     * @return the start marker index (inclusive)
     */
    public int start() {
        return mStart;
    }

    /**
     * Returns the end marker index (exclusive).
     * @return the end marker index (exclusive)
     */
    public int end() {
        return mEnd;
    }

    /**
     * Returns the number of haplotypes with ancestry probabilities
     * that have been added with the {@code this.addSampleData()} method.
     * @return the number of haplotypes with ancestry probabilities
     * that have been added with the {@code this.addSampleData()} method
     */
    public int hapCnt() {
        return hapCnt;
    }

    /**
     * Adds the specified sample data to the VCF records for markers with index
     * between {@code this.start()} (inclusive) and {@code this.end()}
     * (exclusive). The contract for this method is undefined
     * if the ancestry probabilities are not added in the order that target
     * samples appear in this.targRefGT(), or if
     * {@code this.par().probs() == true} and any ancestry probability is less
     * than 0.0 or greater than 1.0 for any marker between {@code this.start()}
     * (inclusive) and {@code this.end()} (exclusive).
     * @param probs1 the ancestry probabilities for the first haplotype
     * @param probs2 the ancestry probabilities for the second haplotype
     * @param probFormatter object for formatting probabilities
     * @throws IndexOutOfBoundsException if
     * {@code probs1.size() < (this.end() * this.par().nAnc())}
     * @throws IndexOutOfBoundsException if
     * {@code probs2.size() < (this.end() * this.par().nAnc())}
     * @throws NullPointerException if {@code probs1 == null || probs2 == null}
     */
    public void addSampleData(FloatArray probs1, FloatArray probs2,
            ProbFormatter probFormatter) {
        int hap1 = hapCnt++;
        int hap2 = hapCnt++;
        for (int m=mStart, pStart=probStart; m<mEnd; ++m, pStart+=nAnc) {
            int pEnd = pStart + nAnc;
            StringBuilder sb = recData[m-mStart];
            sb.append(Const.tab);
            sb.append(targRefGT.allele(m, hap1));
            sb.append(Const.phasedSep);
            sb.append(targRefGT.allele(m, hap2));
            sb.append(Const.colon);
            sb.append(maxOffset(probs1, pStart, pEnd));
            sb.append(Const.colon);
            sb.append(maxOffset(probs2, pStart, pEnd));
            if (printProbs) {
                for (int j=pStart; j<pEnd; ++j) {
                    sb.append(j==pStart ? Const.colon : Const.comma);
                    sb.append(probFormatter.format(probs1.get(j)));
                }
                for (int j=pStart; j<pEnd; ++j) {
                    sb.append(j==pStart ? Const.colon : Const.comma);
                    sb.append(probFormatter.format(probs2.get(j)));
                }
            }
        }
    }

    /**
     * Adds the specified sample data to the VCF records for markers with index
     * between {@code this.start()} (inclusive) and {@code this.end()}
     * (exclusive). The contract for this method is undefined
     * if the ancestry probabilities are not added in the order that target
     * samples appear in {@code this.targRefGT()}.
     * @param anc1 the ancestry of the first haplotype
     * @param anc2 the ancestry of the second haplotype
     * @throws IndexOutOfBoundsException if
     * {@code (anc.size() < this.end() || anc2.size() < this.end())}
     * @throws NullPointerException if {@code (anc1 == null || anc2 == null)}
     */
    public void addSampleData(IntArray anc1, IntArray anc2) {
        int hap1 = hapCnt++;
        int hap2 = hapCnt++;
        for (int m=mStart; m<mEnd; ++m) {
            StringBuilder sb = recData[m-mStart];
            sb.append(Const.tab);
            sb.append(targRefGT.allele(m, hap1));
            sb.append(Const.phasedSep);
            sb.append(targRefGT.allele(m, hap2));
            sb.append(Const.colon);
            sb.append(anc1.get(m));
            sb.append(Const.colon);
            sb.append(anc2.get(m));
        }
    }

    private static int maxOffset(FloatArray ancProbs, int start, int end) {
        int maxIndex = start;
        for (int j=start; j<end; ++j) {
            if (ancProbs.get(j)>ancProbs.get(maxIndex)) {
                maxIndex = j;
            }
        }
        return maxIndex - start;
    }

    /**
     * GZIP-compresses the VCF records and returns the compressed
     * records as an {@code UnsignedByteArray}.
     * @return an unsigned byte array containing compressed VCF records
     * @throws IllegalStateException if
     * {@code this.hapCnt() != (this.sampleData().targSamples().size() << 1)}
     * when this method is invoked
     */
    public UnsignedByteArray toUnsignedByteArray() {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter out = new PrintWriter(
                new BGZIPOutputStream(baos, false))) {
            printRecs(out);
        }
        return new UnsignedByteArray(baos);
    }

    /**
     * Prints the VCF records to the specified {@code PrintWriter}.
      *@param out the {@code PrintWriter} to which the VCF record will be
     * printed
     * @throws IllegalStateException if
     * {@code this.hapCnt() != (this.sampleData().targSamples().size() << 1)}
     * when this method is invoked
     * @throws NullPointerException if {@code out == null}
     */
    private void printRecs(PrintWriter out) {
        if (hapCnt != nTargHaps) {
            throw new IllegalStateException(String.valueOf(hapCnt));
        }
        for (int m=mStart; m<mEnd; ++m) {
            MarkerUtils.printFirst8Fields(targRefGT.marker(m), out);
            out.print("\tGT:AN1:AN2");           // FORMAT
            if (printProbs) {
                out.print(":ANP1:ANP2");
            }
            out.println(recData[m-mStart]);
        }
    }

    /**
     * <p>Class {@code ProbFormatter} converts probabilities to strings.</p>
     *
     * <p>Instances of class {@code ProbFormatter} are immutable.</p>
     *
     */
    public static class ProbFormatter {

        private final String[] probs;
        private final int factor;

        /**
         * Constructs a new {@code ProbFormatter} instance.
         * @param decimalPlaces the number of decimal places after rounding
         * @throws IllegalArgumentException if
         * {@code (decimalPlaces != 2) || (decimalPlaces != 3)}
         */
        public ProbFormatter(int decimalPlaces) {
            switch (decimalPlaces) {
                case 2:
                    {
                        DecimalFormat df = new DecimalFormat("#.##");
                        this.probs = IntStream.range(0, 101)
                                .parallel()
                                .mapToObj(j -> df.format(0.01*j))
                                .toArray(String[]::new);
                        this.factor = 100;
                        break;
                    }
                case 3:
                    {
                        DecimalFormat df = new DecimalFormat("#.###");
                        this.probs = IntStream.range(0, 1001)
                                .parallel()
                                .mapToObj(j -> df.format(0.001*j))
                                .toArray(String[]::new);
                        this.factor = 1000;
                        break;
                    }
                default:
                    throw new IllegalArgumentException(String.valueOf(decimalPlaces));
            }
        }

        /**
         * Returns a string representation of the specified probability.
         * The format of the string representation is unspecified and subject
         * to change.
         * @param prob a probability
         * @return a string representation of the specified probability
         * @throws IllegalArgumentException if {@code prob < 0f || prob > 1f}
         */
        public String format(float prob) {
            if (prob<0f || prob>1f) {
                throw new IllegalArgumentException(String.valueOf(prob));
            }
            return probs[(int) Math.rint(factor*prob)];
        }
    }
}

/*
 * Copyright 2021 Brian L. Browning
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
import blbutil.Filter;
import blbutil.StringUtil;
import blbutil.Utilities;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code VcfHeader} represents the Variant Call Format (VCF)
 * meta-information lines and the Variant Call Format header line
 * that precede the first Variant Call Format record.
 * </p>
 * <p>Instances of class {@code VcfHeader} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfHeader  {

    /**
     * A string equal to the first nine tab-delimited fields of a VCF header
     * line that contains sample data.
     */
    public static final String HEADER_PREFIX =
            "#CHROM" + Const.tab + "POS"
            + Const.tab + "ID" + Const.tab + "REF" + Const.tab + "ALT"
            + Const.tab + "QUAL" + Const.tab + "FILTER" + Const.tab + "INFO"
            + Const.tab + "FORMAT";

    private static final int SAMPLE_OFFSET = 9;

    private final String src;
    private final VcfMetaInfo[] metaInfoLines;
    private final int nHeaderFields;
    private final int[] includedIndices;
    private final Samples samples;

    /**
     * Returns a boolean array whose {@code k}-th value is {@code true}
     * if the FORMAT field for the {@code k}-th sample in a VCF record
     * contains an allele separator character and returns {@code false}
     * otherwise.  The contract for this method is undefined of the
     * specified string is not a properly-formatted VCF record.
     * @param vcfRec a VCF record
     * @return  a boolean array whose {@code k}-th value is {@code true}
     * if the FORMAT field for the {@code k}-th sample does not contain an
     * allele separator character

     * @throws NullPointerException if {@code vcfHeader == null || rec == null}
     */
    public static boolean[] isDiploid(String vcfRec) {
        List<Boolean> list = new ArrayList<>();
        int start = VcfRecGTParser.ninthTabPos(vcfRec) + 1;
        boolean noAlleleSep = true;
        for (int j=start, n=vcfRec.length(); j<n; ++j)  {
            char c = vcfRec.charAt(j);
            if (c==Const.tab) {
                list.add(noAlleleSep==false);
                noAlleleSep = true;
            }
            else if (c==Const.unphasedSep || c==Const.phasedSep) {
                noAlleleSep = false;
            }
        }
        list.add(noAlleleSep==false);
        boolean[] isDiploid = new boolean[list.size()];
        for (int j=0; j<isDiploid.length; ++j) {
            isDiploid[j] = list.get(j);
        }
        return isDiploid;
    }

    /**
     * Returns a VCF header object for the specified VCF meta information lines
     * and header line. The header line must be the last line in the
     * specified {@code lines} array.
     * @param src a string describing the source of the VCF file
     * @param lines the VCF meta-information and header lines
     * @param isDiploid a boolean array whose {@code k}-th value is {@code true}
     * if the {@code k}-th sample is diploid, and is {@code false} if the
     * {@code k}-th sample is haploid
     * @throws IllegalArgumentException if a format error is encountered
     * in a meta-information line or header lines}
     * @throws NullPointerException if
     * {@code src==null || lines == null || isDiploid == nulle}
     */
    public VcfHeader(String src, String[] lines, boolean[] isDiploid) {
        this(src, lines, isDiploid, Filter.acceptAllFilter());
    }

    /**
     * Returns a VCF header object for the specified VCF meta information lines
     * and header line. The header line must be the last line in the
     * specified {@code lines} array.
     * @param src a string describing the source of the VCF file
     * @param lines the VCF meta-information and header lines
     * @param isDiploid a boolean array whose {@code k}-th value is {@code true}
     * if the {@code k}-th sample is diploid, and is {@code false} if the
     * {@code k}-th sample is haploid
     * @param sampleFilter a sample filter
     * @throws IllegalArgumentException if a format error is encountered
     * in a meta-information line or header lines}
     * @throws NullPointerException if
     * {@code src==null || lines == null || isDiploid == null}
     */
    public VcfHeader(String src, String[] lines, boolean[] isDiploid,
            Filter<String> sampleFilter) {
        if (src==null) {
            throw new NullPointerException(String.class.toString());
        }
        if (sampleFilter==null) {
            sampleFilter = Filter.acceptAllFilter();
        }
        checkHeaderLines(lines, src);
        int headerIndex = lines.length-1;
        this.src = src;
        this.metaInfoLines = new VcfMetaInfo[headerIndex];
        for (int j=0; j<headerIndex; ++j) {
            this.metaInfoLines[j] = new VcfMetaInfo(lines[j]);
        }
        String[] headerFields = StringUtil.getFields(lines[headerIndex], Const.tab);
        this.nHeaderFields = headerFields.length;
        this.includedIndices = includedIndices(src, headerFields, sampleFilter);
        this.samples = samples(headerFields, isDiploid, includedIndices);
    }

    private static void checkHeaderLines(String[] lines, String src) {
        if (lines.length==0) {
            String s = Const.nl + Const.nl
                    + "ERROR: Missing the VCF meta information lines and the VCF header line"
                    + Const.nl + "VCF source: " + src
                    + Const.nl;
            throw new IllegalArgumentException(s);
        }
        String line = lines[lines.length-1];
        if (line.startsWith(HEADER_PREFIX) == false) {
                String s = Const.nl + Const.nl
                        + "ERROR: Missing the VCF header line."
                        + Const.nl + "VCF source: " + src
                        + Const.nl + "The VCF header line must immediately follow the meta-information lines."
                        + Const.nl + "The fields of the VCF header line must be tab-delimited and begin with:"
                        + Const.nl + HEADER_PREFIX
                        + Const.nl;
                throw new IllegalArgumentException(s);
        }
    }

    private static int[] includedIndices(String src, String[] headerFields,
            Filter<String> sampleFilter) {
        int nUnfilteredSamples = Math.max(headerFields.length - SAMPLE_OFFSET, 0);
        int[] includedIndices = new int[nUnfilteredSamples];
        int index = 0;
        for (int j=0; j<nUnfilteredSamples; ++j) {
            if (sampleFilter.accept(headerFields[SAMPLE_OFFSET + j])) {
                includedIndices[index++] = j;
            }
        }
        if (index==0) {
            String err = "All samples in the VCF file are excluded";
            String info = Const.nl + "Error      :  " + err
                    + Const.nl     + "File       :  " + src;
            Utilities.exit(new Throwable(err), info);
        }
        if (index < includedIndices.length) {
            includedIndices = Arrays.copyOf(includedIndices, index);
        }
        return includedIndices;
    }

    private Samples samples(String[] headerFields, boolean[] isDiploid,
            int[] includedIndices) {
        String[] ids = new String[includedIndices.length];
        boolean[] restrictedIsDiploid = new boolean[includedIndices.length];
        for (int j=0; j<ids.length; ++j) {
            ids[j] = headerFields[SAMPLE_OFFSET + includedIndices[j]];
            restrictedIsDiploid[j] = isDiploid[includedIndices[j]];
        }
        return Samples.fromIds(ids, restrictedIsDiploid);
    }

    /**
     * Returns the source from which data are read.  The string representation
     * of the source is undefined and subject to change.
     * @return the source from which data are read
     */
    public String src() {
        return src;
    }

    /**
     * Returns the number of VCF meta-information lines. VCF meta-information
     * lines are lines that precede the VCF header line. A VCF meta-information
     * line must begin with "##".
     *
     * @return the number of VCF meta-information lines
     */
     public int nMetaInfoLines() {
         return metaInfoLines.length;
     }

    /**
      * Returns the specified VCF meta-information line.

      * @param index a VCF meta-information line index
      * @return the specified VCF meta-information line
      *
      * @throws IndexOutOfBoundsException if
      * {@code index < 0 || index >= this.nMetaInfoLines()}
      */
     public VcfMetaInfo metaInfoLine(int index) {
         return metaInfoLines[index];
     }

     /**
      * Returns the number of fields in the VCF header line before sample
      * exclusions.
      * @return the number of fields in the VCF header line before sample
      * exclusions
      */
     public int nHeaderFields() {
         return nHeaderFields;
     }

     /**
      * Returns the number of samples before sample exclusions.
      * @return the number of samples before sample exclusions
      */
     public int nUnfilteredSamples() {
         return Math.max(0, nHeaderFields - SAMPLE_OFFSET);
     }

     /**
      * Returns the index of the specified sample in the original
      * list of samples before sample exclusions.
      * @param sample a sample index
      * @return the index of the specified sample in the original
      * list of samples before sample exclusions
      * @throws IndexOutOfBoundsException if
      * {@code sample < 0 || sample >= this.size()}
      */
     public int unfilteredSampleIndex(int sample) {
         return includedIndices[sample];
     }

     /**
      * Returns the number of samples after sample exclusions.
      * @return the number of samples after sample exclusions
      */
     public int nSamples() {
         return samples.size();
     }

    /**
     * Return the list of samples after sample exclusions.
     * @return the list of samples after sample exclusions
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns {@code this.sample().ids()}.
     * @return {@code this.sample().ids()}
     */
    public String[] sampleIds() {
        return samples.ids();
    }

    /**
     * Returns the VCF meta-information lines and the VCF header line after
     * applying sample exclusions.
     * @return the VCF meta-information lines and the VCF header line after
     * applying sample exclusions.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(400);
        for (int j=0; j<metaInfoLines.length; ++j) {
            sb.append(metaInfoLines[j]);
            sb.append(Const.nl);
        }
        String[] sampleIds = samples.ids();
        sb.append(HEADER_PREFIX);
        for (String id : sampleIds) {
            sb.append(Const.tab);
            sb.append(id);
        }
        sb.append(Const.nl);
        return sb.toString();
    }
}

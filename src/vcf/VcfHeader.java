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
import blbutil.FileIt;
import blbutil.InputIt;
import blbutil.StringUtil;
import blbutil.Utilities;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

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
     * The VCF meta-information line prefix: "##"
     */
    public static final String META_INFO_PREFIX = "##";

    /**
     * The first nine tab-delimited fields of a VCF header line that contains
     * sample data
     */
    public static final String HEADER_PREFIX =
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    private static final int FIRST_SAMPLE_FIELD = 9;

    private final String source;
    private final String[] metaInfoLines;
    private final int nHeaderFields;
    private final int[] filteredSampleIndices;
    private final Samples samples;

    /**
     * Returns a VCF header object for the specified VCF meta information lines
     * and header line. The header line must be the last line in the
     * specified {@code lines} array.
     * @param src a string describing the source of the VCF file
     * @param lines the VCF meta-information and header lines
     * @param isDiploid a boolean array whose {@code k}-th value is {@code true}
     * if the {@code k}-th sample is diploid, and is {@code false} if the
     * {@code k}-th sample is haploid
     * @throws IllegalArgumentException if
     * {@code (lines[j].trim().startsWith(META_INFO_PREFIX) == false)} for
     * {@code ((0 <= j) && (j < (lines.length - 1)))}
     * @throws IllegalArgumentException if
     * {@code (lines[lines.length-1].trim().startsWith(HEADER_PREFIX) == false)}
     * @throws NullPointerException if
     * {@code source==null || lines == null || isDiploid == null}
     */
    public VcfHeader(String src, String[] lines, boolean[] isDiploid) {
        this(src, lines, isDiploid, FilterUtil.acceptAllPredicate());
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
     * @param sampleFilter a sample filter or {@code null} if there
     * are no sample exclusions
     * @throws IllegalArgumentException if a format error is encountered
     * in a meta-information line or header lines}
     * @throws NullPointerException if
     * {@code source==null || lines == null || isDiploid == null}
     */
    public VcfHeader(String src, String[] lines, boolean[] isDiploid,
            Predicate<String> sampleFilter) {
        if (src==null) {
            throw new NullPointerException(String.class.toString());
        }
        if (sampleFilter==null) {
            sampleFilter = FilterUtil.acceptAllPredicate();
        }
        checkHeaderLine(lines, src);
        int headerIndex = lines.length-1;
        this.source = src;
        this.metaInfoLines = new String[headerIndex];
        for (int j=0; j<headerIndex; ++j) {
            this.metaInfoLines[j] = lines[j].trim();
            if (metaInfoLines[j].startsWith(META_INFO_PREFIX)==false) {
                String s = "Missing initial \"" + META_INFO_PREFIX
                        + "\" in meta-information line: " + metaInfoLines[j];
                throw new IllegalArgumentException(s);
            }
        }
        String[] headerFields = StringUtil.getFields(lines[headerIndex], Const.tab);
        this.nHeaderFields = headerFields.length;
        this.filteredSampleIndices = filteredSampleIndices(src, headerFields, sampleFilter);
        this.samples = samples(headerFields, isDiploid, filteredSampleIndices);
    }

    private static void checkHeaderLine(String[] lines, String src) {
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

    private static int[] filteredSampleIndices(String src, String[] headerFields,
            Predicate<String> sampleFilter) {
        int nUnfilteredSamples = Math.max(headerFields.length - FIRST_SAMPLE_FIELD, 0);
        int[] includedIndices = new int[nUnfilteredSamples];
        int index = 0;
        for (int j=0; j<nUnfilteredSamples; ++j) {
            if (sampleFilter.test(headerFields[FIRST_SAMPLE_FIELD + j])) {
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

    private static Samples samples(String[] headerFields, boolean[] isDiploid,
            int[] includedIndices) {
        String[] ids = new String[includedIndices.length];
        boolean[] restrictedIsDiploid = new boolean[includedIndices.length];
        for (int j=0; j<ids.length; ++j) {
            ids[j] = headerFields[FIRST_SAMPLE_FIELD + includedIndices[j]];
            restrictedIsDiploid[j] = isDiploid[includedIndices[j]];
        }
        return new Samples(ids, restrictedIsDiploid);
    }

    /**
     * Returns the source from which data are read.  The string representation
     * of the source is undefined and subject to change.
     * @return the source from which data are read
     */
    public String source() {
        return source;
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
     public String metaInfoLine(int index) {
         return metaInfoLines[index];
     }

     /**
      * Returns the VCF meta-information lines.
      * @return the VCF meta-information lines
      */
     public String[] metaInfoLines() {
         return metaInfoLines.clone();
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
         return Math.max(0, nHeaderFields - FIRST_SAMPLE_FIELD);
     }

     /**
      * Returns the index of the specified sample in the original
      * list of samples before sample exclusions.
      * @param sample a sample index
      * @return the index of the specified sample in the original
      * list of samples before sample exclusions
      * @throws IndexOutOfBoundsException if
      * {@code sample < 0 || sample >= this.samples().size()}
      */
     public int filteredSampleIndex(int sample) {
         return filteredSampleIndices[sample];
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
     * Returns a VCF meta-information line with the specified key and value.
     * The returned string is not terminated with an end-of-line separator.
     * @param key the meta-information line key
     * @param value the meta-information line value
     * @param quoteValue {@code true} if the value should be immediately
     * preceded and followed with a {@code '"'} character
     * @return a VCF meta-information line with the specified key and value
     * @throws NullPointerException if {@code (key == null) || (value == null))}
     */
    public static String metaInfoLine(String key, String value, boolean quoteValue) {
        if (key==null) {
            throw new NullPointerException("key is null");
        }
        if (value==null) {
            throw new NullPointerException("value is null");
        }
        StringBuilder sb = new StringBuilder();
        sb.append("##");
        sb.append(key);
        sb.append('=');
        if (quoteValue) {
            sb.append('"');
        }
        sb.append(value);
        if (quoteValue) {
            sb.append('"');
        }
        return sb.toString();
    }


    /**
     * Returns an array of length {@code (metaInfo.length + 1)} containing
     * the elements of {@code meta-info} followed by a VCF meta-information
     * with the specified key and value.
     * @param metaInfo an array
     * @param key the meta-information line key
     * @param value the meta-information line value
     * @param quoteValue {@code true} if the value should be immediately
     * preceded and followed with a {@code '"'} character
     * @return an array containing the elements of {@code meta-info}
     * followed by a VCF meta-information with the specified key and value
     * @throws NullPointerException if
     * {@code (meta-info == null) || (key == null) || (value == null)}
     */
    public static String[] addMetaInfoLine(String[] metaInfo, String key,
            String value, boolean quoteValue) {
        String[] appendedMetaInfo = Arrays.copyOf(metaInfo, (metaInfo.length + 1));
        appendedMetaInfo[metaInfo.length] = VcfHeader.metaInfoLine(key, value, quoteValue);
        return appendedMetaInfo;
    }

    /**
     * Reads and returns the VCF meta-information lines and the VCF header
     * line in the order that they are read from the specified
     * {@code FileIt<String>}. Each line, stripped of its
     * end-of-line character, is an element of the returned list.
     * @param it an iterator that returns the lines of a VCF file
     * @return the VCF meta-information lines and the VCF header line
     * @throws NullPointerException if {@code (it == null)}
     * @throws IllegalArgumentException if the lines preceding the VCF header
     * line do not begin with "##" or if the VCF header line does not
     * begin with "#CHROM"
     */
    public static ArrayList<String> readVcfHeader(FileIt<String> it) {
        return readVcfHeader(it, null, null);
    }

    /**
     * Reads and returns the VCF meta-information lines and the VCF header
     * line in the order that they are read from the specified
     * {@code FileIt<String>}. Each line, stripped of its
     * end-of-line character, is an element of the returned list.
     * If {@code (metaInfoVcf != null)}, the meta-information lines read from
     * {@code it} will be replaced with the meta-information lines read from
     * {@code metaInfoVcf}.  If {@code (samplesVcf != null)}, the VCF header
     * line read from {@code it} will be replaced with the VCF header
     * line read from {@code samplesVcf}.
     * @param it an iterator that returns the lines of a VCF file
     * @param metaInfoVcf a VCF file with replacement meta-information lines or
     * {@code null} if the original meta-information lines should be retained
     * @param samplesVcf a VCF file with replacement sample identifiers or
     * {@code null} if the original sample identifiers should be retained
     * @return the VCF meta-information lines and the VCF header line
     * @throws NullPointerException if {@code (it == null)}
     * @throws IllegalArgumentException if the specified {@code FileIt<String>}
     * does not return a line beginning with "#CHROM", or returns a line
     * that does not begin with "##" before returning a line beginning with
     * "#CHROM"
     */
    public static ArrayList<String> readVcfHeader(FileIt<String> it,
            File metaInfoVcf, File samplesVcf) {
        String hash = "##";
        ArrayList<String> lines = new ArrayList<>();
        String line = it.hasNext() ? it.next() : null;
        while (line != null && line.startsWith(hash)) {
            lines.add(line);
            line = it.hasNext() ? it.next() : null;
        }
        if (line == null) {
            throw new IllegalArgumentException("Missing VCF header line [" + it.source() + "]");
        }
        if (line.startsWith("#CHROM") == false) {
            StringBuilder sb = new StringBuilder(128);
            sb.append("VCF meta-information lines are not followed by a VCF header line: ");
            sb.append(line.substring(0, Math.min(80, line.length())));
            sb.append(" ...");
            throw new IllegalArgumentException(sb.toString());
        }

        if (metaInfoVcf!=null) {
            String[] metaInfo = metaInfo(metaInfoVcf);
            lines.clear();
            lines.addAll(Arrays.asList(metaInfo));
        }
        if (samplesVcf!=null) {
            line = headerLine(samplesVcf);
        }
        lines.add(line);
        return lines;
    }

    /**
     * Reads and returns the VCF meta-information lines and the VCF header
     * line in the order that they are read from the specified file.
     * Each line, stripped of its end-of-line character, is an element of the
     * returned list.
     * @param vcfFile an VCF file
     * @return the VCF meta-information lines and the VCF header line
     * @throws NullPointerException if {@code (it == null)}
     * @throws IllegalArgumentException if the lines preceding the VCF header
     * line do not begin with "##" or if the VCF header line does not
     * begin with "#CHROM"
     */
    public static ArrayList<String> readVcfHeader(File vcfFile) {
        int nBufferedBlocks = 16;
        try (FileIt<String> it = InputIt.fromBGZipFile(vcfFile, nBufferedBlocks)) {
            return readVcfHeader(it);
        }
    }

    /**
     * Returns the list of sample identifiers from the specified VCF header
     * file that contains VCF meta-information lines followed by a VCF header
     * line. The VCF file may contain VCF records following the VCF header,
     * but this is not required.
     * @param vcfHeaderFile a file containing VCF header lines
     * @return the sample identifiers from the specified VCF header file
     * @throws IllegalArgumentException if there does not exist a line
     * beginning with "#CHROM" in the file, or if any line preceding
     * the line beginning with "#CHROM" does not begin with "##"
     * @throws NullPointerException if {@code (vcfHeaderFile == null)}
     */
    public static String[] samples(File vcfHeaderFile) {
        String headerLine = headerLine(vcfHeaderFile);
        String[] fields = StringUtil.getFields(headerLine, Const.tab);
        int firstSampleField = 9;
        return Arrays.copyOfRange(fields, firstSampleField, fields.length);
    }

    /**
     * Returns the header line from the specified VCF header
     * file that contains VCF meta-information lines followed by a VCF header
     * line. The VCF file may contain VCF records following the VCF header,
     * but this is not required.
     * @param vcfHeaderFile a file containing VCF header lines
     * @return the header lines from the specified VCF header file
     * @throws IllegalArgumentException if there does not exist a line
     * beginning with "#CHROM" in the file, or if any line preceding
     * the line beginning with "#CHROM" does not begin with "##"
     * @throws NullPointerException if {@code (vcfHeaderFile == null)}
     */
    public static String headerLine(File vcfHeaderFile) {
        ArrayList<String> vcfHeader = VcfHeader.readVcfHeader(vcfHeaderFile);
        return vcfHeader.remove(vcfHeader.size() - 1);
    }

    /**
     * Returns the meta-information lines from the specified VCF header
     * file that contains VCF meta-information lines followed by a VCF header
     * line. The VCF file may contain VCF records following the VCF header,
     * but this is not required.
     * @param vcfHeaderFile a file containing VCF header lines
     * @return the meta-information lines from the specified VCF header file
     * @throws IllegalArgumentException if there does not exist a line
     * beginning with "#CHROM" in the file, or if any line preceding
     * the line beginning with "#CHROM" does not begin with "##"
     * @throws NullPointerException if {@code (vcfHeaderFile == null)}
     */
    public static String[] metaInfo(File vcfHeaderFile) {
        ArrayList<String> vcfHeader = VcfHeader.readVcfHeader(vcfHeaderFile);
        String headerLine = vcfHeader.remove(vcfHeader.size() - 1);
        return vcfHeader.toArray(String[]::new);
    }

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

     * @throws NullPointerException if {@code (vcfRec == null)}
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
     * Returns a {@code String} containing the VCF meta-information lines
     * and the post-sample-filtering VCF header line. Each line in the
     * {@code String} is terminated with a line separator.
     * @return the VCF meta-information lines and the VCF header line after
     * applying sample exclusions
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

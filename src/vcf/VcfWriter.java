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

import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.stream.IntStream;

/**
 * <p>Class {@code VcfWriter} contains static methods for writing data in
 * VCF 4.2 format.
 * </p>
 * <p>Instances of class {@code VcfWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfWriter {

    /**
     * The version of the VCF specification
     */
    public static final String VCF_VERSION = "VCFv4.2";

    private static final String AF_INFO = "##INFO=<ID=AF,Number=A,Type=Float,"
            + "Description=\"Estimated ALT Allele Frequencies\">";
    private static final String DR2_INFO = "##INFO=<ID=DR2,Number=A,Type=Float,"
            + "Description=\"Dosage R-Squared: estimated squared correlation between "
            + "estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">";
    private static final String IMP_INFO = "##INFO=<ID=IMP,Number=0,Type=Flag,"
            + "Description=\"Imputed marker\">";
    private static final String DS_FORMAT = "##FORMAT=<ID=DS,Number=A,Type=Float,"
            +"Description=\"estimated ALT dose [P(RA) + 2*P(AA)]\">";
    private static final String AP1_FORMAT = "##FORMAT=<ID=AP1,Number=A,Type=Float,"
            +"Description=\"estimated ALT dose on first haplotype\">";
    private static final String AP2_FORMAT = "##FORMAT=<ID=AP2,Number=A,Type=Float,"
            +"Description=\"estimated ALT dose on second haplotype\">";
    private static final String GP_FORMAT = "##FORMAT=<ID=GP,Number=G,Type=Float,"
            + "Description=\"Estimated Genotype Probability\">";

    private static final HashMap<String, String> predefinedMetaInfoLines = predefinedMetaInfoLines();

    private static HashMap<String, String> predefinedMetaInfoLines() {
        HashMap<String, String> map = new HashMap<>();
        map.put("fileformat", "##fileformat=" + VCF_VERSION);
        map.put("FORMAT/GT",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        return map;
    }

    private VcfWriter() {
        // private constructor prevents instantiation
    }

    /**
     * Prints the specified VCF meta-information lines to the specified
     * {@code PrintWriter}.
     * @param vcfHeader VCF meta-information lines and header line
     * @param out the {@code PrintWriter} to which the VCF meta-information
     * lines will be written
     * @throws NullPointerException if
     * {@code (vcfHeader == null) || (out == null)}
     */
    public static void copyMetaInfoLines(VcfHeader vcfHeader, PrintWriter out) {
        for (int j=0, n=vcfHeader.nMetaInfoLines(); j<n; ++j) {
            out.println(vcfHeader.metaInfoLine(j));
        }
    }

    /**
     * Prints the specified VCF meta-information lines to the specified
     * {@code PrintWriter}.  The specified {@code key} must equal:
     * "fileformat", "fileDate", or "FORMAT/GT".
     * @param key VCF meta-information line key
     * @param out the {@code PrintWriter} to which the VCF meta-information
     * lines will be written
     * @throws IllegalArgumentException if {@code key} is not equal to:
     * "fileformat", "fileDate", or "FORMAT/GT"
     * @throws NullPointerException if {@code (key == null) || (out == null)}
     */
    public static void writeMetaInfoLine(String key, PrintWriter out) {
        if (key.equals("fileDate")) {
            out.print(key);
            out.print('=');
            out.println(now());
        }
        String line = predefinedMetaInfoLines.get(key);
        if (line!=null) {
            out.println(line);
        }
        else {
            throw new IllegalArgumentException(key);
        }
    }

    /**
     * Writes VCF meta-information lines for INFO and FORMAT fields that
     * are generated for imputed genetic markers to the
     * specified {@code PrintWriter}.  Meta-information lines for
     * INFO/AF, INFO/DR2, INFO/IMP and FORMAT/DS will be printed
     * @param ap {@code true} if the meta-information lines
     * for the FORMAT/AP1 and FORMAT/AP2 subfields should also be printed and
     * {@code false} otherwise
     * @param gp {@code true} if the meta-information line for the
     * FORMAT/GP subfield should also be printed and {@code false} otherwise
     * @param out the {@code PrintWriter} to which the VCF meta-information
     * lines will be written
     * @throws NullPointerException if {@code (out == null)}
     */
    public static void writeImputatedMarkerMetaInfoLines(boolean ap, boolean gp,
            PrintWriter out) {
        out.println(AF_INFO);
        out.println(DR2_INFO);
        out.println(IMP_INFO);
        out.println(DS_FORMAT);
        if (ap) {
            out.println(AP1_FORMAT);
            out.println(AP2_FORMAT);
        }
        if (gp) {
            out.println(GP_FORMAT);
        }
    }

    /**
     * Prints the VCF header line for the specified list of samples to the
     * specified {@code PrintWriter}.  The line is terminated with the
     * platform end-of-line character or characters.
     * @param sampleIds the list of sample identifiers
     * @param out the {@code PrintWriter} to which the VCF header line
     * will be written
     * @throws NullPointerException if {@code (out == null)}
     * @throws NullPointerException if
     * {@code sampleIds == null}, or if {@code (sampleIds[j] == null)} for any
     * {@code j} satisfying {@code (0 <= j && j < sampleIds.length)}
     */
    public static void writeHeaderLine(String[] sampleIds, PrintWriter out) {
        out.print(VcfHeader.HEADER_PREFIX);
        for (String id : sampleIds) {
            if (id==null) {
                throw new NullPointerException("id==null");
            }
            out.print('\t');
            out.print(id);
        }
        out.println();
    }

    /**
     * Returns a timestamp in the format "yyyyMMdd".
     * @return a timestamp in the format "yyyyMMdd"
     */
    public static String now() {
        String dateFormat = "yyyyMMdd";
        Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat(dateFormat);
        return sdf.format(cal.getTime());
    }

    /**
     * Writes the data in phased genotypes for the specified markers
     * to the specified {@code PrintWriter}.
     * @param phasedTarg the estimated phased genotypes
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param out the {@code PrintWriter} to which VCF records will be written
     * @throws IllegalArgumentException if
     * {@code (phasedTarg.isPhased() == false)}
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > phasedTarg.nMarkers())}
     * @throws NullPointerException if
     * {@code ((phasedTarg == null) || (out == null))}
     */
    public static void appendRecords(GT phasedTarg, int start, int end,
            PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        if (phasedTarg.isPhased()==false) {
            throw new IllegalArgumentException("unphased genotypes");
        }
        Samples samples = phasedTarg.samples();
        Markers markers = phasedTarg.markers();
        VcfRecBuilder[] recBuilders = recBuilders(markers, samples.size(), start, end);
        for (int s=0, n=phasedTarg.nSamples(); s<n; ++s) {
            if (samples.isDiploid(s)) {
                int h1 = s<<1;
                int h2 = h1 | 0b1;
                for (int m=start; m<end; ++m) {
                    int a1 = phasedTarg.allele(m, h1);
                    int a2 = phasedTarg.allele(m, h2);
                    recBuilders[m-start].addSampleData(a1, a2);
                }
            }
            else {
                for (int m=start; m<end; ++m) {
                    int a1 = phasedTarg.allele(m, (s<<1));
                    recBuilders[m-start].addSampleData(a1);
                }
            }
        }
        for (VcfRecBuilder vrb : recBuilders) {
            vrb.writeRec(out);
        }
    }

    /**
     * Writes the data in phased genotypes for the specified markers
     * to the specified {@code PrintWriter}.
     * @param samples the list of samples
     * @param recs the estimated phased genotypes
     * @param out the {@code PrintWriter} to which VCF records will be written
     * @throws IllegalArgumentException if there exists a {@code j} such that
     * {@code ((0 <= j) && (j < recs.length) && (recs[j].isPhased() == false))}
     * @throws NullPointerException if
     * {@code ((samples == null) || (recs == null) || (out == null))}
     * @throws NullPointerException if any element of {@code recs} is
     * {@code null}
     */
    public static void appendRecords(Samples samples, GTRec[] recs,
            PrintWriter out) {
        int nSamples = samples.size();
        for (GTRec rec : recs) {
            if (rec.isPhased()==false || rec.samples().equals(samples)==false) {
                throw new IllegalArgumentException("inconsistent data");
            }
            VcfRecBuilder recBuilder = new VcfRecBuilder(rec.marker(), nSamples);
            for (int s=0; s<nSamples; ++s) {
                if (samples.isDiploid(s)) {
                    int h1 = s<<1;
                    int h2 = h1 | 0b1;
                    int a1 = rec.get(h1);
                    int a2 = rec.get(h2);
                    recBuilder.addSampleData(a1, a2);
                }
                else {
                    int a1 = rec.get(s<<1);
                    recBuilder.addSampleData(a1);
                }
            }
            recBuilder.writeRec(out);
        }
    }

    private static VcfRecBuilder[] recBuilders(Markers markers, int nSamples,
            int start, int end) {
        return IntStream.range(start, end)
                .mapToObj(m -> new VcfRecBuilder(markers.marker(m), nSamples))
                .toArray(VcfRecBuilder[]::new);
    }
}

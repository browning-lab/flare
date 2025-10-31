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

import beagleutil.ChromIds;
import blbutil.Const;
import blbutil.StringUtil;
import blbutil.Utilities;
import ints.IntList;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code MarkerUtils} contains static helper methods for the
 * {@code Marker} class.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class MarkerUtils  {

    private static final String[] EMPTY_ID_ARRAY = new String[0];

    private MarkerUtils() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns the number of distinct genotypes, which equals
     * {@code nAlleles*(1 + nAlleles())/2}.
     * @param nAlleles the number of distinct alleles
     * @return the number of distinct genotypes
     */
    public static int nGenotypes(int nAlleles) {
        return (nAlleles*(nAlleles + 1)) >> 1;
    }

    /**
     * Returns {@code (marker.chromID() + ':' + marker.pos())}.
     * @param marker a marker
     * @return {@code (marker.chromID() + ':' + marker.pos())}
     * @throws NullPointerException if {@code marker == null}
     */
    public static String coordinate(Marker marker) {
        StringBuilder sb = new StringBuilder(marker.chromID());
        sb.append(Const.colon);
        sb.append(marker.pos());
        return sb.toString();
    }

    /**
     * Returns
     * {@code (marker.chromID() + ':' + marker.pos() + ':'
     * + marker.alleles().replace(Const.tab, Const.colon))}.
     * @param marker a marker
     * @return {{@code (marker.chromID() + ':' + marker.pos() + ':'
     * + marker.alleles().replace(Const.tab, Const.colon))}
     * @throws NullPointerException if {@code marker == null}
     */
    public static String coordinateAndAlleles(Marker marker) {
        StringBuilder sb = new StringBuilder(marker.chromID());
        sb.append(Const.colon);
        sb.append(marker.pos());
        sb.append(Const.colon);
        sb.append(marker.alleles().replace(Const.tab, Const.colon));
        return sb.toString();
    }

    /**
     * Returns the list of identifiers in the VCF ID field. An array of
     * length 0 is returned if {@code (marker.hasIdData() == false)}.
     * @param marker a marker
     * @return the list of identifiers in the VCF ID field
     * @throws NullPointerException if {@code marker == null}
     */
    public static String[] ids(Marker marker) {
        String id = marker.id();
        if (id.equals(Const.MISSING_DATA_STRING)) {
            return EMPTY_ID_ARRAY;
        }
        else {
            return StringUtil.getFields(marker.id(), Const.semicolon);
        }
    }

    /**
     * Returns the list of marker alleles in the VCF REF and ALT fields in
     * order of their appearance in the VCF record.
     * @param marker a marker
     * @return the list of marker alleles
     * @throws NullPointerException if {@code marker == null}
     */
    public static String[] alleles(Marker marker) {
        String refAlt = marker.alleles();
        int nAlleles = marker.nAlleles();
        String[] alleles = new String[nAlleles];
        int start = 0;
        int endIndex = refAlt.indexOf(Const.tab);
        alleles[0] = refAlt.substring(start, endIndex);
        for (int j=1; j<alleles.length; ++j) {
            int startIndex = endIndex+1;
            endIndex = refAlt.indexOf(Const.comma, startIndex);
            alleles[j] = (endIndex<0) ? refAlt.substring(startIndex)
                    : refAlt.substring(startIndex, endIndex);
        }
        return alleles;
    }

    /**
     * Returns {@code s.substring(0, maxLength)} if
     * {@code s.length() > maxLength} and returns {@code s} otherwise.
     * @param s a string
     * @param maxLength a maximum length
     * @return {@code s.substring(0, maxLength)} if
     * {@code s.length() > maxLength} and returns {@code s} otherwise
     */
    static String truncate(String s, int maxLength) {
        return s.substring(0, Math.min(s.length(), maxLength));
    }

    /**
     * Returns a list of indices of the first 8 tabs (in increasing order) in
     * the specified VCF record
     * @param vcfRec a VCF record
     * @return a list of indices of the first 8 tabs in the specified VCF
     * record
     * @throws IllegalArgumentException if the specified VCF record does
     * not contain at least 8 tab characters
     * @throws NullPointerException if {@code vcfRec == null}
     */
    static IntList first8TabIndices(String vcfRec) {
        int nTabs = 8;
        IntList indices = new IntList(nTabs);
        int index = vcfRec.indexOf(Const.tab);
        for (int j=0; j<nTabs && index>=0; ++j) {
            indices.add(index);
            index = vcfRec.indexOf(Const.tab, index+1);
        }
        if (indices.size() < nTabs) {
            String s = "VCF record does not contain at least " + nTabs + " tabs:"
                    + truncate(vcfRec, 200);
            throw new IllegalArgumentException(s);
        }
        return indices;
    }

    /**
     * Returns the index of the specified chromosome
     * @param vcfRec a VCF record
     * @param chrom a chromosome identifier from the VCF record
     * @return the index of the specified chromosome
     * @throws IllegalArgumentException if
     * {@code chrom.isEmpty()}
     * @throws IllegalArgumentException if
     * {@code (chrom.length()==1 && chrom.charAt(0)=='.')}
     * @throws IllegalArgumentException if {@code chrom} contains
     * a whitespace character.
     * @throws IndexOutOfBoundsException if
     * {@code ChromIds.instance().getIndex(chrom) >= Marker.MAX_CHROMOSOMES
     */
    static short chromIndex(String vcfRec, String chrom) {
        if (chrom.isEmpty()
                || (chrom.length()==1 && chrom.charAt(0)==Const.MISSING_DATA_CHAR)) {
            String s = "ERROR: missing chromosome: " + truncate(vcfRec, 80);
            throw new IllegalArgumentException(s);
        }
        for (int j=0, n=chrom.length(); j<n; ++j) {
            char c = chrom.charAt(j);
            if (Character.isWhitespace(c)) {
                String s = "ERROR: CHROM field contains whitespace: "
                        + truncate(vcfRec, 80);
                Utilities.exit(s);
            }
        }
        int chrIndex = ChromIds.instance().getIndex(chrom);
        if (chrIndex>=Marker.MAX_CHROMOSOMES) {
            throw new IndexOutOfBoundsException(String.valueOf(chrIndex));
        }
        return (short) chrIndex;
    }

    /**
     * Returns a sorted {@code String[]} array whose elements are
     * REF and ALT fields for single nucleotide and monomorphic variants.
     * @return an sorted {@code String[]} array whose elements are
     * REF and ALT fields for single nucleotide and monomorphic variants
     */
    static String[] snvPerms() {
        char[] bases = new char[] {'*', 'A', 'C', 'G', 'T'};
        Arrays.sort(bases);
        List<String> perms = new ArrayList<>(120);
        permute(new char[0], bases, perms);
        // Add markers without ALT allele
        perms.add("A\t.");
        perms.add("C\t.");
        perms.add("G\t.");
        perms.add("T\t.");
        String[] array = perms.toArray(new String[0]);
        Arrays.sort(array);
        return array;
    }

    private static void permute(char[] start, char[] end,
            List<String> perms) {
        if (end.length==0 && start[0]!='*') {
            StringBuilder sb = new StringBuilder();
            sb.append(start[0]);
            for (int j=1; j<start.length; ++j) {
                sb.append(j==1 ? Const.tab : Const.comma);
                sb.append(start[j]);
            }
            perms.add(sb.toString());
        }
        else {
            for (int j=0; j<end.length; ++j) {
                char[] newStart = Arrays.copyOf(start, start.length + 1);
                newStart[start.length] = end[j];

                char[] newEnd = new char[end.length - 1];
                if (j>0) {
                    System.arraycopy(end, 0, newEnd, 0, j);
                }
                if (j<newEnd.length) {
                    System.arraycopy(end, j+1, newEnd, j, (newEnd.length - j));
                }
                permute(newStart, newEnd, perms);
            }
        }
    }

    /**
     * Appends the first eight tab-delimited fields of a VCF record for the
     * specified marker (the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO
     * fields) to the specified {@code StringBuilder}.
     * @param marker a marker storing the first 8 VCF fields
     * @param sb the {@code StringBuilder} to be appended
     * @throws NullPointerException if {@code (marker == null) || (sb == null)}
     */
    public static void appendFirst8Fields(Marker marker, StringBuilder sb) {
        sb.append(marker.chromID());
        sb.append('\t');
        sb.append(marker.pos());
        sb.append('\t');
        sb.append(marker.fields());
    }

    /**
     * Appends the first eight tab-delimited fields of a VCF record for the
     * specified marker (the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO
     * fields) to the specified {@code StringBuilder}. The specified INFO/AN
     * and INFO/AC fields will be added to the beginning of the VCF INFO field.
     * If the INFO/AN or iNFO/AC fields exist in {@code this.info()}, the
     * existing INFO/AN or INFO/AC field will be deleted.
     *
     * @param marker a marker storing the first 8 VCF fields
     * @param sb the {@code StringBuilder} to be appended
     * @param an the total number of alleles in called genotypes
     * @param alleleCounts an array of length {@code this.nAlleles()} whose
     * {@code k-th} entry is the allele count in called genotypes
     * for the {@code k}-th allele
     * @throws IllegalArgumentException if
     * {@code (this.nAlleles() != alleleCounts.length)}
     * @throws NullPointerException if
     * {@code ((sb == null) || (alleleCounts == null))}
     */
    public static void appendFirst8Fields(Marker marker, int an,
            int[] alleleCounts, StringBuilder sb) {
        sb.append(marker.chromID());
        sb.append('\t');
        sb.append(marker.pos());
        sb.append('\t');
        String fields = marker.fields();
        int lastTabP1 = fields.lastIndexOf(Const.tab) + 1;
        sb.append(fields, 0, lastTabP1);
        appendInfo(marker, an, alleleCounts, sb);
    }


    private static void appendInfo(Marker marker, int an, int[] alleleCounts, StringBuilder sb) {
        if (marker.nAlleles()!=alleleCounts.length) {
            throw new IllegalArgumentException(Arrays.toString(alleleCounts));
        }
        appendCounts(alleleCounts, sb, an);
        String info = marker.info();
        int start = 0;
        while (start<info.length()) {
            int end = info.indexOf(';', start);
            if (end==-1) {
                end = info.length();
            }
            if (start<end) {
                String field = info.substring(start, end).trim();
                if (field.length()>0
                        && field.equals(".") == false
                        && field.startsWith("AN=")==false
                        && field.startsWith("AC=")==false) {
                    sb.append(';');
                    sb.append(field);
                }
            }
            start = end + 1;
        }
    }

    private static void appendCounts(int[] alleleCounts, StringBuilder sb, int an) {
        sb.append("AN=");
        sb.append(an);
        sb.append(";AC=");
        for (int j=1; j<alleleCounts.length; ++j) {   // begin with first ALT allele
            if (j>1){
                sb.append(',');
            }
            sb.append(alleleCounts[j]);
        }
    }

    /*
     * Prints the first eight tab-delimited fields of a VCF record for the
     * specified marker (the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO
     * fields) to the specified {@code StringBuilder}.
     * @param marker a marker storing the first 8 VCF fields
     * @param out the object to which the marker fields will be printed
     * @throws NullPointerException if {@code (marker == null) || (out == null)}
     */
    public static void printFirst8Fields(Marker marker, PrintWriter out) {
        out.print(marker.chromID());
        out.print('\t');
        out.print(marker.pos());
        out.print('\t');
        out.print(marker.fields());
    }

    /**
     * Prints the first eight tab-delimited fields of a VCF record for this
     * marker (the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO fields)
     * to the specified {@code PrintWriter}.  The specified INFO/AN and
     * INFO/AC fields will be added to the beginning of the VCF INFO field.
     * If the INFO/AN or iNFO/AC fields exist in {@code this.info()}, the
     * existing INFO/AN or INFO/AC field will be deleted.
     *
     * @param marker a marker storing the first 8 VCF fields
     * @param an the total number of alleles in called genotypes
     * @param alleleCounts an array of length {@code this.nAlleles()} whose
     * {@code k-th} entry is the allele count in called genotypes
     * for the {@code k}-th allele
     * @param out the object to which the marker fields will be printed
     * @throws IllegalArgumentException if
     * {@code (this.nAlleles() != alleleCounts.length)}
     * @throws NullPointerException if
     * {@code ((alleleCounts == null)) || (out == null)}
     */
    public static void printFirst8Fields(Marker marker, int an,
            int[] alleleCounts, PrintWriter out) {
        out.print(marker.chromID());
        out.print('\t');
        out.print(marker.pos());
        out.print('\t');
        String fields = marker.fields();
        int lastTabP1 = fields.lastIndexOf(Const.tab) + 1;
        out.write(fields, 0, lastTabP1);
        printInfo(marker, an, alleleCounts, out);
    }

    private static void printInfo(Marker marker, int an, int[] alleleCounts, PrintWriter out) {
        if (marker.nAlleles()!=alleleCounts.length) {
            throw new IllegalArgumentException(Arrays.toString(alleleCounts));
        }
        printCounts(an, alleleCounts, out);
        String info = marker.info();
        int start = 0;
        while (start<info.length()) {
            int end = info.indexOf(';', start);
            if (end==-1) {
                end = info.length();
            }
            if (start<end) {
                String field = info.substring(start, end).trim();
                if (field.length()>0
                        && field.equals(".") == false
                        && field.startsWith("AN=")==false
                        && field.startsWith("AC=")==false) {
                    out.print(';');
                    out.print(field);
                }
            }
            start = end + 1;
        }
    }

    private static void printCounts(int an, int[] alleleCounts, PrintWriter out) {
        out.print("AN=");
        out.print(an);
        out.print(";AC=");
        for (int j=1; j<alleleCounts.length; ++j) {   // begin with first ALT allele
            if (j>1){
                out.append(',');
            }
            out.print(alleleCounts[j]);
        }
    }
}

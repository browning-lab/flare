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
import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
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
     * Returns {@code (marker.chrom() + ':' + marker.pos())}.
     * @param marker a marker
     * @return {@code (marker.chrom() + ':' + marker.pos())}
     * @throws NullPointerException if {@code marker == null}
     */
    public static String coordinate(Marker marker) {
        StringBuilder sb = new StringBuilder(marker.chrom());
        sb.append(Const.colon);
        sb.append(marker.pos());
        return sb.toString();
    }

    /**
     * Returns
     * {@code (marker.chrom() + ':' + marker.pos() + ':'
     * + marker.refAlt().replace(Const.tab, Const.colon))}.
     * @param marker a marker
     * @return {{@code (marker.chrom() + ':' + marker.pos() + ':'
     *          + marker.refAlt().replace(Const.tab, Const.colon))}
     * @throws NullPointerException if {@code marker == null}
     */
    public static String coordinateAndAlleles(Marker marker) {
        StringBuilder sb = new StringBuilder(marker.chrom());
        sb.append(Const.colon);
        sb.append(marker.pos());
        sb.append(Const.colon);
        sb.append(marker.refAlt().replace(Const.tab, Const.colon));
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
        if (marker.hasIdData()) {
            return StringUtil.getFields(marker.id(), Const.semicolon);
        }
        else {
            return EMPTY_ID_ARRAY;
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
        String refAlt = marker.refAlt();
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
     * Prints the first seven tab-delimited VCF fields for the specified
     * marker to the specified {@code PrintWriter}. The delimiter preceding
     * the 8-th VCF field (the INFO field) is not appended.
     * @param marker a marker
     * @param out the {@code PrintWriter} that will receive output
     * @throws NullPointerException if {@code (marker == null) || (out == null)}
     */
    public static void printFirst7Fields(Marker marker, PrintWriter out) {
        out.print(marker.chrom());
        out.print(Const.tab);
        out.print(marker.pos());
        out.print(Const.tab);
        out.print(marker.id());
        out.print(Const.tab);
        out.print(marker.refAlt());
        out.print(Const.tab);
        out.print(marker.qual());
        out.print(Const.tab);
        out.print(marker.filter());
    }

    /**
     * Appends the first seven tab-delimited VCF fields for the specified
     * marker to the specified {@code StringBuilder}. The delimiter preceding
     * the 8-th field (the INFO field) is not appended.
     * @param marker a marker
     * @param sb the {@code StringBuilder} that will be appended
     * @throws NullPointerException if {@code (marker == null) || (sb == null)}
     */
    public static void appendFirst7Fields(Marker marker, StringBuilder sb) {
        sb.append(marker.chrom());
        sb.append(Const.tab);
        sb.append(marker.pos());
        sb.append(Const.tab);
        sb.append(marker.id());
        sb.append(Const.tab);
        sb.append(marker.refAlt());
        sb.append(Const.tab);
        sb.append(marker.qual());
        sb.append(Const.tab);
        sb.append(marker.filter());
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
            String s = "VCF record does not contain " + nTabs + " tabs:"
                    + truncate(vcfRec, 800);
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
     * {@code ChromIds.instance().getIndex(chrom) >= Short.MAX_VALUE}.
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
        if (chrIndex>=Short.MAX_VALUE) {
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
        // next add monomorphic markers
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

}

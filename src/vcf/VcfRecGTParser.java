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

import blbutil.BitArray;
import blbutil.Const;
import blbutil.StringUtil;
import ints.IntList;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code VcfRecGTParser} parses VCF records and extracts the GT format
 * field.  If one allele in a diploid genotype is missing, then both alleles
 * are set to missing.
 * </p>
 * <p>Instances of class {@code VcfRecGTParser} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRecGTParser {

    private final VcfHeader vcfHeader;
    private final Samples samples;
    private final String vcfRec;
    private final Marker marker;
    private final int nAlleles;
    private final int nSamples;
    private final int ninthTabPos;

    private static int majAllele(IntList[] hapLists) {
        int majAllele = 0;
        for (int j=1; j<hapLists.length; ++j) {
            if (hapLists[j].size()>hapLists[majAllele].size()) {
                majAllele = j;
            }
        }
        return majAllele;
    }

    /**
     * Constructs a new {@code VcfRecGTParser} object from the specified VCF
     * record.
     * @param vcfHeader the VCF meta-information lines and header line
     * @param vcfRec the VCF record
     * @throws IllegalArgumentException if {@code vcfHeader.size() == 0}
     * @throws IllegalArgumentException if a format error is detected in the
     * {@code vcfRecord}
     * @throws NullPointerException if
     * {@code vcfHeader == null || vcfRec == null}
     */
    public VcfRecGTParser(VcfHeader vcfHeader, String vcfRec) {
        if (vcfHeader.nSamples()==0) {
            throw new IllegalArgumentException("nSamples==0");
        }
        this.vcfHeader = vcfHeader;
        this.samples = vcfHeader.samples();
        this.vcfRec = vcfRec;
        this.marker = new BasicMarker(vcfRec);
        this.nAlleles = marker.nAlleles();
        this.nSamples = vcfHeader.nSamples();
        this.ninthTabPos = ninthTabPos(vcfRec);
    }

    static int ninthTabPos(String vcfRec) {
        int pos = -1;
        for (int j=0; j<9; ++j) {
            pos = vcfRec.indexOf(Const.tab, pos + 1);
            if (pos == -1) {
                throw new IllegalArgumentException(
                        "VCF record format error: " + vcfRec);
            }
        }
        return pos;
    }

    /**
     * Returns the VCF meta-information lines and header line for the backing
     * VCF record
     * @return the VCF meta-information lines and header line
     */
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    /**
     * Returns the backing VCF record.
     * @return the backing VCF record
     */
    public String vcfRecord() {
        return vcfRec;
    }

    /**
     * Returns the marker.
     * @return the marker
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns {@code this.marker().nAlleles()}.
     * @return the number of alleles
     */
    public int nAlleles() {
        return nAlleles;
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return vcfHeader.samples();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return nSamples;
    }

    /**
     * Stores the genotypes genotypes in the specified BitLists.  The contract
     * for this method is unspecified if any bit of the specified
     * {@code alleles} and {@code isMissing} parameters is set when this method
     * is invoked.
     * @param alleles a BitArray in which the allele for each haplotype is stored
     * @param isMissing a BitArray whose {@code k}-th bit will be set
     * if any allele of the {@code k}-th sample is missing
     * @return {@code true} if all genotypes are phased and non-missing
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws IllegalArgumentException if
     * {@code isMissing.size() != this.samples().size()}
     * @throws IllegalArgumentException if
     * {@code alleles.size() != 2*this.samples().size()*this.marker.bitsPerAllele()}
     * @throws NullPointerException if
     * {@code alleles == null || isMissing == null}
     */
    public boolean storeAlleles(BitArray alleles, BitArray isMissing) {
        this.samples.size();
        int bitsPerAllele = marker.bitsPerAllele();
        if (isMissing.size()!=nSamples) {
            throw new IllegalArgumentException(String.valueOf(isMissing.size()));
        }
        if (alleles.size()!=((nSamples<<1)*bitsPerAllele)) {
            throw new IllegalArgumentException(String.valueOf(alleles.size()));
        }
        int pos = ninthTabPos;
        int unfilt = -1;
        boolean isPhased = true;
        for (int s=0; s<nSamples; ++s) {
            if (pos == -1) {
                throwFieldCountError(vcfHeader, vcfRec);
            }
            int nextUnfiltered = vcfHeader.unfilteredSampleIndex(s);
            while (++unfilt < nextUnfiltered) {
                pos = vcfRec.indexOf(Const.tab, pos + 1);
                if (pos == -1) {
                    throwFieldCountError(vcfHeader, vcfRec);
                }
            }
            int alStart = pos+1;
            int alEnd1 = alEnd1(vcfRec, alStart);
            if (alStart==alEnd1) {
                throwIllegalArgException("missing data", s, false);
            }
            int alEnd2 = alEnd2(vcfRec, alEnd1);
            boolean isDiploid = alEnd1!=alEnd2;
            if (isDiploid!=samples.isDiploid(s)) {
                haploidDiploidError(s, false, isDiploid);
            }
            int a1 = parseAllele(alStart, alEnd1);
            int a2 =  alEnd1==alEnd2 ? a1 : parseAllele(alEnd1 + 1, alEnd2);
            if (isDiploid && vcfRec.charAt(alEnd1)==Const.unphasedSep) {
                isPhased = false;
            }
            if (a1 == -1 || a2 == -1) {
                isPhased = false;
                isMissing.set(s);
            }
            else {
                int h1 = s<<1;
                int h2 = h1 | 0b1;
                storeAllele(alleles, h1, bitsPerAllele, a1);
                storeAllele(alleles, h2, bitsPerAllele, a1);
            }
            pos = vcfRec.indexOf(Const.tab, alEnd2);
        }
        return isPhased;
    }

    /**
     * Stores the genotypes and per-genotype phase data in the specified arrays
     * and returns {@code} true if all genotypes are phased and non-missing.
     * @param alleles an array in which alleles will be stored
     * @param isPhased a boolean array whose {@code k}-th element is
     * {@code true} if the genotype of the {@code k}-th sample is haploid
     * uses the phased allele separator
     * @return {@code true} if all genotypes are phased and non-missing
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws IllegalArgumentException if
     * {@code alleles.length != 2*this.samples().size()}
     * @throws IllegalArgumentException if
     * {@code isPhased.length != this.samples().size()}
     * @throws NullPointerException if any parameter is
     * {@code alleles == null || isPhased == null}
     */
    public boolean storeAlleles(int[] alleles, boolean[] isPhased) {
        if (alleles.length != (nSamples<<1)) {
            throw new IllegalArgumentException(String.valueOf(alleles.length));
        }
        if (isPhased.length != nSamples) {
            throw new IllegalArgumentException(String.valueOf(isPhased.length));
        }
        int pos = ninthTabPos;
        int unfilt = -1;
        boolean allPhased = true;
        for (int s=0; s<nSamples; ++s) {
            if (pos == -1) {
                throwFieldCountError(vcfHeader, vcfRec);
            }
            int nextUnfiltered = vcfHeader.unfilteredSampleIndex(s);
            while (++unfilt < nextUnfiltered) {
                pos = vcfRec.indexOf(Const.tab, pos + 1);
                if (pos == -1) {
                    throwFieldCountError(vcfHeader, vcfRec);
                }
            }
            int alStart = pos+1;
            int alEnd1 = alEnd1(vcfRec, alStart);
            if (alStart==alEnd1) {
                throwIllegalArgException("missing data", s, false);
            }
            int alEnd2 = alEnd2(vcfRec, alEnd1);
            boolean isDiploid = alEnd1!=alEnd2;
            if (isDiploid!=samples.isDiploid(s)) {
                haploidDiploidError(s, false, isDiploid);
            }
            int h1 = s<<1;
            int h2 = h1 | 0b1;
            int a1 = parseAllele(alStart, alEnd1);
            int a2 = alEnd1==alEnd2 ? a1 : parseAllele(alEnd1 + 1, alEnd2);
            if ((a1 == -1)^(a2 == -1)) {
                a1 = -1;
                a2 = -1;
            }
            alleles[h1] = a1;
            alleles[h2] = a2;
            isPhased[s] = isDiploid==false || vcfRec.charAt(alEnd1)==Const.phasedSep;
            allPhased &= isPhased[s];
            pos = vcfRec.indexOf(Const.tab, alEnd2);
        }
        return allPhased;
    }

    public HapListRep hapListRep() {
        IntList[] hapLists = IntStream.range(0, nAlleles)
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        IntList missList = new IntList();
        boolean isPhased = true;
        int tabIndex = ninthTabPos;
        int unfilt = -1;
        for (int s=0; s<nSamples; ++s) {
            if (tabIndex == -1) {
                throwFieldCountError(vcfHeader, vcfRec);
            }
            int nextUnfiltered = vcfHeader.unfilteredSampleIndex(s);
            while (++unfilt < nextUnfiltered) {
                tabIndex = vcfRec.indexOf(Const.tab, tabIndex + 1);
                if (tabIndex == -1) {
                    throwFieldCountError(vcfHeader, vcfRec);
                }
            }
            int alStart = tabIndex+1;
            int alEnd1 = alEnd1(vcfRec, alStart);
            if (alStart==alEnd1) {
                throwIllegalArgException("missing data", s, false);
            }
            int alEnd2 = alEnd2(vcfRec, alEnd1);
            boolean isDiploid = alEnd1!=alEnd2;
            isPhased &= (isDiploid==false || vcfRec.charAt(alEnd1)==Const.phasedSep);
            if (isDiploid!=samples.isDiploid(s)) {
                haploidDiploidError(s, false, isDiploid);
            }
            int h1 = s<<1;
            int h2 = h1 | 0b1;
            int a1 = parseAllele(alStart, alEnd1);
            int a2 =  alEnd1==alEnd2 ? a1 : parseAllele(alEnd1 + 1, alEnd2);
            if (a1<0 || a2<0) {
                missList.add(s);
                isPhased = false;
            }
            else {
                hapLists[a1].add(h1);
                hapLists[a2].add(h2);
            }
            tabIndex = vcfRec.indexOf(Const.tab, alEnd2);
        }
        return new HapListRep(this, missList, hapLists, isPhased);
    }

    private static void storeAllele(BitArray alleles, int hap, int bitsPerAllele,
            int allele) {
        int index = hap*bitsPerAllele;
        int mask = 1;
        for (int k=0; k<bitsPerAllele; ++k) {
            if ((allele & mask)==mask) {
                alleles.set(index);
            }
            ++index;
            mask <<= 1;
        }
    }

    /**
     * Returns the list of phased alleles in the backing VCF record.
     * @return the list of phased alleles in the backing VCF record
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased or missing genotype
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     */
    private int[] phasedAlleles() {
        int[] alleles = new int[2*nSamples];
        int tabIndex = ninthTabPos;
        int unfilt = -1;
        for (int s=0, hap=0; s<nSamples; ++s) {
            if (tabIndex == -1) {
                throwFieldCountError(vcfHeader, vcfRec);
            }
            int nextUnfiltered = vcfHeader.unfilteredSampleIndex(s);
            while (++unfilt < nextUnfiltered) {
                tabIndex = vcfRec.indexOf(Const.tab, tabIndex + 1);
                if (tabIndex == -1) {
                    throwFieldCountError(vcfHeader, vcfRec);
                }
            }
            int alStart = tabIndex+1;
            int alEnd1 = alEnd1(vcfRec, alStart);
            if (alStart==alEnd1) {
                throwIllegalArgException("missing data", s, true);
            }
            int alEnd2 = alEnd2(vcfRec, alEnd1);
            boolean isDiploid = alEnd1!=alEnd2;
            int a1 = parseAllele(alStart, alEnd1);
            int a2 =  alEnd1==alEnd2 ? a1 : parseAllele(alEnd1 + 1, alEnd2);
            if (isDiploid!=samples.isDiploid(s)) {
                haploidDiploidError(s, true, isDiploid);
            }
            if ((isDiploid && vcfRec.charAt(alEnd1)!=Const.phasedSep)
                    || (a1==-1) || (a2==-1)) {
                throwIllegalArgException("unphased or missing genotype", s, true);
            }
            alleles[hap++] = a1;
            alleles[hap++] = a2;
            tabIndex = vcfRec.indexOf(Const.tab, alEnd2);
        }
        return alleles;
    }

    /* returns exclusive end */
    private static int alEnd1(String rec, int start) {
        if (start==rec.length()) {
            throwGTFormatError(rec, rec.length());
        }
        int index = start;
        while (index < rec.length()) {
            char c = rec.charAt(index);
            if (c == Const.unphasedSep || c == Const.phasedSep
                    || c == Const.tab || c == Const.colon) {
                return index;
            }
            ++index;
        }
        return index;
    }

    /* returns exclusive end */
    private static int alEnd2(String rec, int start) {
        int index = start;
        while (index < rec.length()) {
            char c = rec.charAt(index);
            if (c == Const.colon || c == Const.tab) {
                return index;
            }
            ++index;
        }
        return index;
    }

    private int haploidDiploidError(int sample, boolean isRef, boolean isDiploid) {
        String err = (isRef ? "Reference sample " : "Sample ")
                + vcfHeader.samples().id(sample)
                + " has an inconsistent number of alleles. "
                + "The first genotype is "
                + (samples.isDiploid(sample) ? "diploid" : "haploid")
                + ", but the genotype at position "
                + marker.chrom() + ":" + marker.pos()
                + " is " + (isDiploid ? "diploid" : "haploid");
        throw new IllegalArgumentException(err);
    }

    private int throwIllegalArgException(String msg, int sample, boolean isRef) {
        String err = "ERROR: " + msg
                + (isRef ? " for reference sample " : " for sample ")
                + vcfHeader.samples().id(sample)
                + " at marker [" + marker + "]";
        throw new IllegalArgumentException(err);
    }

    private int parseAllele(int start, int end) {
        if (start==end) {
            String s = "ERROR: Missing sample allele: " + vcfRec;
            throw new IllegalArgumentException(s);
        }
        int al;
        if (start + 1 == end) {
            char c = vcfRec.charAt(start);
            if (c=='.') {
                return -1;
            }
            else {
                al = (c - '0');
            }
        }
        else {
            al = Integer.parseInt(vcfRec.substring(start, end));
        }
        if (al < 0 || al >= nAlleles) {
            String strAllele = vcfRec.substring(start, end);
            String s = "ERROR: Invalid allele [" + strAllele + "] at character "
                    + start + " in record \"" + marker + "\t...\"";
            throw new IllegalArgumentException(s);
        }
        return al;
    }

    private static void throwGTFormatError(String rec, int index) {
        StringBuilder sb = new StringBuilder(1000);
        sb.append("ERROR: genotype is missing allele separator:");
        sb.append(Const.nl);
        sb.append(rec.substring(0, index));
        sb.append(Const.nl);
        sb.append("Exiting Program");
        sb.append(Const.nl);
        throw new IllegalArgumentException(sb.toString());
    }

    private static void throwFieldCountError(VcfHeader vcfHeader, String vcfRec) {
        String src = vcfHeader.src();
        String[] fields = StringUtil.getFields(vcfRec, Const.tab);
        StringBuilder sb = new StringBuilder(1000);
        sb.append("ERROR: CF header line has ");
        sb.append(vcfHeader.nHeaderFields());
        sb.append(" fields, but data line has ");
        sb.append(fields.length);
        sb.append(" fields");
        sb.append(Const.nl);
        sb.append("File source: ");
        sb.append(src);
        sb.append(Const.nl);
        sb.append(Arrays.toString(fields));
        sb.append(Const.nl);
        throw new IllegalArgumentException(sb.toString());
    }

    /**
     * Returns an array of length {@code this.nAlleles()} whose
     * {@code k}-th element is the list of haplotype indices carrying
     * the {@code k}-th allele if {@code k} is a non-major allele,
     * and whose {@code k}-th element is {@code null} if {@code k} is
     * the major allele.  If there is more than one allele with maximal count,
     * the allele with maximal count having the smallest index is defined to
     * be the major allele.
     * @return the indices of the haplotypes carrying each non-major allele
     * @throws IllegalArgumentException if a format error is detected in
     * the specified VCF record or if the specified VCF header is
     * inconsistent with the specified VCF header.
     *
     * @throws NullPointerException if {@code vcfRec == null || rec == null}
     */
    public int[][] nonMajRefIndices() {
        int[] alleles = phasedAlleles();
        int[] alCnts = new int[nAlleles];
        for (int a : alleles) {
            ++alCnts[a];
        }
        int majAl = 0;
        for (int j=1; j<nAlleles; ++j) {
            if (alCnts[j]>alCnts[majAl]) {
                majAl = j;
            }
        }
        int[][] nonMajIndices = new int[nAlleles][];
        for (int al=0; al<nAlleles; ++al) {
            nonMajIndices[al] = al==majAl ? null : new int[alCnts[al]];
        }
        Arrays.fill(alCnts, 0);
        for (int j=0; j<alleles.length; ++j) {
            int al = alleles[j];
            if (al!=majAl) {
                nonMajIndices[al][alCnts[al]++] = j;
            }
        }
        return nonMajIndices;
    }

    public static class HapListRep {

        private final VcfRecGTParser recParser;
        private final IntList missingSamples;
        private final IntList[] hapLists;
        private final boolean isPhased;
        private final int majorAllele;

        private HapListRep(VcfRecGTParser recParser, IntList missingSamples,
                IntList[] hapLists, boolean isPhased) {
            this.recParser = recParser;
            this.missingSamples = missingSamples;
            this.hapLists = hapLists;
            this.isPhased = isPhased;
            this.majorAllele = majAllele(hapLists);
        }

        /**
         * Returns the non-major allele count.
         * @return the non-major allele count
         */
        public int nonmajorAlleleCnt() {
            int cnt = - hapLists[majorAllele].size();
            for (IntList list : hapLists) {
                cnt += list.size();
            }
            return cnt;
        }

        /**
         * Returns a list of samples with missing indices
         * @return a list of samples with missing indices
         */
        public int[] missingSamples() {
            return missingSamples.toArray();
        }

        /**
         * Returns the major allele. The major allele is the allele
         * with largest number of copies.  If there are multiple major alleles
         * the major allele with smallest index is returned.
         * @return the major allele
         */
        public int majorAllele() {
            return majorAllele;
        }

        /**
         * Returns an array of length {@code this.marker().nAlleles()} that
         * contains lists of haplotypes carrying each allele.
         * The increasing list of haplotype indices carrying non-major allele
         * {@code j} is stored in the {@code j-th} element of the returned array.
         * if (@code setMajorToNull == true} then element of the returned
         * array corresponding to {@code this.majorAllele()} will be {@code null}
         * @param setMajorToNull {@code true} if the element of the returned
         * array corresponding to {@code this.majorAllele()} will be {@code null}
         * @return list of haplotypes carrying each allele.
         */
        public int[][] hapLists(boolean setMajorToNull) {
            int[][] lists = new int[hapLists.length][];
            for (int j=0; j<lists.length; ++j) {
                if (setMajorToNull==false || j!=majorAllele) {
                    lists[j] = hapLists[j].toArray();
                }
            }
            return lists;
        }

        /**
         * Returns the marker.
         * @return the marker
         */
        public Marker marker() {
            return recParser.marker();
        }

        /**
         * Returns the samples.
         * @return the samples
         */
        public Samples samples() {
            return recParser.samples();
        }

        /**
         * Returns {@code true} if all genotypes are phased and there
         * are no missing alleles, and returns {@code false} otherwise.
         * @return {@code true} if all genotypes are phased and nonomissing
         */
        public boolean isPhased() {
            return isPhased;
        }
    }
}

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

/**
 * <p>Class {@code BitArrayGT} represents genotypes for a list of samples
 * at a single marker. Instances of class {@code BitArrayGTRec} store
 * haplotype alleles and flags to indicate missing genotypes in bit sets.  All
 * genotypes are considered to be unphased if any sample has an
 * unphased or missing genotype.t</p>
 *
 * <p>Instances of class {@code BitArrayGTRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BitArrayGTRec implements GTRec {

    private final int bitsPerAllele;
    private final Marker marker;
    private final Samples samples;
    private final boolean isPhased;

    private final BitArray isMissing;
    private final BitArray alleles;

    @Override
    public long estBytes() {
        int overhead = 2*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 4*8;  // assume 8 bytes per reference
        estBytes += isMissing.estBytes();
        estBytes += alleles.estBytes();
        return estBytes;
    }

    /**
     * Constructs a new {@code BitArrayGT} instance representing
     * the specified VCF record's GT format field data.
     *
     * @param recParser the VCF record genotype data
     * @throws IllegalArgumentException if a format error is detected
     * in the VCF record
     * @throws NullPointerException if {@code recParser == null}
     */
    public BitArrayGTRec(VcfRecGTParser recParser) {
        int nSamples = recParser.samples().size();
        int nHaps = nSamples<<1;
        this.bitsPerAllele = recParser.marker().bitsPerAllele();
        this.marker = recParser.marker();
        this.samples = recParser.samples();
        BitArray alleleList = new BitArray(nHaps*bitsPerAllele);
        BitArray isMissingList = new BitArray(nSamples);
        this.isPhased = recParser.storeAlleles(alleleList, isMissingList);
        this.isMissing = isMissingList;
        this.alleles = alleleList;
    }

    /**
     * Constructs a new {@code BitArrayGT} instance representing
     * the specified VCF record's GT format field data.
     *
     * @param hlr the VCF record genotype data
     * @throws IllegalArgumentException if a format error is detected
     * in the VCF record
     * @throws NullPointerException if {@code hlr == null}
     */
    public BitArrayGTRec(VcfRecGTParser.HapListRep hlr) {
        int nSamples = hlr.samples().size();
        int nHaps = nSamples<<1;
        this.bitsPerAllele = hlr.marker().bitsPerAllele();
        this.marker = hlr.marker();
        this.samples = hlr.samples();
        this.isPhased = hlr.isPhased();
        this.isMissing = new BitArray(nSamples);
        this.alleles = new BitArray(nHaps*bitsPerAllele);
        boolean setMajorToNull = false;
        int[][] hapLists = hlr.hapLists(setMajorToNull);
        int[] missingSamples = hlr.missingSamples();
        for (int al=0; al<hapLists.length; ++al) {
            int[] list = hapLists[al];
            for (int h : list) {
                storeAllele(alleles, h, bitsPerAllele, al);
            }
        }
        for (int s : missingSamples) {
            isMissing.set(s);
        }
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

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int size() {
        return 2*samples.size();
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public boolean isPhased() {
        return isPhased;
    }


    @Override
    public boolean isPhased(int sample) {
        return isPhased;
    }

    @Override
    public int get(int hap) {
        return isMissing.get(hap>>1) ? -1 : allele(hap);
    }

    private int allele(int hap) {
        int start = bitsPerAllele*hap;
        int end = start + bitsPerAllele;
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            if (alleles.get(j)) {
                allele += mask;
            }
            mask <<= 1;
        }
        return allele;
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field. The returned VCF record
     * will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field
     */
    @Override
    public String toString() {
        return GTRec.toVcfRec(this);
    }

//    /**
//     * Constructs a new {@code BitArrayGT} instance representing
//     * the specified VCF record's GT format field data.
//     *
//     * @param rec a VCF file record.
//     * @param fam parent-offspring relationships.
//     * @param usePhase {@code true} if phase information in the specified
//     * VCF file record will be used, and {@code false} if phase
//     * information in the specified VCF file record will be ignored.
//     *
//     * @throws IllegalArgumentException if
//     * {@code rec.size()==0|| rec.samples().equals(fam.samples())==false}.
//     * @throws IllegalArgumentException if the VCF record does not have a
//     * GT format field.
//     * @throws NullPointerException if {@code rec==null || fam==null}.
//     */
//    public BitArrayGT(VcfRecord rec, NuclearFamilies fam, boolean usePhase) {
//        this(rec);
//        if (rec.samples().equals(fam.samples())==false) {
//            throw new IllegalArgumentException("inconsistent samples");
//        }
//        setBits(rec, usePhase, bitsPerAllele, allele1, allele2, isMissing1,
//                isMissing2, isPhased);
//        removeMendelianInconsistencies(rec, fam, isPhased, isMissing1,
//                isMissing2);
//    }
//
//    /*
//     * Sets phase to unknown for all parent-offspring relationships, and sets
//     * all genotypes in a duo or trio genotypes to missing if a Mendelian
//     * inconsistency is found.
//     */
//    private static void removeMendelianInconsistencies(VcfRecord rec,
//            NuclearFamilies fam, BitSet isPhased, BitSet isMissing1,
//            BitSet isMissing2) {
//        for (int j=0, n=fam.nDuos(); j<n; ++j) {
//            int p = fam.duoParent(j);
//            int o = fam.duoOffspring(j);
//            isPhased.clear(p);
//            isPhased.clear(o);
//            if (duoIsConsistent(rec, p, o) == false) {
//                logDuoInconsistency(rec, p, o);
//                isMissing1.set(p);
//                isMissing2.set(p);
//                isMissing1.set(o);
//                isMissing2.set(o);
//            }
//        }
//        for (int j=0, n=fam.nTrios(); j<n; ++j) {
//            int f = fam.trioFather(j);
//            int m = fam.trioMother(j);
//            int o = fam.trioOffspring(j);
//            isPhased.clear(f);
//            isPhased.clear(m);
//            isPhased.clear(o);
//            if (trioIsConsistent(rec, f, m, o) == false) {
//                logTrioInconsistency(rec, f, m, o);
//                isMissing1.set(f);
//                isMissing2.set(f);
//                isMissing1.set(m);
//                isMissing2.set(m);
//                isMissing1.set(o);
//                isMissing2.set(o);
//            }
//        }
//    }
//
//    private static boolean duoIsConsistent(VcfRecord rec, int parent,
//            int offspring) {
//        int p1 = rec.gt(parent, 0);
//        int p2 = rec.gt(parent, 1);
//        int o1 = rec.gt(offspring, 0);
//        int o2 = rec.gt(offspring, 1);
//        boolean alleleMissing = (p1<0 || p2<0 || o1<0 || o2<0);
//        return (alleleMissing || p1==o1 || p1==o2 || p2==o1 || p2==o2);
//    }
//
//    private static boolean trioIsConsistent(VcfRecord rec, int father,
//            int mother, int offspring) {
//        int f1 = rec.gt(father, 0);
//        int f2 = rec.gt(father, 1);
//        int m1 = rec.gt(mother, 0);
//        int m2 = rec.gt(mother, 1);
//        int o1 = rec.gt(offspring, 0);
//        int o2 = rec.gt(offspring, 1);
//        boolean fo1 = (o1<0 || f1<0 || f2<0 || o1==f1 || o1==f2);
//        boolean mo2 = (o2<0 || m1<0 || m2<0 || o2==m1 || o2==m2);
//        if (fo1 && mo2) {
//            return true;
//        }
//        else {
//            boolean fo2 = (o2<0 || f1<0 || f2<0 || o2==f1 || o2==f2);
//            boolean mo1 = (o1<0 || m1<0 || m2<0 || o1==m1 || o1==m2);
//            return (fo2 && mo1);
//        }
//    }
//
//    private static void logDuoInconsistency(VcfRecord rec, int parent,
//            int offspring) {
//        StringBuilder sb = new StringBuilder(80);
//        sb.append("WARNING: Inconsistent duo genotype set to missing");
//        sb.append(Const.tab);
//        sb.append(rec.marker());
//        sb.append(Const.colon);
//        sb.append(rec.samples().id(parent));
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(offspring));
//        main.Logger.getInstance().println(sb.toString());
//    }
//
//    private static void logTrioInconsistency(VcfRecord rec, int father,
//            int mother, int offspring) {
//        StringBuilder sb = new StringBuilder(80);
//        sb.append("WARNING: Inconsistent trio genotype set to missing");
//        sb.append(Const.tab);
//        sb.append(rec.marker());
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(father));
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(mother));
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(offspring));
//        main.Logger.getInstance().println(sb.toString());
//    }
//
//    /* Following code is copied from vcf.XBasicGT1 class thaat stored
//       unphased genotypes for a sample [XBasicGT1 deleted on 01Feb21] */
//    private boolean useConstraintAndGL(GT gl, int sample, int father,
//            int mother) {
//        int nMarkers = markers.nMarkers();
//        int phasedCnt = 0;
//        for (int m=0; m<nMarkers; ++m) {
//            int a1 = gl.allele1(m, sample);
//            int a2 = gl.allele2(m, sample);
//            int tr1 = transmittedAllele(gl, m, father);
//            int tr2 = transmittedAllele(gl, m, mother);
//            if ((tr1>=0 || tr2>=0) && isUnphasedConsistent(a1, a2, tr1, tr2)) {
//                if (isPhasedConsistent(a1, a2, tr1, tr2)==false) {
//                    int tmp = a1;
//                    a1 = a2;
//                    a2 = tmp;
//                }
//                if (tr1 != -1) {
//                    a1 = tr1;
//                }
//                if (tr2 != -1) {
//                    a2 = tr2;
//                }
//                if (a1>=0 && a2>=0) {
//                    ++phasedCnt;
//                }
//            }
//            setBits(markers, m, a1, allele1, missing1);
//            setBits(markers, m, a2, allele2, missing2);
//        }
//        return phasedCnt==nMarkers;
//}
//
//    private static boolean isUnphasedConsistent(int a1, int a2, int b1, int b2) {
//        return isPhasedConsistent(a1, a2, b1, b2)
//                || isPhasedConsistent(a1, a2, b2, b1);
//    }
//
//    private static boolean isPhasedConsistent(int a1, int a2, int b1, int b2) {
//        return (isConsistent(a1, b1) && isConsistent(a2, b2));
//    }
//
//    private static boolean isConsistent(int a1, int b1) {
//        return (a1==-1 || b1==-1 || a1==b1);
//    }
//
//    private static int transmittedAllele(GT gl, int marker, int sample) {
//        if (sample==-1) {
//            return -1;
//        }
//        int a1 = gl.allele1(marker, sample);
//        return a1>=0 && a1==gl.allele2(marker, sample) ? a1 : -1;
//    }
//
//    private static void setBits(Markers markers, int marker, int allele,
//            BitSet alleles, BitSet missing) {
//        if (allele == -1) {
//            missing.set(marker);
//        }
//        else {
//            int mask = 1;
//            int start = markers.sumHapBits(marker);
//            int end = markers.sumHapBits(marker+1);
//            for (int i=start; i<end; ++i) {
//                boolean b = (allele & mask)==mask;
//                alleles.set(i, b);
//                mask <<= 1;
//            }
//        }
//    }
}

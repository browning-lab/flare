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

import ints.IntArray;
import ints.IntList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.Samples;
import vcf.MapRefGTRec;

/**
 * <p>Class {@code SeqCoder3} compresses a sequence of allele-coded
 * {@code RefGTRec} objects.  The class is designed for use with bref v3 format.
 * Compression is performed by storing the list of distinct allele sequences
 * and the allele sequence carried by each haplotype.
 * </p>
 * <p>Class {@code SeqCoder3} is not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SeqCoder3 {

    /**
     * The maximum number of alleles are that permitted in order
     * for the {@code add()} method to return {@code true}
     */
    public static final int MAX_NALLELES = 255;

    /**
     * The major allele frequency threshold for allele coding.
     * Sequence coding should only be applied if the major allele frequency
     * less than or equal to this threshold.
     */
    public static final float COMPRESS_FREQ_THRESHOLD = 0.995f;

    private final Samples samples;
    private final int maxNSeq;
    private final List<RefGTRec> recs;
    private final int[] hap2Seq;
    private final IntList seq2Cnt;  // for processiing allele-coded records
    private final List<IntList> seq2AlleleSeqMap;

    /**
     * Constructs a new {@code SeqCoder3} for the specified samples.
     * @param samples the list of samples whose data will be compressed
     * @throws NullPointerException if {@code samples == null}
     */
    public SeqCoder3(Samples samples) {
        this(samples, defaultMaxNSeq(samples.size()));
    }

    /**
     * Constructs a new {@code SeqCoder3} for the specified samples.
     * @param samples the list of samples whose data will be compressed
     * @param maxNSeq the maximum number of distinct allele sequences
     * permitted if the {@code add()} method returns {@code true}
     * @throws NullPointerException if {@code samples == null}
     * @throws IllegalArgumentException
     * {@code maxNSeq < 0 || maxNSeq >= Character.MAX_VALUE}
     */
    public SeqCoder3(Samples samples, int maxNSeq) {
        if (maxNSeq < 0 || maxNSeq >= Character.MAX_VALUE) {
            throw new IllegalArgumentException(String.valueOf(maxNSeq));
        }
        this.samples = samples;
        this.maxNSeq = maxNSeq;
        this.recs = new ArrayList<>(100);
        this.hap2Seq = new int[2*samples.size()];
        this.seq2Cnt = new IntList(3*maxNSeq/2 + 1);
        this.seq2AlleleSeqMap = new ArrayList<>(maxNSeq + 1);
        initialize();
    }

    /**
     * Returns the default maximum number of sequences for the specified
     * number of samples.  The default value is equal to
     * {@code (int) Math.min((long) Math.pow(2, 2*Math.log10(size) + 1),
     * Character.MAX_VALUE)}
     * @param nSamples the number of samples
     * @return the default maximum number of sequences for the specified
     * number of samples
     * @throws IllegalArgumentException if {@code size < 1}
     */
    public static int defaultMaxNSeq(int nSamples) {
        if (nSamples < 1) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        if (nSamples==1) {
            return 3;
        }
        else {
            double exponent = 2*Math.log10(nSamples) + 1;
            long maxNSeq = (long) Math.floor(Math.pow(2.0, exponent));
            if (maxNSeq > Character.MAX_VALUE) {
                return Character.MAX_VALUE;
            }
            else {
                return (int) maxNSeq;
            }
        }
    }

    /**
     * Returns the list of samples whose phased genotype data will be compressed.
     * @return the list of samples whose phased genotype data will be compressed
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns the number of compressed {@code RefGTRec} objects.
     * @return the number of compressed {@code RefGTRec} objects
     */
    public int nRecs() {
        return recs.size();
    }

    /**
     * Returns the maximum number of distinct allele sequences.
     * @return the maximum number of distinct allele sequences
     */
    public int maxNSeq() {
        return maxNSeq;
    }

    /**
     * Attempts to add the specified {@code RefGTRec} object to the list of
     * compressed {@code RefGTRec} objects, and returns {@code true}
     * if the {@code RefGTRec} object was added.
     *
     * @param rec reference genotypes for a marker
     * @return {@code true} if the specified {@code RefGTRec} object was
     * added to the list of compressed markers
     *
     * @throws IllegalArgumentException if
     * {@code rec.samples().equals(this.samples()) == false}
     * @throws IllegalArgumentException if
     * {@code rec.isAlleleCoded() == false}
     * @throws NullPointerException if {@code rec == null}
     */
    public boolean add(RefGTRec rec) {
        if (rec.samples().equals(samples)==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (rec.isAlleleRecord()==false) {
            throw new IllegalArgumentException(rec.getClass().toString());
        }
        boolean success = setAlleleMap(rec);
        if (success) {
            recs.add(rec);
            int nullRow = rec.nullRow();
            for (int a=0, n=rec.marker().nAlleles(); a<n; ++a) {
                if (a!=nullRow) {
                    int nCopies = rec.alleleCount(a);
                    for (int c=0; c<nCopies; ++c) {
                        int h = rec.nonNullRowHap(a, c);
                        int oldSeq = hap2Seq[h];
                        IntList list = seq2AlleleSeqMap.get(oldSeq);
                        int index=0;
                        while (index<list.size() && list.get(index)!=a) {
                            index+=2;
                        }
                        int newSeq = list.get(index+1);
                        if (newSeq != oldSeq) {
                            while (newSeq >= seq2Cnt.size()) {
                                seq2Cnt.add(0);
                            }
                            hap2Seq[h] = newSeq;
                            seq2Cnt.decrementAndGet(oldSeq);
                            seq2Cnt.incrementAndGet(newSeq);
                        }
                    }
                }
            }
        }
        assert seq2Cnt.size()==seq2AlleleSeqMap.size();
        return success;
    }

    private void clearSeq2AlleleMap() {
        for (int j=0, n=seq2AlleleSeqMap.size(); j<n; ++j) {
            seq2AlleleSeqMap.get(j).clear();
        }
    }

    private boolean setAlleleMap(RefGTRec rec) {
        assert seq2Cnt.size()==seq2AlleleSeqMap.size();
        int nStartSeq = seq2Cnt.size();
        int[] seq2NonMajorCnt = new int[nStartSeq];
        clearSeq2AlleleMap();
        int nAlleles = rec.marker().nAlleles();
        int nullRow = rec.nullRow();
        for (int a=0; a<nAlleles; ++a) {
            if (a!=nullRow) {
                int nCopies = rec.alleleCount(a);
                for (int c=0; c<nCopies; ++c) {
                    int h = rec.nonNullRowHap(a, c);
                    int seq = hap2Seq[h];
                    ++seq2NonMajorCnt[seq];
                    IntList list = seq2AlleleSeqMap.get(seq);
                    if (list.isEmpty()) {
                        list.add(a);
                        list.add(seq);
                    }
                    else {
                        int index=0;
                        while (index<list.size() && list.get(index)!=a) {
                            index+=2;
                        }
                        if (index==list.size()) {
                            list.add(a);
                            list.add(seq2AlleleSeqMap.size());
                            seq2AlleleSeqMap.add(new IntList(4));
                        }
                    }
                }
            }
        }
        addMajorAllele(seq2NonMajorCnt, nullRow);
        if (seq2AlleleSeqMap.size() > maxNSeq) {
            seq2AlleleSeqMap.subList(nStartSeq, seq2AlleleSeqMap.size()).clear();
            return false;
        }
        else {
            return true;
        }
    }

    private void addMajorAllele(int[] seq2NonNullRowCnt, int nullRowAllele) {
        for (int seq=0; seq<seq2NonNullRowCnt.length; ++seq) {
            if (seq2NonNullRowCnt[seq] < seq2Cnt.get(seq)) {
                IntList list = seq2AlleleSeqMap.get(seq);
                if (list.isEmpty()) {
                    list.add(nullRowAllele);
                    list.add(seq);
                }
                else {
                    // assign nullRow allele the existing sequence index
                    list.add(list.get(0));
                    assert list.get(1)==seq;
                    list.add(seq2AlleleSeqMap.size());
                    list.set(0, nullRowAllele);
                    seq2AlleleSeqMap.add(new IntList(4));
                }
            }
        }
    }

    /**
     * Returns and clears the stored list of compressed {@code RefGTRec}
     * objects.
     *
     * @return the list of compressed {@code RefGTRec} objects
     */
    public List<RefGTRec> getCompressedList() {
        if (recs.isEmpty()) {
            return Collections.emptyList();
        }
        List<RefGTRec> list = new ArrayList<>(recs.size());
        int[] seq2Hap = seq2FirstHap();
        IntArray hap2seq =  IntArray.packedCreate(hap2Seq, seq2Hap.length);
        for (int j=0, n=recs.size(); j<n; ++j) {
            RefGTRec rec = recs.get(j);
            Marker m = rec.marker();
            IntArray seq2allele = seq2Allele(rec, seq2Hap);
            list.add(new MapRefGTRec(m, samples, hap2seq, seq2allele));
        }
        initialize();
        return list;
    }

    private int[] seq2FirstHap() {
        int[] seqToFirstHap = new int[seq2AlleleSeqMap.size()];
        Arrays.fill(seqToFirstHap, -1);
        for (int h=0; h<hap2Seq.length; ++h) {
            int seq = hap2Seq[h];
            if (seqToFirstHap[seq] == -1) {
                seqToFirstHap[seq] = h;
            }
        }
        return seqToFirstHap;
    }

    private IntArray seq2Allele(RefGTRec rec, int[] seq2Hap) {
        int[] seq2Allele = new int[seq2Hap.length];
        for (int j=0; j<seq2Allele.length; ++j) {
            seq2Allele[j] = rec.get(seq2Hap[j]);
        }
        return IntArray.packedCreate(seq2Allele, rec.marker().nAlleles());
    }

    /**
     * Clears the list of compressed {@code RefGTRec} objects.
     */
    private void initialize() {
        recs.clear();
        seq2Cnt.clear();
        seq2AlleleSeqMap.clear();

        // initialize with empty sequence (seq index is 0)
        Arrays.fill(hap2Seq, 0);
        seq2Cnt.add(hap2Seq.length);
        seq2AlleleSeqMap.add(new IntList(4));
    }
}
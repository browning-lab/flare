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

import blbutil.Const;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.FilterUtils;
import blbutil.InputIt;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import bref.Bref3It;
import bref.SeqCoder3;
import ints.IndexArray;
import ints.WrappedIntArray;
import java.io.Closeable;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import vcf.FilterUtil;
import vcf.GTRec;
import vcf.IntArrayRefGTRec;
import vcf.Marker;
import vcf.RefGT;
import vcf. RefGTRec;
import vcf.RefIt;
import vcf.Samples;

/**
 * <p>Class {@code AdmixReader} reads and merges reference and target
 * VCF files with phased genotype data.</p>
 *
 * <p>Instances of class {@code AdmixReader} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AdmixReader implements Closeable {

    private final Filter<String> gtSampleFilter;
    private final SampleFileIt<RefGTRec> targIt;
    private final SampleFileIt<RefGTRec> refIt;
    private final Samples targSamples;
    private final Samples refSamples;
    private final Samples targRefSamples;
    private final int minMac;
    private RefGTRec nextRefRec;
    private RefGTRec nextTargRec;

    private final List<RefGTRec> lowFreqBuffer;
    private final List<RefGTRec> allRecs;
    private final SeqCoder3 seqCoder;
    private final int maxSeqCodedAlleles;
    private final int maxSeqCodingMajorCnt;

    private int markerCnt = 0;

   /**
     * Constructs and returns a new {@code AdmixReader} instance for the
     * specified command line parameters.
     * @param par the command line parameters
     * @return a new {@code AdmixReader} instance
     * @throws NullPointerException if {@code par == null}
     */
    public static AdmixReader instance(AdmixPar par) {
        AdmixReader reader = new AdmixReader(par);
        boolean success = reader.initRecs();
        if (success==false) {
            String nl = blbutil.Const.nl;
            String s = "Error: Reference and target VCF files have no markers in common."
                    + nl + "  Do the two VCF files share some identical markers that appear"
                    + nl + "  in the same order? Identical markers haave identical CHROM,"
                    + nl + "  POS, REF, and ALT fields.";
            throw new IllegalStateException(s);
        }
        return reader;
    }

    private AdmixReader(AdmixPar par) {
        Filter<String> refSampFilt = FilterUtils.includeFilter(
                AdmixUtils.readMap(par.ref_panel()).keySet()
        );
        boolean includeFilter = par.gt_samples()!=null
                && par.gt_samples().startsWith("^")==false;
        Filter<Marker> markerFilter = FilterUtil.markerFilter(par.excludemarkers());
        this.gtSampleFilter = FilterUtil.sampleFilter(par.gtSamplesFile(),
                includeFilter);
        this.targIt = refIt(par.gt(), gtSampleFilter, markerFilter, par.nthreads());
        this.refIt = refIt(par.ref(), refSampFilt, markerFilter, par.nthreads());
        this.targSamples = targIt.samples();
        this.refSamples = refIt.samples();
        this.targRefSamples = targRefSamples(refSamples, targSamples);
        this.minMac = minMac(par, refSamples);

        this.lowFreqBuffer = new ArrayList<>();
        this.allRecs = new ArrayList<>();
        this.seqCoder = new SeqCoder3(targRefSamples);
        this.maxSeqCodedAlleles = Math.min(seqCoder.maxNSeq(), SeqCoder3.MAX_NALLELES);
        this.maxSeqCodingMajorCnt = maxSeqCodingMajorCnt(targRefSamples);
    }

    private static SampleFileIt<RefGTRec> refIt(File refFile,
            Filter<String> sampleFilter, Filter<Marker> markerFilter,
            int nThreads) {
        String filename = refFile.toString();
        SampleFileIt<RefGTRec> refIt;
        if (filename.endsWith(".bref3")) {
            refIt = new Bref3It(refFile, sampleFilter, markerFilter);
        }
        else {
            int nBufferedBlocks = nThreads << 3;
            FileIt<String> it = InputIt.fromBGZipFile(refFile, nBufferedBlocks);
            refIt = RefIt.create(it, sampleFilter, markerFilter);
        }
        return refIt;
    }

    private static Samples targRefSamples(Samples refSamples, Samples targSamples) {
        int nTargSamples = targSamples.size();
        int nAllSamples = nTargSamples + refSamples.size();
        int[] idIndex = new int[nAllSamples];
        boolean[] isDiploid = new boolean[nAllSamples];
        for (int j=0; j<nTargSamples; ++j) {
            idIndex[j] = targSamples.idIndex(j);
            isDiploid[j] = targSamples.isDiploid(j);
        }
        for (int j=nTargSamples; j<nAllSamples; ++j) {
            int targIndex = j - nTargSamples;
            idIndex[j] = refSamples.idIndex(targIndex);
            isDiploid[j] = refSamples.isDiploid(targIndex);
        }
        checkForDuplicates(refSamples, targSamples);
        return new Samples(idIndex, isDiploid);
    }

    private static void checkForDuplicates(Samples refSamples,
            Samples targSamples) {
        HashSet<String> targIds = new HashSet<>(Arrays.asList(targSamples.ids()));
        String[] refIds = refSamples.ids();
        for (int j=0; j<refIds.length; ++j) {
            if (targIds.contains(refIds[j])) {
                String err = "A sample cannot be both a reference sample and a study sample";
                String info = Const.nl + "Error      :  " + err
                        + Const.nl     + "Sample     :  " + refIds[j];
                Utilities.exit(new Throwable(err), info);
            }
        }
    }

    private static int minMac(AdmixPar par, Samples refSamples) {
        if (par.array()==false && (par.min_mac() >= refSamples.size())) {
            String err = "The min-mac parameter must be less than the number of reference samples";
            String info = Const.nl + "Error       :  " + err
                    + Const.nl     + "min-mac     :  " + par.min_mac()
                    + Const.nl     + "ref samples :  " + refSamples.size();
            Utilities.exit(new Throwable(err), info);
        }
        int nRefHaps = refSamples.size()<<1;
        int mafMinMac = (int) Math.ceil(par.min_maf()*nRefHaps);
        return par.array()==true ? mafMinMac : Math.max(mafMinMac, par.min_mac());
    }

    private int maxSeqCodingMajorCnt(Samples samples) {
        int nHaps = samples.size() << 1;
        return (int) Math.floor(nHaps*SeqCoder3.COMPRESS_FREQ_THRESHOLD);
    }

    private boolean initRecs() {
        if (refIt.hasNext()) {
            nextRefRec = refIt.next();
        }
        else {
            return false;
        }
        if (targIt.hasNext()) {
            nextTargRec = targIt.next();
        }
        else {
            return false;
        }
        advanceRefItToTargChr();
        boolean success = advanceRefItToTargPos();
        if (success==false) {
            success = updateNextRecs();
        }
        return success;
    }

    private void advanceRefItToTargChr() {
        int targChr = nextTargRec.marker().chromIndex();
        while (nextRefRec.marker().chromIndex()!=targChr && refIt.hasNext()) {
            nextRefRec = refIt.next();
        }
    }

    /*
     * Advances next ref record until the ref record marker equals the
     * target record marker, or the ref record chrom differs from
     * target record chrom, or the ref record position is greater than target
     * record position, or there are no more ref records.  Returns
     * {@code true} if a {@coe nextRefRec.marker().equals(nextTargRec.marker())}.
     */
    private boolean advanceRefItToTargPos() {
        Marker targMarker = nextTargRec.marker();
        int targChr = targMarker.chromIndex();
        int targPos = targMarker.pos();
        Marker refMarker = nextRefRec.marker();
        while (refMarker.chromIndex()==targChr
                && (refMarker.pos()<targPos
                        || (refMarker.pos()==targPos && targMarker.equals(refMarker)==false))
                && refIt.hasNext()) {
            nextRefRec = refIt.next();
            refMarker = nextRefRec.marker();
        }
        return targMarker.equals(refMarker);
    }

    private boolean updateNextRecs() {
        boolean success = false;
        int targChr = nextTargRec.marker().chromIndex();
        while (targIt.hasNext() && success==false) {
            nextTargRec = targIt.next();
            if (nextTargRec.marker().chromIndex()!=targChr) {
                targChr = nextTargRec.marker().chromIndex();
                advanceRefItToTargChr();
            }
            success = advanceRefItToTargPos();
        }
        if (success==false) {
            nextRefRec = null;
            nextTargRec = null;
        }
        return success;
    }

    /**
     * Returns the sample filter.
     * @return the sample filter
     */
    public Filter<String> gtSampleFilter() {
        return gtSampleFilter;
    }

    /**
     * Returns the list of reference samples.
     * @return the list of reference samples
     */
    public Samples refSamples() {
        return refSamples;
    }

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    public Samples targSamples() {
        return targSamples;
    }

    /**
     * Returns the list of reference and target samples.
     * The reference samples precede the target samples.
     * @return the list of reference and target samples
     */
    public Samples allSamples() {
        return targRefSamples;
    }

    /**
     * Returns the cumulative number of markers returned by previous
     * invocations of {@code this.nextChrom()}.
     * @return the cumulative number of markers
     */
    public int nMarkers() {
        return markerCnt;
    }

    /**
     * Closes this {@code AdmixReader} and releases system resources held
     * by it. If this {@code AdmixReader} is already closed, then invoking
     * this method has no effect.
     */
    @Override
    public void close() {
        refIt.close();
        targIt.close();
        nextRefRec = null;
        nextTargRec = null;
    }

    /**
     * Reads and returns the phased data for the next chromosome. Returns
     * {@code Optional.empty()} if {@code this.close()} has previously been
     * invoked or if there is no more phased data.
     * @return the phased data for the next chromosome
     */
    public Optional<RefGT> nextChrom() {
        if (nextTargRec==null) {
            return Optional.empty();
        }
        int[] scratch = new int[targRefSamples.size()<<1];
        int chr = nextTargRec.marker().chromIndex();
        boolean success = true;
        while (success==true && nextTargRec.marker().chromIndex()==chr) {
            int[] refAlCnts = nextRefRec.alleleCounts();
            Arrays.sort(refAlCnts);
            int mac = refAlCnts.length==1 ? 0 : refAlCnts[refAlCnts.length-2];
            if (mac>=minMac) {
                RefGTRec combinedRec = combine(nextTargRec, nextRefRec, scratch);
                addToBuffer(combinedRec);
            }
            success = updateNextRecs();
        }
        flushCompressedRecords();
        if (allRecs.isEmpty()) {
            return Optional.empty();
        }
        else {
            RefGT refGT = new RefGT(allRecs.toArray(new RefGTRec[0]));
            allRecs.clear();
            markerCnt += refGT.nMarkers();
            return Optional.of(refGT);
        }
    }

    private RefGTRec combine(RefGTRec nextTargRec, RefGTRec nextRefRec,
            int[] alleles) {
        int nTargHaps = nextTargRec.size();
        int nRefHaps = nextRefRec.size();
        for (int h=0; h<nTargHaps; ++h) {
            alleles[h] = nextTargRec.get(h);
        }
        for (int h=0; h<nRefHaps; ++h) {
            alleles[nTargHaps+h] = nextRefRec.get(h);
        }
        Marker mkr = nextRefRec.marker();
        IndexArray ia = new IndexArray(new WrappedIntArray(alleles), mkr.nAlleles());
        return new IntArrayRefGTRec(mkr, targRefSamples, ia);
    }

    private void addToBuffer(RefGTRec rec) {
        if (lowFreqBuffer.size()==Integer.MAX_VALUE) {
           flushCompressedRecords();
        }
        rec = RefGTRec.alleleRefGTRec(rec);
        int[] alCnts = rec.alleleCounts();
        if (alCnts[rec.majorAllele()]>maxSeqCodingMajorCnt
               || alCnts.length>maxSeqCodedAlleles) {
            lowFreqBuffer.add(rec);
        }
        else {
           boolean success = seqCoder.add(rec);
           if (success == false) {
               flushCompressedRecords();
               success = seqCoder.add(rec);
               assert success;
           }
           lowFreqBuffer.add(null);
       }
    }

    private void flushCompressedRecords() {
        List<RefGTRec> list = seqCoder.getCompressedList();
        int index = 0;
        for (int j=0, n=lowFreqBuffer.size(); j<n; ++j) {
            GTRec rec = lowFreqBuffer.get(j);
            if (rec==null) {
                lowFreqBuffer.set(j, list.get(index++));
            }
        }
        allRecs.addAll(lowFreqBuffer);
        lowFreqBuffer.clear();
    }
}

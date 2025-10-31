package admix;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Interface {@code ParamsInterface} represents the analysis parameters
 * for a local ancestry inference analysis.</p>
 *
 * <p>Implementations of interface {@code ParamsInterface} are required to be
 * immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface ParamsInterface {

    /**
     * Returns the reference and target sample metadata.
     * @return the reference and target sample metadata
     */
    SampleData sampleData();

    /**
     * Returns the number of generations to the admixture event.
     * @return the number of generations to the admixture event
     */
    double T();

    /**
     * Returns  an array of length {@code this.sampleData().nAnc()}
     * whose {@code i}-th entry is the proportion of ancestry {@code i}
     * in the target samples.
     * @return an array of length {@code this.sampleData().nAnc()}
     * whose {@code i}-th entry is the proportion of ancestry {@code i}
     * in the target samples
     */
    double[] studyMu();

    /**
     * Returns the copying probability for the specified reference panel
     * conditional on the specified ancestry.
     * @param anc the ancestry index
     * @param panel the reference panel index
     * @return the copying probability for the specified reference panel
     * conditional on the specified ancestry
     * @throws IndexOutOfBoundsException if
     * {@code anc < 0 || anc >= this.sampleData().nAnc()}
     * @throws IndexOutOfBoundsException if
     * {@code panel < 0 || panel >= this.sampleData().nRefPanels()}
     */
    double p(int anc, int panel);

    /**
     * Returns the copying probability for each reference panel conditional
     * on each ancestry.  The returned matrix will have
     * {@code this.sampleData().nAnc()} rows and
     * {@code this.sampleData().nRefPanels()} columns.  The {@code [i][j]}
     * entry of the returned array is the copying probability for
     * reference panel {@code j} conditional on ancestry {@code i}.
     * @return the copying probability for each reference panel conditional
     * on each ancestry
     */
    double[][] p();

    /**
     * Returns the miscopy probabilities for each ancestry and reference panel.
     * The returned array will have {@code this.sampleData().nAnc()} rows and
     * {@code this.sampleData().nRefPanels()} columns.
     * The {@code [i][j]} entry of the returned array is the miscopy
     * probability for ancestry {@code i} and reference panel {@code j}.
     * @return the miscopy probabilities
     */
    double[][] theta();

    /**
     * Returns the intensity of the exponential IBD segment length
     * distribution for IBD segments having the specified ancestry
     * and having a common ancestor before the admixture event.
     * @param anc an ancestry index
     * @return the intensity of the exponential IBD segment length
     * distribution for IBD segments having the specified ancestry
     * and having a common ancestor before the admixture event
     * @throws IndexOutOfBoundsException if
     * {@code anc < 0 || anc >= this.nAnc()}
     */
    double rho(int anc);

    /**
     * Returns an array of length {@code this.sampleData().nAnc()}
     * whose {@code i}-th element is the intensity of the
     * exponential IBD segment length distribution for IBD segments
     * having a common ancestor with ancestry {@code i} before
     * the admixture event.
     *
     * @return an array of length {@code this.sampleData().nAnc()}
     * whose {@code i}-th element is the intensity of the
     * exponential IBD segment length distribution for IBD segments
     * having a common ancestor with ancestry {@code i} before
     * the admixture event.
     */
    double[] rho();

    /**
     * Returns a string representation of {@code this}.  The exact details
     * of the representation are unspecified and subject to change.  However,
     * it is guaranteed that a file consisting of the returned string will
     * be a correctly-formatted analysis parameters file.
     *
     * @return a string representation of {@code this}
     */
    @Override
    String toString();

    /**
     * Returns an array of length {@code nAnc} with each element
     * equal to {@code (1.0/nAnc)}.
     * @param nAnc the number of ancestries
     * @return an array of length {@code nAnc} with each element
     * equal to {@code (1.0/nAnc)}.
     * @throws IllegalArgumentException if {@code (nAnc <= 0)}
     */
    public static double[] defaultMu(int nAnc) {
        if (nAnc<=0) {
            throw new IllegalArgumentException(String.valueOf(nAnc));
        }
        double[] mu = new double[nAnc];
        Arrays.fill(mu, 1.0f/nAnc);
        return mu;
    }

    /**
     * Returns the default {@code theta} values. The returned array will have
     * one row per ancestry and one column per reference panel.
     * @param sampleData reference and target sample metadata
     * @return the default {@code theta} values
     * @throws NullPointerException if {@code (sampleData == null)}
     */
    public static double[][] defaultTheta(SampleData sampleData) {
        int nAnc = sampleData.nAnc();
        int nRefHaps = sampleData.nRefHaps();
        int nRefPanels = sampleData.nRefPanels();
        double theta0 = liStephensPMismatch(nRefHaps);
        double[] row = IntStream.range(0, nRefPanels)
                .mapToDouble(j -> theta0)
                .toArray();
        return IntStream.range(0, nAnc)
                .mapToObj(j -> row)
                .toArray(double[][]::new);
    }

    /**
     * <p>Return an approximation to the allele mismatch probability suggested
     * by Li and Stephens.  The approximation uses a Riemann sum approximation
     * of the natural log function.</p>
     *
     * <p>Refs:
     * Li N, Stephens M. Genetics 2003 Dec;165(4):2213-33 and
     * Marchini J, Howie B. Myers S, McVean G, Donnelly P. 2007;39(7):906-13.</p>
     *
     * @param nHaps the number of haplotypes
     * @return an approximation to the Li and Stephens allele mismatch
     * probability
     * @throws IllegalArgumentException if {@code nHaps < 1}
     */
    public static double liStephensPMismatch(int nHaps) {
        if (nHaps<1) {
            throw new IllegalArgumentException(String.valueOf(nHaps));
        }
        double lambda = 1.0/((Math.log(nHaps) + 0.5));
        return lambda/(2.0*(lambda + nHaps));
    }
}

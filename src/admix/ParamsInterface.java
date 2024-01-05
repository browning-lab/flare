package admix;

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
     * Returns the fixed parameters.
     * @return the fixed parameters
     */
    FixedParams fixedParams();

    /**
     * Returns the number of generations to the admixture event.
     * @return the number of generations to the admixture event
     */
    double T();

    /**
     * Returns  an array of length {@code this.fixedParams().nAnc()}
     * whose {@code i}-th entry is the proportion of ancestry {@code i}
     * in the target samples.
     * @return an array of length {@code this.fixedParams().nAnc()}
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
     * {@code anc < 0 || anc >= this.fixedParams().nAnc()}
     * @throws IndexOutOfBoundsException if
     * {@code panel < 0 || panel >= this.fixedParams().nRefPanels()}
     */
    double p(int anc, int panel);

    /**
     * Returns the copying probability for each reference panel conditional
     * on each ancestry.  The returned matrix will have
     * {@code this.fixedParams().nAnc()} rows and
     * {@code this.fixedParams().nRefPanels()} columns.  The {@code [i][j]}
     * entry of the returned array is the copying probability for
     * reference panel {@code j} conditional on ancestry {@code i}.
     * @return the copying probability for each reference panel conditional
     * on each ancestry
     */
    double[][] p();

    /**
     * Returns the miscopy probabilities for each ancestry and reference panel.
     * The returned array will have {@code this.fixedParams().nAnc()} rows and
     * {@code this.fixedParams().nRefPanels()} columns.
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
     * Returns an array of length {@code this.fixedParams().nAnc()}
     * whose {@code i}-th element is the intensity of the
     * exponential IBD segment length distribution for IBD segments
     * having a common ancestor with ancestry {@code i} before
     * the admixture event.
     *
     * @return an array of length {@code this.fixedParams().nAnc()}
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
}

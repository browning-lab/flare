package admix;

import blbutil.Const;
import ints.IntArray;

/**
 * <p>Class {@code ParamUtils} contains static methods for
 * calculating values from {@code ParamsInterface} instances.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ParamUtils {

    private ParamUtils() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns an array with length {@code params.mu().length} and whose
     * {@code i}-th entry is
     * {@code (1.0 - params.mu()[i]) == 0.0 ? 0.0 : 1.0/(1.0 - params.mu()[i])}.
     *
     * @param params the parameters for a local ancestry analysis
     * @return array with length {@code params.mu().length} and whose
     * {@code i}-th entry is
     * {@code (1.0 - params.mu()[i]) == 0.0 ? 0.0 : 1.0/(1.0 - params.mu()[i])}
     * @throws NullPointerException if {@code params == null}
     */
    public static double[] inv1Mmu(ParamsInterface params) {
        double[] modMu = params.studyMu();
        for (int i=0; i<modMu.length; ++i) {
            double den = 1.0 - modMu[i];
            modMu[i] = den==0.0 ? 0.0 : 1.0/den;
        }
        return modMu;
    }

    /**
     * Returns an array with the same size as {@code params.p()} and whose
     * {@code [i][j]}-th entry is
     * {@code params.p()[i][j] * (1.0/params.fixedParams().nRefPanelHaps().get(j))}.
     * @param params the parameters for a local ancestry analysis
     *
     * @return an array with the same size as {@code params.p()} and whose
     * {@code [i][j]}-th entry is
     * {@code params.p()[i][j] * (1/params.fixedParams().nRefPanelHaps().get(j))}
     * @throws NullPointerException if {@code  params == null}
     */
    public static double[][] q(ParamsInterface params) {
        IntArray nPanelHaps = params.fixedParams().nPanelHaps();
        double[][] modP = params.p();
        double[] invNPanelHaps = new double[nPanelHaps.size()];
        for (int j=0, n=nPanelHaps.size(); j<n; ++j) {
            invNPanelHaps[j] = 1.0/nPanelHaps.get(j);
        }
        for (int i=0; i<modP.length; ++i) {
            for (int j=0; j<invNPanelHaps.length; ++j) {
                modP[i][j] *=invNPanelHaps[j];
            }
        }
        return modP;
    }

    /**
     * Returns a two-dimensional array obtained that is obtained by taking
     * the matrix product of {@code ParamUtils.q(params)}
     * and {@code params.mu()}.
     * The {@code (i, j)}-th element of the returned array is
     * {@code params.mu()[i]*params.p()[i][j]/params.fixedParams().nPanelHaps()[j]}
     * @param params the parameters for a local ancestry analysis
     *
     * @return a two-dimensional array obtained that is obtained by taking
     * the matrix product of {@code ParamUtils.q(params)}
     * and {@code params.mu()}
     * @throws NullPointerException if {@code params == null}
     */
   public  static double[][] qMu(ParamsInterface params) {
        IntArray nPanelHaps = params.fixedParams().nPanelHaps();
        double[][] modP = params.p();
        double[] mu = params.studyMu();
        double[] invNPanelHaps = new double[nPanelHaps.size()];
        for (int j=0, n=nPanelHaps.size(); j<n; ++j) {
            invNPanelHaps[j] = 1.0/nPanelHaps.get(j);
        }
        for (int i=0; i<modP.length; ++i) {
            for (int j=0; j<invNPanelHaps.length; ++j) {
                modP[i][j] *= (mu[i]*invNPanelHaps[j]);
            }
        }
        return modP;
    }

    /**
     * Returns an array with the same dimensions as the {@code params.p()}
     * array and whose {@code [i][j]}-th entry is
     * {@code 1.0/(1.0 - ParamUtils.q(params)[i][j])}.
     * @param params the parameters for a local ancestry analysis
     *
     * @return an array with the same dimensions as the {@code params.p()}
     * array and whose {@code [i][j]}-th entry is
     * {@code 1.0/(1.0 - ParamUtils.q(params)[i][j])}
     * @throws NullPointerException if {@code params == null}
     */
    public static double[][] inv1Mq(ParamsInterface params) {
        double[][] modP = params.p();
        IntArray nPanelHaps = params.fixedParams().nPanelHaps();
        for (int i=0; i<modP.length; ++i) {
            for (int j=0, n=nPanelHaps.size(); j<n; ++j) {
                modP[i][j] = nPanelHaps.get(j)/(nPanelHaps.get(j) - modP[i][j]);
            }
        }
        return modP;
    }

    /**
     * Returns a three dimensional array whose whose {@code [i][j]}-th entry is
     * equal to
     * {@code new double[] {1.0 - params.theta()[i][j], params.theta()[i][j]}}.
     * @param params the parameters for a local ancestry analysis
     *
     * @return a three dimensional array whose whose {@code [i][j]}-th entry is
     * equal to
     * {@code new double[] {1.0 - params.theta()[i][j], params.theta()[i][j]}}
     * @throws NullPointerException if {@code params == null}
     */
    public static double[][][] pObserved(ParamsInterface params) {
        int nAnc = params.fixedParams().nAnc();
        int nRefPanels = params.fixedParams().nRefPanels();
        double[][] theta = params.theta();
        double[][][] mismatchProbs = new double[nAnc][nRefPanels][];
        for (int i=0; i<nAnc; ++i) {
            for (int j=0; j<nRefPanels; ++j) {
                double p = theta[i][j];
                mismatchProbs[i][j] = new double[] {1.0 - p, p};
            }
        }
        return mismatchProbs;
    }

    /**
     * Returns a string representation of the specified {@code params}.
     * The exact details of the representation are unspecified and
     * subject to change.  However, it is guaranteed that a file
     * consisting of the returned string will be a correctly-formatted
     * analysis parameters file.
     * @param params the parameters for a local ancestry analysis
     *
     * @return a string representation of the specified {@code params}
     * @throws NullPointerException if {@code params == null}
     */
    public static String toString(ParamsInterface params) {
        FixedParams fixedParams = params.fixedParams();
        StringBuilder sb = new StringBuilder(1<<9);
        sb.append("# list of ancestries");
        for (int i=0, n=fixedParams.nAnc(); i<n; ++i) {
            sb.append(i==0 ? Const.nl : Const.tab);
            sb.append(fixedParams.ancId(i));
        }
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("# list of reference panels");
        for (int j=0, n=fixedParams.nRefPanels(); j<n; ++j) {
            sb.append(j==0 ? Const.nl : Const.tab);
            sb.append(fixedParams.refPanelId(j));
        }
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("# T: number of generations since admixture");
        sb.append(Const.nl);
        sb.append((float) params.T());
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("# mu[i]: proportion of target genotypes with ancestry i");
        sb.append(Const.nl);
        AdmixUtils.appendLineWithVector(sb, params.studyMu());
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("# p[i][j]: probability that a model state haplotype is in reference panel j");
        sb.append(Const.nl);
        sb.append("#          when the model state ancestry is i");
        sb.append(Const.nl);
        AdmixUtils.appendLinesWithMatrix(sb, params.p());
        sb.append(Const.nl);
        sb.append("# theta[i][j]: probability that a model state haplotype and the target");
        sb.append(Const.nl);
        sb.append("#              haplotype carry different alleles when the model state haplotype");
        sb.append(Const.nl);
        sb.append("#              is in reference panel j and the model state ancestry is i");
        sb.append(Const.nl);
        AdmixUtils.appendLinesWithMatrix(sb, params.theta());
        sb.append(Const.nl);
        sb.append("# rho[i]: rate of the exponential IBD segment cM-length distribution when the");
        sb.append(Const.nl);
        sb.append("#         most recent common ancestor precedes admixture and has ancestry i");
        sb.append(Const.nl);
        AdmixUtils.appendLineWithVector(sb, params.rho());
        return sb.toString();
    }
}

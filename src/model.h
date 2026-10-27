#ifndef model_H
#define model_H
#include <RcppEigen.h>

namespace model{

  // Items, thresholds and cells of every pair, in the order of the frequency table:
  // pairs (k, l < k) by k, cells (sk, sl) by sk
  struct Layout{
    int p, P, nthr, ncell;
    Eigen::VectorXi cat, tk;
    std::vector<int> pk, pl, off;
    Layout(const Eigen::Ref<const Eigen::VectorXd> C_VEC);
  };

  // Model structure, parsed once from the R list: free entries of the pattern
  // matrices and their positions in theta (blocks: thresholds, loadings, latent
  // correlations, log-sds, regressions, residual covariances)
  struct Spec{
    Layout L;
    int p, q, corrflag, nload, ncorr, nvar, nreg, nres, d, ol, oc, ov, ob, orc;
    bool direct;
    Eigen::MatrixXd constrmat, regmat, rescov;
    Eigen::VectorXd constrlogsd;
    std::vector<std::vector<int>> lidx, lfac;
    std::vector<int> ra, ca, sfree, bi, bj;
    Eigen::MatrixXi resid;
    Spec(const Rcpp::List MODEL, const Eigen::Ref<const Eigen::VectorXd> C_VEC);
  };

  // Model quantities at theta: loadings, latent covariance (LKJ factor U and its
  // derivatives alpha, or direct Sig), total effects A = (I - B)^{-1}, reduced
  // loadings Lt = Lam A, residual covariances, and for each item the Jacobian of
  // its reduced loadings with respect to its free loadings and the regressions
  struct Setup{
    const Spec &S;
    Eigen::VectorXd tau, sd;
    Eigen::MatrixXd Lam, Z, C, P, U, alpha, Sig, A, Lt, res;
    std::vector<Eigen::MatrixXd> J;
    Setup(const Eigen::Ref<const Eigen::VectorXd> THETA, const Spec &SPEC);
  };

}

#endif

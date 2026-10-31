#ifndef model_H
#define model_H
#include <RcppEigen.h>

namespace model{

  struct Layout{
    int p, P, nthr, ncell;
    Eigen::VectorXi cat, tk;
    std::vector<int> pk, pl, off;
    Layout(const Eigen::Ref<const Eigen::VectorXd> C_VEC);
  };

  struct Spec{
    Layout L;
    int p, q, corrflag, nload, ncorr, nvar, nreg, nres, nint, nmean, nscale, d, ol, oc, ov, ob, orc, on, oa, os;
    bool direct, thetapar, plain;
    Eigen::MatrixXd constrmat, regmat, rescov;
    Eigen::VectorXd constrlogsd, constrint, constrmean, constrlogscale;
    std::vector<std::vector<int>> lidx, lfac;
    std::vector<int> ra, ca, sfree, bi, bj, nidx, mfree, scidx;
    Eigen::MatrixXi resid;
    Spec(const Rcpp::List MODEL, const Eigen::Ref<const Eigen::VectorXd> C_VEC);
  };

  struct Setup{
    const Spec &S;
    Eigen::VectorXd tau, sd, kappa, mu, v, th, s, taus;
    Eigen::MatrixXd Lam, Z, C, P, U, alpha, Sig, A, Lt, res;
    std::vector<Eigen::MatrixXd> J;
    Setup(const Eigen::Ref<const Eigen::VectorXd> THETA, const Spec &SPEC);
    double rc(const int K, const int L) const { return res.size() > 0 ? res(K, L) : 0; }
  };

}

#endif

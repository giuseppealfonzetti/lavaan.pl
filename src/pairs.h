#ifndef pairs_H
#define pairs_H
#include <RcppEigen.h>

#include "genericUtils.h"
#include "bivariateProbs.h"

namespace pairs{

  struct Setup{
    int p, q, nthr, nload, ncorr, nvar, nreg, nres, m;
    Eigen::MatrixXd Lam, Z, C, P, U, alpha;
    Eigen::VectorXd sd, tau;
    std::vector<int> ra, ca, sfree, bi, bj;
    // Jacobians of the reduced loadings Lt = Lam (I - B)^{-1}, one per item,
    // over the m model parameters following the thresholds
    std::vector<Eigen::MatrixXd> JL;
    // Direct parametrisation of the latent covariances (CORRFLAG = 2)
    bool direct;
    Eigen::MatrixXd Sig;
    // Structural regressions B and residual covariances of the items,
    // with the position of the free ones among the m model parameters
    Eigen::MatrixXd A, Lt, res;
    Eigen::MatrixXi resid;

    Setup(
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const int CORRFLAG = 1,
      const Eigen::Ref<const Eigen::MatrixXd> REGMAT = Eigen::MatrixXd(),
      const Eigen::Ref<const Eigen::MatrixXd> RESCOV = Eigen::MatrixXd());
  };

  struct Pair{
    int k, l, ck, cl, tk, tl, mt, offset;
    double rho;
    Eigen::VectorXd x, y, X, Y, grho;

    Pair(
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const int K,
      const int L,
      const double RHO);

    Pair(
      const Setup &S,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const int K,
      const int L,
      const bool GRADFLAG);
  };

  double pattern_derivatives(
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Pair &PR,
      const int SK,
      const int SL,
      const int ORDER,
      Eigen::VectorXd &G1,
      Eigen::MatrixXd &G2);

  void regression_curvature(const Setup &S, const int K, const Eigen::Ref<const Eigen::VectorXd> W, Eigen::MatrixXd &HRHO);

  Eigen::MatrixXd rho_hessian(const Setup &S, const Pair &PR);

  // Univariate log-likelihood of the observed items, with
  // frequencies UFREQ weighted by the number of missing items
  // of each unit (available-cases term)
  void univariate_contribution(
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> UFREQ,
      const int ORDER,
      double &LL,
      Eigen::VectorXd &GRADIENT,
      Eigen::MatrixXd &HESSIAN,
      Eigen::MatrixXd &USCORE);

  // Single pair contribution to
  // 1. Log-likelihood
  // 2. Gradient
  void pair_contribution(
      const Setup &S,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const int K,
      const int L,
      const int GRADFLAG,
      double &LL,
      Eigen::VectorXd &GRADIENT);

}

#endif

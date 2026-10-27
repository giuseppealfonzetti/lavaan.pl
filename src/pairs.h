#ifndef pairs_H
#define pairs_H
#include <RcppEigen.h>

#include "model.h"
#include "bivariateNormal.h"

namespace pairs{

  // Bivariate normal distribution function of a pair, and its derivatives with
  // respect to (a, b, rho), on the grid of thresholds (a_0, ..., a_ck) x (b_0, ..., b_cl),
  // with a_0 = b_0 = -Inf and a_ck = b_cl = +Inf. Each cell probability is a signed
  // sum of four grid corners.
  struct Grid{
    int ck, cl, order;
    Eigen::MatrixXd F, Fa, Fb, Fr, Faa, Fbb, Fab, Far, Fbr, Frr;
    Grid(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TK, const int TL, const int CK, const int CL,
         const double RHO, const int ORDER);
    // Probability of cell (SK, SL), with its gradient G1 and Hessian G2 with respect to
    // (thresholds of k, thresholds of l, rho)
    double cell(const int SK, const int SL, Eigen::VectorXd &G1, Eigen::MatrixXd &G2) const;
  };

  // Pair r of the model at a Setup: implied correlation rho and its gradient grho
  // over the active parameters idx (free loadings of k and l, latent correlations,
  // log-sds, regressions, residual covariance of the pair)
  struct Pair{
    int r, k, l, ck, cl, tk, tl, mt, off, nk, lo, mb, s;
    double rho;
    std::vector<int> idx;
    Eigen::MatrixXd Dk, Dl;
    Eigen::VectorXd x, y, X, Y, grho;
    Pair(const model::Setup &ST, const int R, const bool GRAD);
  };

  // Hessian of rho over the active parameters of the pair
  Eigen::MatrixXd rho_hessian(const model::Setup &ST, const Pair &PR);

  // Pairwise log-likelihood of pair R and, with ORDER = 1, its gradient
  void pair_contribution(const model::Setup &ST, const Eigen::Ref<const Eigen::VectorXd> FREQ, const int R,
                         const int ORDER, double &LL, Eigen::VectorXd &GRADIENT);

  // Univariate log-likelihood of the observed items, with frequencies UFREQ
  // weighted by the number of missing items of each unit (available-cases term)
  void univariate_contribution(const Eigen::Ref<const Eigen::VectorXd> TAU, const model::Layout &L,
                               const Eigen::Ref<const Eigen::VectorXd> UFREQ, const int ORDER, double &LL,
                               Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &USCORE);

}

#endif

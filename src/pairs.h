#ifndef pairs_H
#define pairs_H
#include <RcppEigen.h>

#include "model.h"
#include "bivariateNormal.h"

namespace pairs{

  struct Grid{
    int ck, cl, order;
    Eigen::MatrixXd F, Fa, Fb, Fr, Faa, Fbb, Fab, Far, Fbr, Frr;
    Grid(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TK, const int TL, const int CK, const int CL,
         const double RHO, const int ORDER);
    double cell(const int SK, const int SL, Eigen::VectorXd &G1, Eigen::MatrixXd &G2) const;
  };

  struct Pair{
    int r, k, l, ck, cl, tk, tl, mt, off, nk, lo, mb, s, s0, pa, psk, psl;
    bool plain;
    double rho, c;
    std::vector<int> idx;
    Eigen::MatrixXd Dk, Dl, Ju;
    Eigen::VectorXd x, y, X, Y, grho, gc, gmk, gml, gvk, gvl, gsk, gsl;
    Pair(const model::Setup &ST, const int R, const bool GRAD, const int K = -1);
    std::vector<int> local() const;
    Eigen::VectorXd lift(const Eigen::Ref<const Eigen::VectorXd> V) const;
  };

  Eigen::MatrixXd rho_hessian(const model::Setup &ST, const Pair &PR);

  Eigen::MatrixXd curvature(const model::Setup &ST, const Pair &PR, const Eigen::Ref<const Eigen::VectorXd> GU);

  void pair_contribution(const model::Setup &ST, const Eigen::Ref<const Eigen::VectorXd> FREQ, const int R,
                         const int ORDER, double &LL, Eigen::VectorXd &GRADIENT);

  void univariate_contribution(const Eigen::Ref<const Eigen::VectorXd> TAU, const model::Layout &L,
                               const Eigen::Ref<const Eigen::VectorXd> UFREQ, const int ORDER, double &LL,
                               Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &USCORE,
                               const Eigen::Ref<const Eigen::VectorXd> ITEMW = Eigen::VectorXd());

  void univariate(const model::Setup &ST, const Eigen::Ref<const Eigen::VectorXd> UFREQ, const int ORDER, double &LL,
                  Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &USCORE,
                  const Eigen::Ref<const Eigen::VectorXd> ITEMW = Eigen::VectorXd());

}

#endif

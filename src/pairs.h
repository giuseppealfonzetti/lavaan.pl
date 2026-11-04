#ifndef pairs_H
#define pairs_H
#include <RcppEigen.h>

#include "model.h"
#include "bivariateNormal.h"

namespace pairs{

  struct Data{
    const Eigen::Map<Eigen::VectorXd> freq, ufreq, stratum, w;
    const Eigen::Map<Eigen::MatrixXd> x, y;
    const bool tables;
    Data(const Rcpp::List D);
  };

  struct Reduced{
    const Eigen::Ref<const Eigen::VectorXd> tau;
    const Eigen::Ref<const Eigen::MatrixXd> pis;
    const Eigen::Ref<const Eigen::VectorXd> loc, scale;
    Reduced(const model::Setup &ST): tau(ST.taus), pis(ST.pis), loc(ST.loc), scale(ST.s){}
    Reduced(const Eigen::Ref<const Eigen::VectorXd> TAU, const Eigen::Ref<const Eigen::MatrixXd> PIS,
            const Eigen::Ref<const Eigen::VectorXd> LOC, const Eigen::Ref<const Eigen::VectorXd> SCALE):
      tau(TAU), pis(PIS), loc(LOC), scale(SCALE){}
  };

  struct Work{
    Eigen::VectorXd g, gc, wk, wl;
    Eigen::MatrixXd M, G2;
  };

  struct Grid{
    int ck, cl, order;
    double rho, r2, sr;
    Eigen::VectorXd a, b;
    Eigen::MatrixXd F, Fa, Fb, Fr, Faa, Fbb, Fab, Far, Fbr, Frr;
    Grid(const int CK, const int CL, const double RHO, const int ORDER);
    Grid(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TK, const int TL, const int CK, const int CL,
         const double RHO, const int ORDER, const double SHK = 0, const double SHL = 0);
    void fill(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TK, const int TL, const double SHK, const double SHL,
              const int S0, const int S1, const int T0, const int T1);
    double cell(const int SK, const int SL, Eigen::VectorXd &G1, Eigen::MatrixXd &G2) const;
  };

  struct Pair{
    int r, k, l, ck, cl, tk, tl, mt, ml, lk, ll, nsk, nsl, m, off, nk, lo, mb, s, s0, pa, psk, psl, pg, pkk, pkl, nx;
    bool plain, numk, numl;
    double rho, c;
    std::vector<int> idx;
    Eigen::MatrixXd Dk, Dl, Ju, Gk, Gl;
    Eigen::VectorXd x, y, X, Y, grho, gc, gmk, gml, gvk, gvl, gsk, gsl;
    Pair(const model::Setup &ST, const int R, const bool GRAD, const int K = -1);
    Pair(const model::Layout &L, const int R, const double RHO, const int NX);
    std::vector<int> local() const;
    Eigen::VectorXd lift(const Eigen::Ref<const Eigen::VectorXd> V) const;
    void size_reduced(const int NX);
    void slots(const Eigen::Ref<const Eigen::VectorXd> X, const double YK, const double YL, Eigen::VectorXd &WK, Eigen::VectorXd &WL) const;
    void add_score(const double N, const Eigen::Ref<const Eigen::VectorXd> G, const Eigen::Ref<const Eigen::VectorXd> WK,
                   const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::VectorXd &OUT) const;
    void add_info(const double N, const Eigen::Ref<const Eigen::MatrixXd> M, const Eigen::Ref<const Eigen::VectorXd> WK,
                  const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::MatrixXd &OUT) const;
    void add_scale(const double N, const Eigen::Ref<const Eigen::VectorXd> SCALE, Eigen::VectorXd &GU, Eigen::MatrixXd *M) const;
  };

  Eigen::MatrixXd rho_hessian(const model::Setup &ST, const Pair &PR);

  Eigen::MatrixXd curvature(const model::Setup &ST, const Pair &PR, const Eigen::Ref<const Eigen::VectorXd> GU);

  double unit(const Reduced &R, const Data &D, const Pair &PR, Grid &G, const int U, const int ORDER, Work &W);

  void moments(const model::Layout &L, const Reduced &R, const Data &D, const Pair &PR, const int ORDER, double &LL,
               Eigen::VectorXd &GU, Eigen::MatrixXd &M, Eigen::MatrixXd *SC = nullptr);

  void pair_contribution(const model::Setup &ST, const Data &D, const int R, const int ORDER, double &LL,
                         Eigen::VectorXd &GRADIENT);

  void univariate_contribution(const Eigen::Ref<const Eigen::VectorXd> TAU, const model::Layout &L,
                               const Eigen::Ref<const Eigen::VectorXd> UFREQ, const int ORDER, double &LL,
                               Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &USCORE,
                               const Eigen::Ref<const Eigen::VectorXd> ITEMW = Eigen::VectorXd());

  void univariate(const model::Setup &ST, const Data &D, const int ORDER, double &LL, Eigen::VectorXd &GRADIENT,
                  Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &USCORE,
                  const Eigen::Ref<const Eigen::VectorXd> ITEMW = Eigen::VectorXd());

}

#endif

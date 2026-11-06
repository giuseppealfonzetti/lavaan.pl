#ifndef pairs_H
#define pairs_H
#include <RcppEigen.h>

#include "model.h"
#include "bivariate_normal.h"

namespace pairs{

  enum class Order {value, gradient, hessian};

  struct GroupData{
    const Eigen::Map<Eigen::VectorXd> freq, ufreq, stratum, weights;
    const Eigen::Map<Eigen::MatrixXd> strata, responses;
    const bool tables;
    GroupData(const Rcpp::List DATA);
  };

  struct Reduced{
    const Eigen::Ref<const Eigen::VectorXd> tau;
    const Eigen::Ref<const Eigen::MatrixXd> slope;
    const Eigen::Ref<const Eigen::VectorXd> location, scale;
    Reduced(const model::Matrices &PAR): tau(PAR.tau_std), slope(PAR.pi_std), location(PAR.location), scale(PAR.scale){}
    Reduced(const Eigen::Ref<const Eigen::VectorXd> TAU, const Eigen::Ref<const Eigen::MatrixXd> SLOPE,
            const Eigen::Ref<const Eigen::VectorXd> LOCATION, const Eigen::Ref<const Eigen::VectorXd> SCALE):
      tau(TAU), slope(SLOPE), location(LOCATION), scale(SCALE){}
  };

  struct ReducedIndex{
    const model::Items &items;
    const int nx, size;
    ReducedIndex(const model::Items &ITEMS, const int NX):
      items(ITEMS), nx(NX), size(ITEMS.nthr + NX * ITEMS.p + 2 * ITEMS.nnum){}
    int slope(const int C, const int K) const { return items.nthr + C * items.p + K; }
    int location(const int K) const { return items.nthr + nx * items.p + items.numeric_index(K); }
    int scale(const int K) const { return location(K) + items.nnum; }
    int rho(const int R) const { return size + R; }
  };

  struct UnitWork{
    Eigen::VectorXd score, cell_score, weights_k, weights_l;
    Eigen::MatrixXd info, cell_hessian;
  };

  struct Grid{
    int cat_k, cat_l;
    Order order;
    double rho, var_cond, sd_cond;
    Eigen::VectorXd bound_k, bound_l, density_k, density_l;
    Eigen::MatrixXd cdf, cdf_a, cdf_b, cdf_rho, cdf_aa, cdf_bb, cdf_ab, cdf_arho, cdf_brho, cdf_rhorho;
    Grid(const int CAT_K, const int CAT_L, const double RHO, const Order ORDER);
    Grid(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TAU_K, const int TAU_L, const int CAT_K, const int CAT_L,
         const double RHO, const Order ORDER, const double SHIFT_K = 0, const double SHIFT_L = 0);
    void fill(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TAU_K, const int TAU_L, const double SHIFT_K,
              const double SHIFT_L, const int S0, const int S1, const int T0, const int T1);
    double cell_probability(const int SK, const int SL, Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN) const;
  };

  struct Pair{
    int index, k, l, cat_k, cat_l, tau_offset_k, tau_offset_l, nthr, nloc_k, nloc_l, nloc, nslope_k, nslope_l, nreduced,
      cell_offset, nload_k, first_corr, first_beta, nlocal, ncov_local, first_alpha, scale_pos_k, scale_pos_l, first_gamma,
      first_kappa_k, first_kappa_l, nx;
    bool standardised, numeric_k, numeric_l;
    double rho, cov;
    std::vector<int> theta_index;
    Eigen::MatrixXd d_lambda_row_k, d_lambda_row_l, jacobian, d_pi_k, d_pi_l;
    Eigen::VectorXd row_k, row_l, chol_row_k, chol_row_l, d_cov, d_cov_full, d_mu_k, d_mu_l, d_var_k, d_var_l, d_scale_k,
      d_scale_l;
    Pair(const model::Matrices &PAR, const int INDEX, const bool DERIVATIVES, const int ITEM = -1);
    Pair(const model::Items &ITEMS, const int INDEX, const double RHO, const int NX);
    std::vector<int> parameter_indices() const;
    Eigen::VectorXd to_local(const Eigen::Ref<const Eigen::VectorXd> V) const;
    void set_reduced_dims(const int NX);
    void covariate_weights(const Eigen::Ref<const Eigen::VectorXd> X, const double YK, const double YL, Eigen::VectorXd &WK,
                           Eigen::VectorXd &WL) const;
    void add_score(const double N, const Eigen::Ref<const Eigen::VectorXd> G, const Eigen::Ref<const Eigen::VectorXd> WK,
                   const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::VectorXd &OUT) const;
    void add_info(const double N, const Eigen::Ref<const Eigen::MatrixXd> M, const Eigen::Ref<const Eigen::VectorXd> WK,
                  const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::MatrixXd &OUT) const;
    void add_log_scale(const double N, const Eigen::Ref<const Eigen::VectorXd> SCALE, Eigen::VectorXd &SCORE,
                       Eigen::MatrixXd *INFO) const;
  };

  Eigen::MatrixXd cov_hessian(const model::Matrices &PAR, const Pair &PAIR);

  Eigen::MatrixXd reduced_curvature(const model::Matrices &PAR, const Pair &PAIR, const Eigen::Ref<const Eigen::VectorXd> SCORE);

  double unit_loglik(const Reduced &REDUCED, const GroupData &DATA, const Pair &PAIR, Grid &GRID, const int UNIT, const Order ORDER,
                     UnitWork &WORK);

  void pair_loglik(const model::Items &ITEMS, const Reduced &REDUCED, const GroupData &DATA, const Pair &PAIR, const Order ORDER,
                   double &LL, Eigen::VectorXd &SCORE, Eigen::MatrixXd &INFO, Eigen::MatrixXd *CELL_SCORES = nullptr);

  void add_pair_gradient(const model::Matrices &PAR, const GroupData &DATA, const int INDEX, const Order ORDER, double &LL,
                         Eigen::VectorXd &GRADIENT);

  void univariate_margins(const Eigen::Ref<const Eigen::VectorXd> TAU, const model::Items &ITEMS,
                          const Eigen::Ref<const Eigen::VectorXd> UFREQ, const Order ORDER, double &LL,
                          Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &SCORES,
                          const Eigen::Ref<const Eigen::VectorXd> ITEM_WEIGHTS = Eigen::VectorXd());

  void add_univariate(const model::Matrices &PAR, const GroupData &DATA, const Order ORDER, double &LL, Eigen::VectorXd &GRADIENT,
                      Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &SCORES,
                      const Eigen::Ref<const Eigen::VectorXd> ITEM_WEIGHTS = Eigen::VectorXd());

}

#endif

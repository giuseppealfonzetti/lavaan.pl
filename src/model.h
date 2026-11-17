#ifndef model_H
#define model_H
#include <RcppEigen.h>

namespace model{

  enum class LatCov {none, lkj, direct};

  struct Items{
    int p, npairs, nthr, ncell, nnum;
    Eigen::VectorXi cat, tau_offset, numeric_index;
    std::vector<int> pair_k, pair_l, cell_offset;
    Items(const Eigen::Ref<const Eigen::VectorXd> CAT);
  };

  struct Blocks{
    int lambda, corr, logsd, beta, rescov, nu, alpha, scale, gamma, kappa;
  };

  struct Structure{
    Items items;
    int p, q, nx, d, d_user, nimplied, nload, ncorr, nvar, nreg, nres, nint, nmean, nscale, ngam, nkap;
    Blocks first;
    LatCov latcov;
    bool theta_parameterization, standardised;
    std::vector<bool> scale_from_rsd;
    Eigen::MatrixXd lambda, beta, theta, gamma, kappa;
    Eigen::VectorXd logsd, nu, alpha, logscale;
    std::vector<std::vector<int>> lambda_index, lambda_factor, kappa_index, kappa_covariate;
    std::vector<int> corr_row, corr_col, logsd_factor, beta_row, beta_col, nu_index, alpha_factor, scale_index,
      gamma_factor, gamma_covariate, implied_item, implied_slot, user_index, ext_to_user;
    Eigen::MatrixXi rescov_index;
    Structure(const Rcpp::List MODEL, const Eigen::Ref<const Eigen::VectorXd> CAT);
  };

  struct Matrices{
    const Structure &structure;
    Eigen::VectorXd tau, sd, alpha, mu, common_var, residual_var, scale, tau_std, location;
    Eigen::VectorXd user_vector;
    Eigen::MatrixXd lambda, lkj_z, lkj_w, lkj_cumw, chol, dchol, psi, ib_inv, lambda_ib_inv, theta, gamma, pi, pi_std,
      implied_jacobian;
    std::vector<Eigen::MatrixXd> lambda_jacobian;
    Matrices(const Eigen::Ref<const Eigen::VectorXd> PARAMETER_VECTOR, const Structure &STRUCTURE);
    double rescov(const int K, const int L) const { return theta.size() > 0 ? theta(K, L) : 0; }
    Eigen::MatrixXd to_user(const Eigen::Ref<const Eigen::MatrixXd> M) const;
    Eigen::MatrixXd hessian_to_user(const Eigen::Ref<const Eigen::MatrixXd> H, const Eigen::Ref<const Eigen::VectorXd> G) const;
  };

  struct Groups{
    std::vector<Structure> structure;
    std::vector<int> parameter_offset, pair_offset;
    int ngroups, d, npairs;
    double n;
    Groups(const Rcpp::List MODEL, const Rcpp::List DATA);
    std::vector<Matrices> matrices(const Eigen::Ref<const Eigen::VectorXd> PARAMETER_VECTOR) const;
  };

}

#endif

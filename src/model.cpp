#include "model.h"

model::Items::Items(const Eigen::Ref<const Eigen::VectorXd> CAT): p(CAT.size()), npairs(0), nthr(0), ncell(0), nnum(0){
  cat = CAT.cast<int>();
  tau_offset.resize(p);
  numeric_index = Eigen::VectorXi::Constant(p, -1);
  for(int k = 0; k < p; k++){
    tau_offset(k) = nthr;
    nthr += cat(k) - 1;
    if(cat(k) == 1) numeric_index(k) = nnum++;
  }
  for(int k = 1; k < p; k++){
    for(int l = 0; l < k; l++){
      pair_k.push_back(k);
      pair_l.push_back(l);
      cell_offset.push_back(ncell);
      if(cat(k) > 1 && cat(l) > 1) ncell += cat(k) * cat(l);
      npairs++;
    }
  }
}

model::Structure::Structure(const Rcpp::List MODEL, const Eigen::Ref<const Eigen::VectorXd> CAT): items(CAT){
  lambda = Rcpp::as<Eigen::MatrixXd>(MODEL["LAMBDA"]);
  logsd = Rcpp::as<Eigen::VectorXd>(MODEL["LOGSD"]);
  beta = Rcpp::as<Eigen::MatrixXd>(MODEL["BETA"]);
  theta = Rcpp::as<Eigen::MatrixXd>(MODEL["THETA"]);
  nu = Rcpp::as<Eigen::VectorXd>(MODEL["NU"]);
  alpha = Rcpp::as<Eigen::VectorXd>(MODEL["ALPHA"]);
  logscale = Rcpp::as<Eigen::VectorXd>(MODEL["LOGSCALE"]);
  gamma = Rcpp::as<Eigen::MatrixXd>(MODEL["GAMMA"]);
  kappa = Rcpp::as<Eigen::MatrixXd>(MODEL["KAPPA"]);
  theta_parameterization = Rcpp::as<std::string>(MODEL["PARAMETERIZATION"]) == "theta";
  const std::string lc = Rcpp::as<std::string>(MODEL["LATCOV"]);
  latcov = lc == "direct" ? LatCov::direct : lc == "lkj" ? LatCov::lkj : LatCov::none;
  p = lambda.rows();
  q = lambda.cols();

  first.lambda = items.nthr;
  lambda_index.assign(p, std::vector<int>());
  lambda_factor.assign(p, std::vector<int>());
  nload = 0;
  for(int h = 0; h < q; h++){
    for(int j = 0; j < p; j++){
      if(!std::isfinite(lambda(j, h))){
        lambda_index[j].push_back(first.lambda + nload);
        lambda_factor[j].push_back(h);
        nload++;
      }
    }
  }

  first.corr = first.lambda + nload;
  for(int c = 1; c < q && latcov != LatCov::none; c++){
    for(int r = 0; r < c; r++){
      corr_row.push_back(r);
      corr_col.push_back(c);
    }
  }
  ncorr = corr_row.size();

  first.logsd = first.corr + ncorr;
  for(int h = 0; h < q; h++) if(!std::isfinite(logsd(h))) logsd_factor.push_back(h);
  nvar = logsd_factor.size();

  first.beta = first.logsd + nvar;
  for(int j = 0; j < beta.cols(); j++){
    for(int i = 0; i < beta.rows(); i++){
      if(!std::isfinite(beta(i, j))){
        beta_row.push_back(i);
        beta_col.push_back(j);
      }
    }
  }
  nreg = beta_row.size();

  first.rescov = first.beta + nreg;
  nres = 0;
  rescov_index = Eigen::MatrixXi::Constant(p, p, -1);
  for(int c = 1; c < theta.cols(); c++){
    for(int r = 0; r < c; r++){
      if(!std::isfinite(theta(r, c))){
        rescov_index(r, c) = rescov_index(c, r) = first.rescov + nres;
        nres++;
      }
    }
  }
  first.nu = first.rescov + nres;
  nu_index.assign(p, -1);
  nint = 0;
  for(int j = 0; j < p; j++) if(!std::isfinite(nu(j))) nu_index[j] = first.nu + nint++;
  first.alpha = first.nu + nint;
  for(int h = 0; h < q; h++) if(!std::isfinite(alpha(h))) alpha_factor.push_back(h);
  nmean = alpha_factor.size();
  first.scale = first.alpha + nmean;
  scale_index.assign(p, -1);
  nscale = 0;
  for(int j = 0; j < p; j++) if(!std::isfinite(logscale(j))) scale_index[j] = first.scale + nscale++;
  first.gamma = first.scale + nscale;
  nx = gamma.cols();
  for(int c = 0; c < nx; c++){
    for(int h = 0; h < q; h++){
      if(!std::isfinite(gamma(h, c))){
        gamma_factor.push_back(h);
        gamma_covariate.push_back(c);
      }
    }
  }
  ngam = gamma_factor.size();
  first.kappa = first.gamma + ngam;
  kappa_index.assign(p, std::vector<int>());
  kappa_covariate.assign(p, std::vector<int>());
  nkap = 0;
  for(int c = 0; c < nx; c++){
    for(int j = 0; j < p; j++){
      if(!std::isfinite(kappa(j, c))){
        kappa_index[j].push_back(first.kappa + nkap++);
        kappa_covariate[j].push_back(c);
      }
    }
  }
  d = first.kappa + nkap;
  scale_from_rsd.assign(p, theta_parameterization);
  for(int j = 0; j < p; j++) if(items.cat(j) == 1) scale_from_rsd[j] = true;
  standardised = !theta_parameterization && items.nnum == 0 && nx == 0 && nint + nmean + nscale == 0 &&
    (nu.array() == 0).all() && (alpha.array() == 0).all() && (logscale.array() == 0).all();
}

model::Matrices::Matrices(const Eigen::Ref<const Eigen::VectorXd> PARAMETER_VECTOR, const Structure &STRUCTURE): structure(STRUCTURE){
  const Structure &S = structure;
  const Eigen::Ref<const Eigen::VectorXd> x = PARAMETER_VECTOR;
  const int p = S.p, q = S.q;
  if(x.size() != S.d) Rcpp::stop("PARAMETER_VECTOR has length %d, the model has %d parameters.", (int)x.size(), S.d);

  tau = x.head(S.items.nthr);
  lambda = S.lambda;
  for(int j = 0; j < p; j++){
    for(std::size_t a = 0; a < S.lambda_index[j].size(); a++) lambda(j, S.lambda_factor[j][a]) = x(S.lambda_index[j][a]);
  }
  sd = S.logsd;
  for(int s = 0; s < S.nvar; s++) sd(S.logsd_factor[s]) = x(S.first.logsd + s);
  sd = sd.array().exp();

  lkj_z = Eigen::MatrixXd::Zero(q, q);
  if(S.latcov != LatCov::direct) for(int a = 0; a < S.ncorr; a++) lkj_z(S.corr_row[a], S.corr_col[a]) = std::tanh(x(S.first.corr + a));
  lkj_w = (1 - lkj_z.array().square()).sqrt().matrix();
  lkj_cumw = Eigen::MatrixXd::Ones(q, q);
  chol = Eigen::MatrixXd::Zero(q, q);
  for(int c = 0; c < q; c++){
    for(int r = 1; r <= c; r++) lkj_cumw(r, c) = lkj_cumw(r - 1, c) * lkj_w(r - 1, c);
    for(int r = 0; r < c; r++) chol(r, c) = lkj_z(r, c) * lkj_cumw(r, c);
    chol(c, c) = lkj_cumw(c, c);
  }
  dchol = Eigen::MatrixXd::Zero(q, S.ncorr);
  for(int a = 0; a < S.ncorr && S.latcov != LatCov::direct; a++){
    const int r = S.corr_row[a], c = S.corr_col[a];
    dchol(r, a) = lkj_cumw(r, c) * lkj_w(r, c) * lkj_w(r, c);
    for(int i = r + 1; i <= c; i++) dchol(i, a) = -chol(i, c) * lkj_z(r, c);
  }
  if(S.latcov == LatCov::direct){
    psi = sd.array().square().matrix().asDiagonal();
    for(int a = 0; a < S.ncorr; a++) psi(S.corr_row[a], S.corr_col[a]) = psi(S.corr_col[a], S.corr_row[a]) = x(S.first.corr + a);
  }else{
    psi = sd.asDiagonal() * chol.transpose() * chol * sd.asDiagonal();
  }

  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(q, q);
  if(S.beta.size() > 0) B = S.beta;
  for(int b = 0; b < S.nreg; b++) B(S.beta_row[b], S.beta_col[b]) = x(S.first.beta + b);
  ib_inv = (Eigen::MatrixXd::Identity(q, q) - B).inverse();
  lambda_ib_inv = lambda * ib_inv;
  common_var = (lambda_ib_inv * psi).cwiseProduct(lambda_ib_inv).rowwise().sum();

  if(S.theta.size() > 0){
    theta = S.theta;
    theta.diagonal().setZero();
    for(int c = 1; c < p; c++){
      for(int r = 0; r < c; r++) if(S.rescov_index(r, c) >= 0) theta(r, c) = theta(c, r) = x(S.rescov_index(r, c));
    }
  }

  tau_std = tau;
  pi_std.resize(p, S.nx);
  if(!S.standardised){
    Eigen::VectorXd nu = S.nu, logscale = S.logscale;
    alpha = S.alpha;
    for(int j = 0; j < p; j++){
      if(S.nu_index[j] >= 0) nu(j) = x(S.nu_index[j]);
      if(S.scale_index[j] >= 0) logscale(j) = x(S.scale_index[j]);
    }
    for(int a = 0; a < S.nmean; a++) alpha(S.alpha_factor[a]) = x(S.first.alpha + a);
    mu = nu + lambda_ib_inv * alpha;
    residual_var = (2 * logscale).array().exp();
    scale = logscale.array().exp();
    for(int j = 0; j < p; j++) if(S.scale_from_rsd[j]) scale(j) = 1 / std::sqrt(common_var(j) + residual_var(j));
    location = -scale.cwiseProduct(mu);
    for(int j = 0; j < p; j++){
      for(int c = 0; c < S.items.cat(j) - 1; c++){
        const int t = S.items.tau_offset(j) + c;
        tau_std(t) = scale(j) * (tau(t) - mu(j));
      }
    }
    gamma = S.gamma;
    Eigen::MatrixXd kappa = S.kappa;
    for(int a = 0; a < S.ngam; a++) gamma(S.gamma_factor[a], S.gamma_covariate[a]) = x(S.first.gamma + a);
    for(int j = 0; j < p; j++){
      for(std::size_t a = 0; a < S.kappa_index[j].size(); a++) kappa(j, S.kappa_covariate[j][a]) = x(S.kappa_index[j][a]);
    }
    pi = lambda_ib_inv * gamma + kappa;
    pi_std = scale.asDiagonal() * pi;
  }

  lambda_jacobian.resize(p);
  for(int j = 0; j < p; j++){
    const int nl = S.lambda_index[j].size();
    lambda_jacobian[j].resize(q, nl + S.nreg);
    for(int a = 0; a < nl; a++) lambda_jacobian[j].col(a) = ib_inv.row(S.lambda_factor[j][a]).transpose();
    for(int b = 0; b < S.nreg; b++) lambda_jacobian[j].col(nl + b) = lambda_ib_inv(j, S.beta_row[b]) * ib_inv.row(S.beta_col[b]).transpose();
  }
}

model::Groups::Groups(const Rcpp::List MODEL, const Rcpp::List DATA): ngroups(MODEL.size()), d(0), npairs(0), n(0){
  structure.reserve(ngroups);
  for(int g = 0; g < ngroups; g++){
    const Rcpp::List data = DATA[g];
    structure.emplace_back(Rcpp::as<Rcpp::List>(MODEL[g]), Rcpp::as<Eigen::VectorXd>(data["CAT"]));
    parameter_offset.push_back(d);
    pair_offset.push_back(npairs);
    d += structure[g].d;
    npairs += structure[g].items.npairs;
    n += Rcpp::as<double>(data["N"]);
  }
}

std::vector<model::Matrices> model::Groups::matrices(const Eigen::Ref<const Eigen::VectorXd> PARAMETER_VECTOR) const{
  if(PARAMETER_VECTOR.size() != d) Rcpp::stop("PARAMETER_VECTOR has length %d, the model has %d parameters.", (int)PARAMETER_VECTOR.size(), d);
  std::vector<Matrices> out;
  out.reserve(ngroups);
  for(int g = 0; g < ngroups; g++) out.emplace_back(PARAMETER_VECTOR.segment(parameter_offset[g], structure[g].d), structure[g]);
  return out;
}

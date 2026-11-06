#include "pairs.h"

pairs::GroupData::GroupData(const Rcpp::List DATA):
  freq(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(DATA["FREQ"])), ufreq(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(DATA["UFREQ"])),
  stratum(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(DATA["STRATUM"])), weights(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(DATA["WEIGHTS"])),
  strata(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(DATA["STRATA"])), responses(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(DATA["Y"])),
  tables(Rcpp::as<bool>(DATA["TABLES"])){}

pairs::Grid::Grid(const int CAT_K, const int CAT_L, const double RHO, const Order ORDER):
  cat_k(CAT_K), cat_l(CAT_L), order(ORDER), rho(RHO), var_cond(1 - RHO * RHO), sd_cond(std::sqrt(1 - RHO * RHO)){
  const double inf = std::numeric_limits<double>::infinity();
  bound_k.resize(cat_k + 1);
  bound_l.resize(cat_l + 1);
  bound_k(0) = bound_l(0) = -inf;
  bound_k(cat_k) = bound_l(cat_l) = inf;
  cdf = Eigen::MatrixXd::Zero(cat_k + 1, cat_l + 1);
  if(order == Order::value) return;
  density_k.resize(cat_k + 1);
  density_l.resize(cat_l + 1);
  cdf_a = cdf_b = cdf_rho = Eigen::MatrixXd::Zero(cat_k + 1, cat_l + 1);
  if(order == Order::hessian) cdf_aa = cdf_bb = cdf_ab = cdf_arho = cdf_brho = cdf_rhorho = Eigen::MatrixXd::Zero(cat_k + 1, cat_l + 1);
}

pairs::Grid::Grid(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TAU_K, const int TAU_L, const int CAT_K, const int CAT_L,
                  const double RHO, const Order ORDER, const double SHIFT_K, const double SHIFT_L): Grid(CAT_K, CAT_L, RHO, ORDER){
  fill(TAU, TAU_K, TAU_L, SHIFT_K, SHIFT_L, 1, cat_k, 1, cat_l);
}

void pairs::Grid::fill(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TAU_K, const int TAU_L, const double SHIFT_K,
                       const double SHIFT_L, const int S0, const int S1, const int T0, const int T1){
  const double inf = std::numeric_limits<double>::infinity();
  bound_k.segment(1, cat_k - 1) = TAU.segment(TAU_K, cat_k - 1).array() - SHIFT_K;
  bound_l.segment(1, cat_l - 1) = TAU.segment(TAU_L, cat_l - 1).array() - SHIFT_L;
  for(int s = S0; s <= S1; s++){
    for(int t = T0; t <= T1; t++){
      if(s == cat_k && t == cat_l) cdf(s, t) = 1;
      else if(s == cat_k) cdf(s, t) = R::pnorm(bound_l(t), 0, 1, 1, 0);
      else if(t == cat_l) cdf(s, t) = R::pnorm(bound_k(s), 0, 1, 1, 0);
      else cdf(s, t) = binorm::pbvnorm(bound_k(s), bound_l(t), rho);
    }
  }
  if(order == Order::value) return;

  for(int s = S0; s <= S1 && s < cat_k; s++) density_k(s) = R::dnorm(bound_k(s), 0, 1, 0);
  for(int t = T0; t <= T1 && t < cat_l; t++) density_l(t) = R::dnorm(bound_l(t), 0, 1, 0);
  for(int s = S0; s <= S1; s++){
    for(int t = T0; t <= T1; t++){
      const bool fa = s < cat_k, fb = t < cat_l;
      const double a = bound_k(s), b = bound_l(t);
      if(fa){
        const double pa = density_k(s);
        const double w = fb ? (b - rho * a) / sd_cond : inf;
        cdf_a(s, t) = pa * R::pnorm(w, 0, 1, 1, 0);
        if(order == Order::hessian) cdf_aa(s, t) = -a * cdf_a(s, t) - (fb ? rho / sd_cond * pa * R::dnorm(w, 0, 1, 0) : 0);
      }
      if(fb){
        const double pb = density_l(t);
        const double w = fa ? (a - rho * b) / sd_cond : inf;
        cdf_b(s, t) = pb * R::pnorm(w, 0, 1, 1, 0);
        if(order == Order::hessian) cdf_bb(s, t) = -b * cdf_b(s, t) - (fa ? rho / sd_cond * pb * R::dnorm(w, 0, 1, 0) : 0);
      }
      if(fa && fb){
        const double f2 = binorm::dbvnorm(a, b, rho);
        cdf_rho(s, t) = f2;
        if(order == Order::hessian){
          cdf_ab(s, t) = f2;
          cdf_arho(s, t) = -f2 * (a - rho * b) / var_cond;
          cdf_brho(s, t) = -f2 * (b - rho * a) / var_cond;
          cdf_rhorho(s, t) = f2 * (rho / var_cond + a * b / var_cond -
                                   rho * (a * a - 2 * rho * a * b + b * b) / (var_cond * var_cond));
        }
      }
    }
  }
}

double pairs::Grid::cell_probability(const int SK, const int SL, Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN) const{
  const double pe = cdf(SK + 1, SL + 1) - cdf(SK + 1, SL) - cdf(SK, SL + 1) + cdf(SK, SL) + 1e-8;
  if(order == Order::value) return pe;

  const int rho_pos = cat_k + cat_l - 2;
  GRADIENT = Eigen::VectorXd::Zero(rho_pos + 1);
  if(order == Order::hessian) HESSIAN = Eigen::MatrixXd::Zero(rho_pos + 1, rho_pos + 1);
  for(int e = 0; e < 2; e++){
    for(int f = 0; f < 2; f++){
      const int s = SK + 1 - e, t = SL + 1 - f;
      const double sign = (e == f) ? 1 : -1;
      const int pos_k = (s > 0 && s < cat_k) ? s - 1 : -1;
      const int pos_l = (t > 0 && t < cat_l) ? cat_k - 2 + t : -1;
      if(pos_k >= 0) GRADIENT(pos_k) += sign * cdf_a(s, t);
      if(pos_l >= 0) GRADIENT(pos_l) += sign * cdf_b(s, t);
      GRADIENT(rho_pos) += sign * cdf_rho(s, t);
      if(order != Order::hessian) continue;
      if(pos_k >= 0){
        HESSIAN(pos_k, pos_k) += sign * cdf_aa(s, t);
        HESSIAN(pos_k, rho_pos) += sign * cdf_arho(s, t);
        HESSIAN(rho_pos, pos_k) += sign * cdf_arho(s, t);
      }
      if(pos_l >= 0){
        HESSIAN(pos_l, pos_l) += sign * cdf_bb(s, t);
        HESSIAN(pos_l, rho_pos) += sign * cdf_brho(s, t);
        HESSIAN(rho_pos, pos_l) += sign * cdf_brho(s, t);
      }
      if(pos_k >= 0 && pos_l >= 0){
        HESSIAN(pos_k, pos_l) += sign * cdf_ab(s, t);
        HESSIAN(pos_l, pos_k) += sign * cdf_ab(s, t);
      }
      HESSIAN(rho_pos, rho_pos) += sign * cdf_rhorho(s, t);
    }
  }
  return pe;
}

static Eigen::VectorXd row_product_gradient(const model::Matrices &PAR, const int FIRST_CORR, const Eigen::MatrixXd &D_A,
                                            const Eigen::MatrixXd &D_B, const Eigen::VectorXd &ROW_A, const Eigen::VectorXd &ROW_B,
                                            const Eigen::VectorXd &CHOL_A, const Eigen::VectorXd &CHOL_B){
  const model::Structure &S = PAR.structure;
  Eigen::VectorXd g;
  if(S.latcov == model::LatCov::direct){
    g = D_A.transpose() * (PAR.psi * ROW_B) + D_B.transpose() * (PAR.psi * ROW_A);
    for(int a = 0; a < S.ncorr; a++) g(FIRST_CORR + a) = ROW_A(S.corr_row[a]) * ROW_B(S.corr_col[a]) + ROW_A(S.corr_col[a]) * ROW_B(S.corr_row[a]);
    for(int v = 0; v < S.nvar; v++){
      const int h = S.logsd_factor[v];
      g(FIRST_CORR + S.ncorr + v) = 2 * PAR.sd(h) * PAR.sd(h) * ROW_A(h) * ROW_B(h);
    }
  }else{
    g = D_A.transpose() * PAR.sd.cwiseProduct(PAR.chol.transpose() * CHOL_B) +
      D_B.transpose() * PAR.sd.cwiseProduct(PAR.chol.transpose() * CHOL_A);
    for(int a = 0; a < S.ncorr; a++){
      g(FIRST_CORR + a) = ROW_A(S.corr_col[a]) * PAR.dchol.col(a).dot(CHOL_B) + ROW_B(S.corr_col[a]) * PAR.dchol.col(a).dot(CHOL_A);
    }
    for(int v = 0; v < S.nvar; v++){
      const int h = S.logsd_factor[v];
      g(FIRST_CORR + S.ncorr + v) = ROW_A(h) * PAR.chol.col(h).dot(CHOL_B) + ROW_B(h) * PAR.chol.col(h).dot(CHOL_A);
    }
  }
  return g;
}

void pairs::Pair::set_reduced_dims(const int NX){
  nx = NX;
  nthr = cat_k + cat_l - 2;
  numeric_k = cat_k == 1;
  numeric_l = cat_l == 1;
  nloc_k = numeric_k ? 1 : cat_k - 1;
  nloc_l = numeric_l ? 1 : cat_l - 1;
  nloc = nloc_k + nloc_l;
  nslope_k = nx + numeric_k;
  nslope_l = nx + numeric_l;
  nreduced = nloc + 1 + nslope_k + nslope_l;
}

pairs::Pair::Pair(const model::Items &ITEMS, const int INDEX, const double RHO, const int NX):
  index(INDEX), k(ITEMS.pair_k[INDEX]), l(ITEMS.pair_l[INDEX]), cat_k(ITEMS.cat(k)), cat_l(ITEMS.cat(l)),
  tau_offset_k(ITEMS.tau_offset(k)), tau_offset_l(ITEMS.tau_offset(l)), cell_offset(ITEMS.cell_offset[INDEX]), nlocal(0),
  standardised(false), rho(RHO){
  set_reduced_dims(NX);
}

pairs::Pair::Pair(const model::Matrices &PAR, const int INDEX, const bool DERIVATIVES, const int ITEM): index(INDEX){
  const model::Structure &S = PAR.structure;
  const model::Items &items = S.items;
  standardised = S.standardised;
  k = ITEM >= 0 ? ITEM : items.pair_k[index];
  l = ITEM >= 0 ? ITEM : items.pair_l[index];
  cat_k = items.cat(k);
  cat_l = items.cat(l);
  tau_offset_k = items.tau_offset(k);
  tau_offset_l = items.tau_offset(l);
  set_reduced_dims(S.nx);
  cell_offset = ITEM >= 0 ? 0 : items.cell_offset[index];
  const int rescov_index = ITEM >= 0 ? -1 : S.rescov_index(k, l);

  if(S.latcov == model::LatCov::direct){
    row_k = PAR.lambda_ib_inv.row(k).transpose();
    row_l = PAR.lambda_ib_inv.row(l).transpose();
    rho = row_k.dot(PAR.psi * row_l) + PAR.rescov(k, l);
  }else{
    row_k = PAR.sd.cwiseProduct(PAR.lambda_ib_inv.row(k).transpose());
    row_l = PAR.sd.cwiseProduct(PAR.lambda_ib_inv.row(l).transpose());
    chol_row_k = PAR.chol * row_k;
    chol_row_l = PAR.chol * row_l;
    rho = chol_row_k.dot(chol_row_l) + PAR.rescov(k, l);
  }
  cov = rho;
  if(!standardised) rho *= PAR.scale(k) * PAR.scale(l);
  if(!DERIVATIVES) return;

  nload_k = S.lambda_index[k].size();
  const int nload_l = S.lambda_index[l].size();
  first_corr = nload_k + nload_l;
  first_beta = first_corr + S.ncorr + S.nvar;
  nlocal = ncov_local = first_beta + S.nreg + (rescov_index >= 0);
  theta_index = S.lambda_index[k];
  theta_index.insert(theta_index.end(), S.lambda_index[l].begin(), S.lambda_index[l].end());
  for(int a = S.first.corr; a < S.first.beta + S.nreg; a++) theta_index.push_back(a);
  if(rescov_index >= 0) theta_index.push_back(rescov_index);

  d_lambda_row_k = Eigen::MatrixXd::Zero(S.q, nlocal);
  d_lambda_row_l = Eigen::MatrixXd::Zero(S.q, nlocal);
  d_lambda_row_k.leftCols(nload_k) = PAR.lambda_jacobian[k].leftCols(nload_k);
  d_lambda_row_l.middleCols(nload_k, nload_l) = PAR.lambda_jacobian[l].leftCols(nload_l);
  d_lambda_row_k.middleCols(first_beta, S.nreg) = PAR.lambda_jacobian[k].rightCols(S.nreg);
  d_lambda_row_l.middleCols(first_beta, S.nreg) = PAR.lambda_jacobian[l].rightCols(S.nreg);
  d_cov = row_product_gradient(PAR, first_corr, d_lambda_row_k, d_lambda_row_l, row_k, row_l, chol_row_k, chol_row_l);
  if(rescov_index >= 0) d_cov(nlocal - 1) = 1;
  if(standardised) return;

  auto add = [&](const int i){
    if(i < 0) return -1;
    theta_index.push_back(i);
    return nlocal++;
  };
  const int nu_pos_k = add(S.nu_index[k]), nu_pos_l = add(S.nu_index[l]);
  first_alpha = nlocal;
  for(int a = 0; a < S.nmean; a++) add(S.first.alpha + a);
  scale_pos_k = add(S.scale_index[k]);
  scale_pos_l = add(S.scale_index[l]);
  first_gamma = nlocal;
  for(int a = 0; a < S.ngam; a++) add(S.first.gamma + a);
  first_kappa_k = nlocal;
  for(const int i : S.kappa_index[k]) add(i);
  first_kappa_l = nlocal;
  for(const int i : S.kappa_index[l]) add(i);

  auto pad = [&](const Eigen::VectorXd &V){
    Eigen::VectorXd o = Eigen::VectorXd::Zero(nlocal);
    o.head(V.size()) = V;
    return o;
  };
  d_cov_full = pad(d_cov);
  d_mu_k = pad(d_lambda_row_k.transpose() * PAR.alpha);
  d_mu_l = pad(d_lambda_row_l.transpose() * PAR.alpha);
  if(nu_pos_k >= 0) d_mu_k(nu_pos_k) = 1;
  if(nu_pos_l >= 0) d_mu_l(nu_pos_l) = 1;
  for(int a = 0; a < S.nmean; a++){
    d_mu_k(first_alpha + a) = PAR.lambda_ib_inv(k, S.alpha_factor[a]);
    d_mu_l(first_alpha + a) = PAR.lambda_ib_inv(l, S.alpha_factor[a]);
  }
  d_scale_k = d_scale_l = Eigen::VectorXd::Zero(nlocal);
  if(S.scale_from_rsd[k]){
    d_var_k = pad(row_product_gradient(PAR, first_corr, d_lambda_row_k, d_lambda_row_k, row_k, row_k, chol_row_k, chol_row_k));
    d_scale_k = d_var_k;
    if(scale_pos_k >= 0) d_scale_k(scale_pos_k) += 2 * PAR.residual_var(k);
    d_scale_k *= -0.5 * std::pow(PAR.scale(k), 3);
  }else if(scale_pos_k >= 0){
    d_scale_k(scale_pos_k) = PAR.scale(k);
  }
  if(S.scale_from_rsd[l]){
    d_var_l = pad(row_product_gradient(PAR, first_corr, d_lambda_row_l, d_lambda_row_l, row_l, row_l, chol_row_l, chol_row_l));
    d_scale_l = d_var_l;
    if(scale_pos_l >= 0) d_scale_l(scale_pos_l) += 2 * PAR.residual_var(l);
    d_scale_l *= -0.5 * std::pow(PAR.scale(l), 3);
  }else if(scale_pos_l >= 0){
    d_scale_l(scale_pos_l) = PAR.scale(l);
  }

  d_pi_k = d_pi_l = Eigen::MatrixXd::Zero(nlocal, nx);
  d_pi_k.topRows(ncov_local) = d_lambda_row_k.transpose() * PAR.gamma;
  d_pi_l.topRows(ncov_local) = d_lambda_row_l.transpose() * PAR.gamma;
  for(int a = 0; a < S.ngam; a++){
    d_pi_k(first_gamma + a, S.gamma_covariate[a]) = PAR.lambda_ib_inv(k, S.gamma_factor[a]);
    d_pi_l(first_gamma + a, S.gamma_covariate[a]) = PAR.lambda_ib_inv(l, S.gamma_factor[a]);
  }
  for(std::size_t a = 0; a < S.kappa_index[k].size(); a++) d_pi_k(first_kappa_k + a, S.kappa_covariate[k][a]) = 1;
  for(std::size_t a = 0; a < S.kappa_index[l].size(); a++) d_pi_l(first_kappa_l + a, S.kappa_covariate[l][a]) = 1;

  const double sk = PAR.scale(k), sl = PAR.scale(l);
  jacobian = Eigen::MatrixXd::Zero(nreduced, nthr + nlocal);
  for(int a = 0; a < cat_k - 1; a++){
    jacobian(a, a) = sk;
    jacobian.row(a).tail(nlocal) = (-sk * d_mu_k + (PAR.tau(tau_offset_k + a) - PAR.mu(k)) * d_scale_k).transpose();
  }
  if(numeric_k) jacobian.row(0).tail(nlocal) = (-sk * d_mu_k - PAR.mu(k) * d_scale_k).transpose();
  for(int b = 0; b < cat_l - 1; b++){
    jacobian(nloc_k + b, cat_k - 1 + b) = sl;
    jacobian.row(nloc_k + b).tail(nlocal) = (-sl * d_mu_l + (PAR.tau(tau_offset_l + b) - PAR.mu(l)) * d_scale_l).transpose();
  }
  if(numeric_l) jacobian.row(nloc_k).tail(nlocal) = (-sl * d_mu_l - PAR.mu(l) * d_scale_l).transpose();
  jacobian.row(nloc).tail(nlocal) = (sl * cov * d_scale_k + sk * cov * d_scale_l + sk * sl * d_cov_full).transpose();
  const int slope_row_k = nloc + 1, slope_row_l = slope_row_k + nslope_k;
  for(int a = 0; a < nx; a++){
    jacobian.row(slope_row_k + a).tail(nlocal) = (sk * d_pi_k.col(a) + PAR.pi(k, a) * d_scale_k).transpose();
    jacobian.row(slope_row_l + a).tail(nlocal) = (sl * d_pi_l.col(a) + PAR.pi(l, a) * d_scale_l).transpose();
  }
  if(numeric_k) jacobian.row(slope_row_k + nx).tail(nlocal) = d_scale_k.transpose();
  if(numeric_l) jacobian.row(slope_row_l + nx).tail(nlocal) = d_scale_l.transpose();
}

std::vector<int> pairs::Pair::parameter_indices() const{
  std::vector<int> out;
  for(int a = 0; a < cat_k - 1; a++) out.push_back(tau_offset_k + a);
  for(int a = 0; a < cat_l - 1; a++) out.push_back(tau_offset_l + a);
  out.insert(out.end(), theta_index.begin(), theta_index.end());
  return out;
}

Eigen::VectorXd pairs::Pair::to_local(const Eigen::Ref<const Eigen::VectorXd> V) const{
  if(!standardised) return jacobian.transpose() * V;
  Eigen::VectorXd out(nthr + nlocal);
  out.head(nthr) = V.head(nthr);
  out.tail(nlocal) = V(nthr) * d_cov;
  return out;
}

void pairs::Pair::covariate_weights(const Eigen::Ref<const Eigen::VectorXd> X, const double YK, const double YL,
                                    Eigen::VectorXd &WK, Eigen::VectorXd &WL) const{
  WK.resize(nslope_k);
  WL.resize(nslope_l);
  WK.head(nx) = -X;
  WL.head(nx) = -X;
  if(numeric_k) WK(nx) = YK;
  if(numeric_l) WL(nx) = YL;
}

void pairs::Pair::add_score(const double N, const Eigen::Ref<const Eigen::VectorXd> G, const Eigen::Ref<const Eigen::VectorXd> WK,
                            const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::VectorXd &OUT) const{
  const int first_slope = nloc + 1;
  OUT.head(first_slope) += N * G;
  if(nslope_k > 0) OUT.segment(first_slope, nslope_k) += N * G.head(nloc_k).sum() * WK;
  if(nslope_l > 0) OUT.segment(first_slope + nslope_k, nslope_l) += N * G.segment(nloc_k, nloc_l).sum() * WL;
}

void pairs::Pair::add_info(const double N, const Eigen::Ref<const Eigen::MatrixXd> M, const Eigen::Ref<const Eigen::VectorXd> WK,
                           const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::MatrixXd &OUT) const{
  const int first_slope = nloc + 1;
  OUT.topLeftCorner(first_slope, first_slope) += N * M;
  if(nslope_k + nslope_l == 0) return;
  const Eigen::VectorXd u = M.leftCols(nloc_k).rowwise().sum(), v = M.middleCols(nloc_k, nloc_l).rowwise().sum();
  const double uu = u.head(nloc_k).sum(), uv = u.segment(nloc_k, nloc_l).sum(), vv = v.segment(nloc_k, nloc_l).sum();
  OUT.block(0, first_slope, first_slope, nslope_k) += N * u * WK.transpose();
  OUT.block(first_slope, 0, nslope_k, first_slope) += N * WK * u.transpose();
  OUT.block(0, first_slope + nslope_k, first_slope, nslope_l) += N * v * WL.transpose();
  OUT.block(first_slope + nslope_k, 0, nslope_l, first_slope) += N * WL * v.transpose();
  OUT.block(first_slope, first_slope, nslope_k, nslope_k) += N * uu * WK * WK.transpose();
  OUT.block(first_slope, first_slope + nslope_k, nslope_k, nslope_l) += N * uv * WK * WL.transpose();
  OUT.block(first_slope + nslope_k, first_slope, nslope_l, nslope_k) += N * uv * WL * WK.transpose();
  OUT.block(first_slope + nslope_k, first_slope + nslope_k, nslope_l, nslope_l) += N * vv * WL * WL.transpose();
}

void pairs::Pair::add_log_scale(const double N, const Eigen::Ref<const Eigen::VectorXd> SCALE, Eigen::VectorXd &SCORE,
                                Eigen::MatrixXd *INFO) const{
  const int row_k = nloc + 1 + nx, row_l = nloc + 1 + nslope_k + nx;
  if(numeric_k){
    if(SCORE.size() > 0) SCORE(row_k) += N / SCALE(k);
    if(INFO) (*INFO)(row_k, row_k) += N / (SCALE(k) * SCALE(k));
  }
  if(numeric_l){
    if(SCORE.size() > 0) SCORE(row_l) += N / SCALE(l);
    if(INFO) (*INFO)(row_l, row_l) += N / (SCALE(l) * SCALE(l));
  }
}

static void beta_curvature(const model::Matrices &PAR, const pairs::Pair &PAIR, const int ITEM, const Eigen::MatrixXd &D,
                           const Eigen::VectorXd &W, Eigen::MatrixXd &HESSIAN){
  const model::Structure &S = PAR.structure;
  const int nload = PAIR.first_corr, first_beta = PAIR.first_beta;
  const Eigen::VectorXd Aw = PAR.ib_inv * W;
  for(int b = 0; b < S.nreg; b++){
    const int i = S.beta_row[b], j = S.beta_col[b];
    const Eigen::VectorXd col = Aw(j) * D.row(i).head(nload).transpose();
    HESSIAN.block(0, first_beta + b, nload, 1) += col;
    HESSIAN.block(first_beta + b, 0, 1, nload) += col.transpose();
    for(int c = 0; c < S.nreg; c++){
      const int r = S.beta_row[c], t = S.beta_col[c];
      HESSIAN(first_beta + b, first_beta + c) += PAR.lambda_ib_inv(ITEM, r) * PAR.ib_inv(t, i) * Aw(j) +
        PAR.lambda_ib_inv(ITEM, i) * PAR.ib_inv(j, r) * Aw(t);
    }
  }
}

static Eigen::MatrixXd row_product_hessian(const model::Matrices &PAR, const pairs::Pair &PAIR, const Eigen::MatrixXd &D_A,
                                           const Eigen::MatrixXd &D_B, const Eigen::VectorXd &ROW_A, const Eigen::VectorXd &ROW_B,
                                           const Eigen::VectorXd &CHOL_A, const Eigen::VectorXd &CHOL_B,
                                           const Eigen::Ref<const Eigen::VectorXd> D_PRODUCT, const int ITEM_A, const int ITEM_B){
  const model::Structure &S = PAR.structure;
  const int q = S.q, first_corr = PAIR.first_corr, ncorr = S.ncorr, nvar = S.nvar;

  if(S.latcov == model::LatCov::direct){
    const Eigen::MatrixXd cross = D_A.transpose() * PAR.psi * D_B;
    Eigen::MatrixXd hessian = cross + cross.transpose();
    for(int a = 0; a < ncorr; a++){
      Eigen::VectorXd ea = Eigen::VectorXd::Zero(q), eb = Eigen::VectorXd::Zero(q);
      ea(S.corr_row[a]) += ROW_B(S.corr_col[a]); ea(S.corr_col[a]) += ROW_B(S.corr_row[a]);
      eb(S.corr_row[a]) += ROW_A(S.corr_col[a]); eb(S.corr_col[a]) += ROW_A(S.corr_row[a]);
      const Eigen::VectorXd col = D_A.transpose() * ea + D_B.transpose() * eb;
      hessian.col(first_corr + a) += col;
      hessian.row(first_corr + a) += col.transpose();
    }
    for(int v = 0; v < nvar; v++){
      const int h = S.logsd_factor[v], pos = first_corr + ncorr + v;
      const double w = 2 * PAR.sd(h) * PAR.sd(h);
      const Eigen::VectorXd col = w * (D_A.row(h).transpose() * ROW_B(h) + D_B.row(h).transpose() * ROW_A(h));
      hessian.col(pos) += col;
      hessian.row(pos) += col.transpose();
      hessian(pos, pos) += 2 * w * ROW_A(h) * ROW_B(h);
    }
    if(S.nreg > 0){
      beta_curvature(PAR, PAIR, ITEM_A, D_A, PAR.psi * ROW_B, hessian);
      beta_curvature(PAR, PAIR, ITEM_B, D_B, PAR.psi * ROW_A, hessian);
    }
    return hessian;
  }

  const Eigen::MatrixXd sd_d_a = PAR.sd.asDiagonal() * D_A, sd_d_b = PAR.sd.asDiagonal() * D_B;
  Eigen::MatrixXd d_chol_a = PAR.chol * sd_d_a, d_chol_b = PAR.chol * sd_d_b;
  for(int a = 0; a < ncorr; a++){
    d_chol_a.col(first_corr + a) = ROW_A(S.corr_col[a]) * PAR.dchol.col(a);
    d_chol_b.col(first_corr + a) = ROW_B(S.corr_col[a]) * PAR.dchol.col(a);
  }
  for(int v = 0; v < nvar; v++){
    d_chol_a.col(first_corr + ncorr + v) = ROW_A(S.logsd_factor[v]) * PAR.chol.col(S.logsd_factor[v]);
    d_chol_b.col(first_corr + ncorr + v) = ROW_B(S.logsd_factor[v]) * PAR.chol.col(S.logsd_factor[v]);
  }
  Eigen::MatrixXd hessian = d_chol_a.transpose() * d_chol_b;
  hessian += hessian.transpose().eval();

  const Eigen::VectorXd dchol_b = PAR.dchol.transpose() * CHOL_B, dchol_a = PAR.dchol.transpose() * CHOL_A;
  const Eigen::VectorXd cholt_b = PAR.chol.transpose() * CHOL_B, cholt_a = PAR.chol.transpose() * CHOL_A;
  for(int a = 0; a < ncorr; a++){
    const Eigen::VectorXd v = sd_d_a.row(S.corr_col[a]).transpose() * dchol_b(a) + sd_d_b.row(S.corr_col[a]).transpose() * dchol_a(a);
    hessian.col(first_corr + a) += v;
    hessian.row(first_corr + a) += v.transpose();
  }
  for(int s = 0; s < nvar; s++){
    const int h = S.logsd_factor[s];
    const Eigen::VectorXd v = sd_d_a.row(h).transpose() * cholt_b(h) + sd_d_b.row(h).transpose() * cholt_a(h);
    hessian.col(first_corr + ncorr + s) += v;
    hessian.row(first_corr + ncorr + s) += v.transpose();
  }
  for(int a = 0; a < ncorr; a++){
    for(int b = a; b < ncorr && S.corr_col[b] == S.corr_col[a]; b++){
      const int c = S.corr_col[a], i0 = S.corr_row[a], i1 = S.corr_row[b];
      Eigen::VectorXd d2_chol = Eigen::VectorXd::Zero(q);
      if(i0 < i1){
        d2_chol(i1) = -PAR.lkj_cumw(i1, c) * PAR.lkj_z(i0, c) * PAR.lkj_w(i1, c) * PAR.lkj_w(i1, c);
        for(int i = i1 + 1; i <= c; i++) d2_chol(i) = PAR.chol(i, c) * PAR.lkj_z(i0, c) * PAR.lkj_z(i1, c);
      }else{
        d2_chol(i0) = -2 * PAR.lkj_cumw(i0, c) * PAR.lkj_z(i0, c) * PAR.lkj_w(i0, c) * PAR.lkj_w(i0, c);
        for(int i = i0 + 1; i <= c; i++) d2_chol(i) = PAR.chol(i, c) * (2 * PAR.lkj_z(i0, c) * PAR.lkj_z(i0, c) - 1);
      }
      const double v = ROW_A(c) * d2_chol.dot(CHOL_B) + ROW_B(c) * d2_chol.dot(CHOL_A);
      hessian(first_corr + a, first_corr + b) += v;
      if(b != a) hessian(first_corr + b, first_corr + a) += v;
    }
  }
  for(int s = 0; s < nvar; s++){
    const int pos = first_corr + ncorr + s;
    hessian(pos, pos) += D_PRODUCT(pos);
    for(int a = 0; a < ncorr; a++){
      if(S.corr_col[a] != S.logsd_factor[s]) continue;
      hessian(first_corr + a, pos) += D_PRODUCT(first_corr + a);
      hessian(pos, first_corr + a) += D_PRODUCT(first_corr + a);
    }
  }
  if(S.nreg > 0){
    beta_curvature(PAR, PAIR, ITEM_A, D_A, PAR.sd.cwiseProduct(cholt_b), hessian);
    beta_curvature(PAR, PAIR, ITEM_B, D_B, PAR.sd.cwiseProduct(cholt_a), hessian);
  }
  return hessian;
}

Eigen::MatrixXd pairs::cov_hessian(const model::Matrices &PAR, const Pair &PAIR){
  return row_product_hessian(PAR, PAIR, PAIR.d_lambda_row_k, PAIR.d_lambda_row_l, PAIR.row_k, PAIR.row_l, PAIR.chol_row_k,
                             PAIR.chol_row_l, PAIR.d_cov, PAIR.k, PAIR.l);
}

Eigen::MatrixXd pairs::reduced_curvature(const model::Matrices &PAR, const Pair &PAIR, const Eigen::Ref<const Eigen::VectorXd> SCORE){
  const model::Structure &S = PAR.structure;
  const int nthr = PAIR.nthr, nlocal = PAIR.nlocal, ncov_local = PAIR.ncov_local;
  const double score_rho = SCORE(PAIR.nloc);
  auto mu_hessian = [&](const int ITEM, const Eigen::MatrixXd &D, const double A, const Eigen::VectorXd &WX){
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nlocal, nlocal);
    Eigen::MatrixXd h0 = Eigen::MatrixXd::Zero(ncov_local, ncov_local);
    const Eigen::VectorXd w = A * PAR.alpha + PAR.gamma * WX;
    beta_curvature(PAR, PAIR, ITEM, D, w, h0);
    H.topLeftCorner(ncov_local, ncov_local) = h0;
    for(int a = 0; a < S.nmean; a++){
      const Eigen::VectorXd col = A * D.row(S.alpha_factor[a]).transpose();
      H.block(0, PAIR.first_alpha + a, ncov_local, 1) += col;
      H.block(PAIR.first_alpha + a, 0, 1, ncov_local) += col.transpose();
    }
    for(int a = 0; a < S.ngam; a++){
      const Eigen::VectorXd col = WX(S.gamma_covariate[a]) * D.row(S.gamma_factor[a]).transpose();
      H.block(0, PAIR.first_gamma + a, ncov_local, 1) += col;
      H.block(PAIR.first_gamma + a, 0, 1, ncov_local) += col.transpose();
    }
    return H;
  };
  auto scale_hessian = [&](const int ITEM, const Eigen::MatrixXd &D, const Eigen::VectorXd &ROW, const Eigen::VectorXd &CHOL_ROW,
                           const Eigen::VectorXd &D_VAR, const int SCALE_POS){
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nlocal, nlocal);
    const double sc = PAR.scale(ITEM);
    if(!S.scale_from_rsd[ITEM]){
      if(SCALE_POS >= 0) H(SCALE_POS, SCALE_POS) = sc;
      return H;
    }
    Eigen::VectorXd gt = D_VAR;
    H.topLeftCorner(ncov_local, ncov_local) = row_product_hessian(PAR, PAIR, D, D, ROW, ROW, CHOL_ROW, CHOL_ROW,
                                                                  D_VAR.head(ncov_local), ITEM, ITEM);
    if(SCALE_POS >= 0){
      gt(SCALE_POS) += 2 * PAR.residual_var(ITEM);
      H(SCALE_POS, SCALE_POS) += 4 * PAR.residual_var(ITEM);
    }
    return Eigen::MatrixXd(0.75 * std::pow(sc, 5) * gt * gt.transpose() - 0.5 * std::pow(sc, 3) * H);
  };

  const Eigen::MatrixXd H_scale_k = scale_hessian(PAIR.k, PAIR.d_lambda_row_k, PAIR.row_k, PAIR.chol_row_k, PAIR.d_var_k, PAIR.scale_pos_k);
  const Eigen::MatrixXd H_scale_l = scale_hessian(PAIR.l, PAIR.d_lambda_row_l, PAIR.row_l, PAIR.chol_row_l, PAIR.d_var_l, PAIR.scale_pos_l);
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(nthr + nlocal, nthr + nlocal);
  auto add_item_terms = [&](const int ITEM, const int TAU_OFFSET, const int CAT, const int FIRST_LOC, const int NLOC,
                            const int FIRST_THR, const int FIRST_SLOPE, const bool NUMERIC, const Eigen::MatrixXd &D_LAMBDA_ROW,
                            const Eigen::VectorXd &D_MU, const Eigen::MatrixXd &D_PI, const Eigen::VectorXd &D_SCALE,
                            const Eigen::MatrixXd &H_SCALE){
    const Eigen::VectorXd score_loc = SCORE.segment(FIRST_LOC, NLOC), score_slope = SCORE.segment(FIRST_SLOPE, S.nx);
    const double sum_loc = score_loc.sum();
    double weight_scale = PAR.pi.row(ITEM).dot(score_slope);
    if(NUMERIC){
      weight_scale += -PAR.mu(ITEM) * score_loc(0) + SCORE(FIRST_SLOPE + S.nx);
    }else{
      weight_scale += score_loc.dot((PAR.tau.segment(TAU_OFFSET, CAT - 1).array() - PAR.mu(ITEM)).matrix());
    }
    Eigen::VectorXd e = Eigen::VectorXd::Zero(nthr + nlocal), g = Eigen::VectorXd::Zero(nthr + nlocal);
    if(!NUMERIC) e.segment(FIRST_THR, CAT - 1) = score_loc;
    e.tail(nlocal) = -sum_loc * D_MU + D_PI * score_slope;
    g.tail(nlocal) = D_SCALE;
    W += g * e.transpose() + e * g.transpose();
    W.bottomRightCorner(nlocal, nlocal) += weight_scale * H_SCALE + PAR.scale(ITEM) * mu_hessian(ITEM, D_LAMBDA_ROW, -sum_loc, score_slope);
  };
  add_item_terms(PAIR.k, PAIR.tau_offset_k, PAIR.cat_k, 0, PAIR.nloc_k, 0, PAIR.nloc + 1, PAIR.numeric_k, PAIR.d_lambda_row_k,
                 PAIR.d_mu_k, PAIR.d_pi_k, PAIR.d_scale_k, H_scale_k);
  add_item_terms(PAIR.l, PAIR.tau_offset_l, PAIR.cat_l, PAIR.nloc_k, PAIR.nloc_l, PAIR.cat_k - 1, PAIR.nloc + 1 + PAIR.nslope_k,
                 PAIR.numeric_l, PAIR.d_lambda_row_l, PAIR.d_mu_l, PAIR.d_pi_l, PAIR.d_scale_l, H_scale_l);
  if(score_rho != 0){
    const double sk = PAR.scale(PAIR.k), sl = PAR.scale(PAIR.l), c = PAIR.cov;
    const Eigen::VectorXd &gk = PAIR.d_scale_k, &gl = PAIR.d_scale_l, &gc = PAIR.d_cov_full;
    Eigen::MatrixXd R = c * (gk * gl.transpose()) + sl * (gk * gc.transpose()) + sk * (gl * gc.transpose());
    R += R.transpose().eval();
    R += c * sl * H_scale_k + c * sk * H_scale_l;
    R.topLeftCorner(ncov_local, ncov_local) += sk * sl * cov_hessian(PAR, PAIR);
    W.bottomRightCorner(nlocal, nlocal) += score_rho * R;
  }
  return W;
}

double pairs::unit_loglik(const Reduced &REDUCED, const GroupData &DATA, const Pair &PAIR, Grid &GRID, const int UNIT, const Order ORDER,
                          UnitWork &WORK){
  const Reduced &R = REDUCED;
  UnitWork &W = WORK;
  const int k = PAIR.k, l = PAIR.l, t = DATA.stratum(UNIT);
  const double yk = DATA.responses(UNIT, k), yl = DATA.responses(UNIT, l);
  const double shift_k = DATA.strata.col(t).dot(R.slope.row(k).transpose()), shift_l = DATA.strata.col(t).dot(R.slope.row(l).transpose());
  if(!PAIR.numeric_k && !PAIR.numeric_l){
    const int sk = yk, sl = yl;
    GRID.fill(R.tau, PAIR.tau_offset_k, PAIR.tau_offset_l, shift_k, shift_l, std::max(sk, 1), sk + 1, std::max(sl, 1), sl + 1);
    const double pe = GRID.cell_probability(sk, sl, W.cell_score, W.cell_hessian);
    if(ORDER != Order::value) W.score = W.cell_score / pe;
    if(ORDER == Order::hessian){
      W.info.noalias() = W.cell_score * W.cell_score.transpose() / (pe * pe);
      W.info -= W.cell_hessian / pe;
    }
    return std::log(pe);
  }
  const double rho = PAIR.rho, r2 = 1 - rho * rho, r = std::sqrt(r2), log2pi = std::log(2 * M_PI);
  if(PAIR.numeric_k && PAIR.numeric_l){
    const double zk = R.scale(k) * yk + R.location(k) - shift_k, zl = R.scale(l) * yl + R.location(l) - shift_l;
    const double q = zk * zk - 2 * rho * zk * zl + zl * zl, r4 = r2 * r2;
    if(ORDER != Order::value){
      W.score.resize(3);
      W.score << -(zk - rho * zl) / r2, -(zl - rho * zk) / r2, (rho + zk * zl) / r2 - rho * q / r4;
    }
    if(ORDER == Order::hessian){
      W.info.resize(3, 3);
      W.info(0, 0) = W.info(1, 1) = 1 / r2;
      W.info(0, 1) = W.info(1, 0) = -rho / r2;
      W.info(0, 2) = W.info(2, 0) = -zl / r2 + 2 * rho * (zk - rho * zl) / r4;
      W.info(1, 2) = W.info(2, 1) = -zk / r2 + 2 * rho * (zl - rho * zk) / r4;
      W.info(2, 2) = -((1 + rho * rho) / r4 + 4 * rho * zk * zl / r4 - q / r4 - 4 * rho * rho * q / (r4 * r2));
    }
    return -log2pi - 0.5 * std::log(r2) - 0.5 * q / r2 + std::log(R.scale(k)) + std::log(R.scale(l));
  }
  const int numeric_item = PAIR.numeric_k ? k : l;
  const int cat = PAIR.numeric_k ? PAIR.cat_l : PAIR.cat_k, tau_offset = PAIR.numeric_k ? PAIR.tau_offset_l : PAIR.tau_offset_k;
  const int category = PAIR.numeric_k ? yl : yk;
  const double z = R.scale(numeric_item) * (PAIR.numeric_k ? yk : yl) + R.location(numeric_item) - (PAIR.numeric_k ? shift_k : shift_l);
  const double shift_other = PAIR.numeric_k ? shift_l : shift_k;
  const bool top = category == cat - 1, bot = category == 0;
  const double hi = top ? 0 : R.tau(tau_offset + category) - shift_other, lo = bot ? 0 : R.tau(tau_offset + category - 1) - shift_other;
  const double uh = (hi - rho * z) / r, ul = (lo - rho * z) / r;
  const double ph = top ? 0 : R::dnorm(uh, 0, 1, 0), pl = bot ? 0 : R::dnorm(ul, 0, 1, 0);
  const double pr = (top ? 1 : R::pnorm(uh, 0, 1, 1, 0)) - (bot ? 0 : R::pnorm(ul, 0, 1, 1, 0)) + 1e-8;
  const double r3 = r2 * r;
  Eigen::Vector4d du, dl;
  du << -rho / r, 0, 1 / r, (rho * hi - z) / r3;
  dl << -rho / r, 1 / r, 0, (rho * lo - z) / r3;
  const Eigen::Vector4d dd = ph * du - pl * dl;
  const int first_thr = PAIR.numeric_k ? 1 : 0, pos_z = PAIR.numeric_k ? 0 : PAIR.nloc_k;
  const int pos[4] = {pos_z, bot ? -1 : first_thr + category - 1, top ? -1 : first_thr + category, PAIR.nloc};
  if(ORDER != Order::value){
    W.score = Eigen::VectorXd::Zero(PAIR.nloc + 1);
    for(int a = 0; a < 4; a++) if(pos[a] >= 0) W.score(pos[a]) += dd(a) / pr;
    W.score(pos_z) -= z;
  }
  if(ORDER == Order::hessian){
    Eigen::Matrix4d H = Eigen::Matrix4d::Zero();
    if(!top){
      H += -ph * uh * du * du.transpose();
      H(0, 3) += -ph / r3;
      H(3, 0) += -ph / r3;
      H(2, 3) += ph * rho / r3;
      H(3, 2) += ph * rho / r3;
      H(3, 3) += ph * (hi / r3 + 3 * rho * (rho * hi - z) / (r3 * r2));
    }
    if(!bot){
      H += pl * ul * dl * dl.transpose();
      H(0, 3) += pl / r3;
      H(3, 0) += pl / r3;
      H(1, 3) += -pl * rho / r3;
      H(3, 1) += -pl * rho / r3;
      H(3, 3) += -pl * (lo / r3 + 3 * rho * (rho * lo - z) / (r3 * r2));
    }
    const Eigen::Matrix4d M4 = dd * dd.transpose() / (pr * pr) - H / pr;
    W.info = Eigen::MatrixXd::Zero(PAIR.nloc + 1, PAIR.nloc + 1);
    for(int a = 0; a < 4; a++){
      for(int b = 0; b < 4; b++) if(pos[a] >= 0 && pos[b] >= 0) W.info(pos[a], pos[b]) += M4(a, b);
    }
    W.info(pos_z, pos_z) += 1;
  }
  return -0.5 * log2pi - 0.5 * z * z + std::log(R.scale(numeric_item)) + std::log(pr);
}

void pairs::pair_loglik(const model::Items &ITEMS, const Reduced &REDUCED, const GroupData &DATA, const Pair &PAIR, const Order ORDER,
                        double &LL, Eigen::VectorXd &SCORE, Eigen::MatrixXd &INFO, Eigen::MatrixXd *CELL_SCORES){
  const int nscore = PAIR.nloc + 1, ncells = PAIR.cat_k * PAIR.cat_l;
  const bool derivatives = ORDER != Order::value, hessian = ORDER == Order::hessian;
  UnitWork W;
  Grid grid(PAIR.cat_k, PAIR.cat_l, PAIR.rho, ORDER);
  if(DATA.tables && !PAIR.numeric_k && !PAIR.numeric_l){
    Eigen::VectorXd g, score_t(derivatives ? nscore : 0), cell_full(CELL_SCORES ? PAIR.nreduced : 0);
    Eigen::MatrixXd cell_hessian, info_t(hessian ? nscore : 0, hessian ? nscore : 0);
    for(int t = 0; t < DATA.strata.cols(); t++){
      grid.fill(REDUCED.tau, PAIR.tau_offset_k, PAIR.tau_offset_l, DATA.strata.col(t).dot(REDUCED.slope.row(PAIR.k).transpose()),
                DATA.strata.col(t).dot(REDUCED.slope.row(PAIR.l).transpose()), 1, PAIR.cat_k, 1, PAIR.cat_l);
      PAIR.covariate_weights(DATA.strata.col(t), 0, 0, W.weights_k, W.weights_l);
      const int o = t * ITEMS.ncell + PAIR.cell_offset;
      if(derivatives) score_t.setZero();
      if(hessian) info_t.setZero();
      for(int sk = 0; sk < PAIR.cat_k; sk++){
        for(int sl = 0; sl < PAIR.cat_l; sl++){
          const int c = sk * PAIR.cat_l + sl;
          const double n = DATA.freq(o + c);
          if(n == 0 && !CELL_SCORES) continue;
          const double pe = grid.cell_probability(sk, sl, g, cell_hessian);
          if(CELL_SCORES){
            cell_full.setZero();
            PAIR.add_score(1 / pe, g, W.weights_k, W.weights_l, cell_full);
            CELL_SCORES->col(t * ncells + c) = PAIR.to_local(cell_full);
          }
          if(n == 0) continue;
          LL += n * std::log(pe);
          if(derivatives) score_t += n / pe * g;
          if(hessian) info_t += n * (g * g.transpose() / (pe * pe) - cell_hessian / pe);
        }
      }
      if(derivatives) PAIR.add_score(1, score_t, W.weights_k, W.weights_l, SCORE);
      if(hessian) PAIR.add_info(1, info_t, W.weights_k, W.weights_l, INFO);
    }
    return;
  }
  for(int u = 0; u < DATA.responses.rows(); u++){
    const double yk = DATA.responses(u, PAIR.k), yl = DATA.responses(u, PAIR.l);
    if(std::isnan(yk) || std::isnan(yl)) continue;
    const double n = DATA.weights(u);
    LL += n * unit_loglik(REDUCED, DATA, PAIR, grid, u, ORDER, W);
    if(!derivatives) continue;
    PAIR.covariate_weights(DATA.strata.col(DATA.stratum(u)), yk, yl, W.weights_k, W.weights_l);
    PAIR.add_score(n, W.score, W.weights_k, W.weights_l, SCORE);
    if(hessian) PAIR.add_info(n, W.info, W.weights_k, W.weights_l, INFO);
    PAIR.add_log_scale(n, REDUCED.scale, SCORE, hessian ? &INFO : nullptr);
  }
}

void pairs::add_pair_gradient(const model::Matrices &PAR, const GroupData &DATA, const int INDEX, const Order ORDER, double &LL,
                              Eigen::VectorXd &GRADIENT){
  const bool gradient = ORDER == Order::gradient;
  const Pair pair(PAR, INDEX, gradient);
  Eigen::VectorXd score = Eigen::VectorXd::Zero(gradient ? pair.nreduced : 0);
  Eigen::MatrixXd info;
  pair_loglik(PAR.structure.items, Reduced(PAR), DATA, pair, ORDER, LL, score, info);
  if(!gradient) return;
  const std::vector<int> I = pair.parameter_indices();
  const Eigen::VectorXd local = pair.to_local(score);
  for(std::size_t a = 0; a < I.size(); a++) GRADIENT(I[a]) += local(a);
}

void pairs::univariate_margins(const Eigen::Ref<const Eigen::VectorXd> TAU, const model::Items &ITEMS,
                               const Eigen::Ref<const Eigen::VectorXd> UFREQ, const Order ORDER, double &LL,
                               Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &SCORES,
                               const Eigen::Ref<const Eigen::VectorXd> ITEM_WEIGHTS){
  if(UFREQ.size() == 0) return;
  for(int k = 0, uk = 0; k < ITEMS.p; uk += ITEMS.cat(k), k++){
    const int ck = ITEMS.cat(k), tk = ITEMS.tau_offset(k);
    const double wk = ITEM_WEIGHTS.size() > 0 ? ITEM_WEIGHTS(k) : 1;
    if(wk == 0) continue;
    for(int a = 0; a < ck; a++){
      const double nc = wk * UFREQ(uk + a);
      if(nc == 0) continue;
      const bool up = a < ck - 1, dw = a > 0;
      const double hi = up ? TAU(tk + a) : 0, lo = dw ? TAU(tk + a - 1) : 0;
      const double pe = (up ? R::pnorm(hi, 0, 1, 1, 0) : 1) - (dw ? R::pnorm(lo, 0, 1, 1, 0) : 0) + 1e-8;
      LL += nc * std::log(pe);
      if(ORDER == Order::value) continue;

      const double ghi = up ? R::dnorm(hi, 0, 1, 0) : 0, glo = dw ? -R::dnorm(lo, 0, 1, 0) : 0;
      if(up) GRADIENT(tk + a) += nc * ghi / pe;
      if(dw) GRADIENT(tk + a - 1) += nc * glo / pe;
      if(SCORES.size() > 0){
        if(up) SCORES(tk + a, uk + a) = ghi / pe;
        if(dw) SCORES(tk + a - 1, uk + a) = glo / pe;
      }
      if(ORDER == Order::hessian){
        if(up) HESSIAN(tk + a, tk + a) += nc * (ghi * ghi / (pe * pe) + hi * ghi / pe);
        if(dw) HESSIAN(tk + a - 1, tk + a - 1) += nc * (glo * glo / (pe * pe) + lo * glo / pe);
        if(up && dw){
          HESSIAN(tk + a, tk + a - 1) += nc * ghi * glo / (pe * pe);
          HESSIAN(tk + a - 1, tk + a) += nc * ghi * glo / (pe * pe);
        }
      }
    }
  }
}

void pairs::add_univariate(const model::Matrices &PAR, const GroupData &DATA, const Order ORDER, double &LL, Eigen::VectorXd &GRADIENT,
                           Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &SCORES, const Eigen::Ref<const Eigen::VectorXd> ITEM_WEIGHTS){
  const model::Structure &S = PAR.structure;
  const model::Items &items = S.items;
  if(S.standardised || DATA.ufreq.size() == 0){
    univariate_margins(PAR.tau_std, items, DATA.ufreq, ORDER, LL, GRADIENT, HESSIAN, SCORES, ITEM_WEIGHTS);
    return;
  }
  const bool derivatives = ORDER != Order::value, hessian = ORDER == Order::hessian;
  const int nthr = items.nthr, ucells = items.cat.sum();
  std::vector<Pair> item_views;
  item_views.reserve(derivatives ? S.p : 0);
  for(int k = 0; k < S.p && derivatives; k++) item_views.emplace_back(PAR, -1, true, k);
  for(int t = 0; t < DATA.strata.cols(); t++){
    const auto xt = DATA.strata.col(t);
    const Eigen::VectorXd shift = PAR.pi_std * xt;
    Eigen::VectorXd tau_t = PAR.tau_std;
    for(int k = 0; k < S.p; k++) tau_t.segment(items.tau_offset(k), items.cat(k) - 1).array() -= shift(k);
    Eigen::VectorXd g = Eigen::VectorXd::Zero(nthr);
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(hessian ? nthr : 0, hessian ? nthr : 0);
    Eigen::MatrixXd us = Eigen::MatrixXd::Zero(SCORES.size() > 0 ? nthr : 0, ucells);
    univariate_margins(tau_t, items, DATA.ufreq.segment(t * ucells, ucells), ORDER, LL, g, h, us, ITEM_WEIGHTS);
    if(!derivatives) continue;
    for(int k = 0, uk = 0; k < S.p; uk += items.cat(k), k++){
      const Pair &item = item_views[k];
      const int m = item.cat_k - 1, tk = item.tau_offset_k;
      const std::vector<int> I = item.parameter_indices();
      Eigen::MatrixXd Jk = item.jacobian.topRows(m);
      if(S.nx > 0) Jk.rowwise() -= Eigen::RowVectorXd(xt.transpose() * item.jacobian.middleRows(item.nloc + 1, S.nx));
      const Eigen::VectorXd gl = Jk.transpose() * g.segment(tk, m);
      for(std::size_t a = 0; a < I.size(); a++) GRADIENT(I[a]) += gl(a);
      for(int c = 0; c < items.cat(k) && us.size() > 0; c++){
        const Eigen::VectorXd sl = Jk.transpose() * us.block(tk, uk + c, m, 1);
        for(std::size_t a = 0; a < I.size(); a++) SCORES(I[a], t * ucells + uk + c) += sl(a);
      }
      if(!hessian) continue;
      Eigen::VectorXd gu = Eigen::VectorXd::Zero(item.nreduced), gt = Eigen::VectorXd::Zero(item.nloc + 1), wk, wl;
      gt.head(m) = g.segment(tk, m);
      item.covariate_weights(xt, 0, 0, wk, wl);
      item.add_score(1, gt, wk, wl, gu);
      const Eigen::MatrixXd HL = Jk.transpose() * h.block(tk, tk, m, m) * Jk - reduced_curvature(PAR, item, gu);
      for(std::size_t a = 0; a < I.size(); a++){
        for(std::size_t b = 0; b < I.size(); b++) HESSIAN(I[a], I[b]) += HL(a, b);
      }
    }
  }
}

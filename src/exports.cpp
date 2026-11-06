#include "workers.h"

using pairs::Order;

template <typename T>
static T get(const Rcpp::List &X, const char *NAME){
  return Rcpp::as<T>(X[NAME]);
}

// [[Rcpp::export]]
Rcpp::List cpp_pl(const Eigen::Map<Eigen::VectorXd> PARAMETER_VECTOR, const Rcpp::List MODEL, const Rcpp::List DATA, const int ORDER = 0){
  const model::Groups groups(MODEL, DATA);
  const std::vector<model::Matrices> pars = groups.matrices(PARAMETER_VECTOR);
  const Order order = static_cast<Order>(ORDER);
  double ll = 0;
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(groups.d);
  for(int g = 0; g < groups.ngroups; g++){
    const pairs::GroupData data(Rcpp::as<Rcpp::List>(DATA[g]));
    std::vector<int> pair_indices(groups.structure[g].items.npairs);
    std::iota(pair_indices.begin(), pair_indices.end(), 0);
    workers::PairLoglik w(pars[g], data, pair_indices, order);
    RcppParallel::parallelReduce(0, pair_indices.size(), w);
    Eigen::MatrixXd hd, ud;
    pairs::add_univariate(pars[g], data, order, w.ll, w.gradient, hd, ud);
    ll += w.ll;
    if(order != Order::value) gradient.segment(groups.parameter_offset[g], groups.structure[g].d) = w.gradient;
  }
  return Rcpp::List::create(Rcpp::Named("pl") = ll, Rcpp::Named("gradient") = gradient);
}

// [[Rcpp::export]]
Rcpp::List cpp_sandwich(const Eigen::Map<Eigen::VectorXd> PARAMETER_VECTOR, const Rcpp::List MODEL, const Rcpp::List DATA,
                        const bool PAIRVAR = false){
  const model::Groups groups(MODEL, DATA);
  const std::vector<model::Matrices> pars = groups.matrices(PARAMETER_VECTOR);
  const int d = groups.d;
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(d, d), J = Eigen::MatrixXd::Zero(d, d);
  Eigen::MatrixXd pairvar = PAIRVAR ? Eigen::MatrixXd::Zero(d, d) : Eigen::MatrixXd();
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
  for(int g = 0; g < groups.ngroups; g++){
    const model::Structure &S = groups.structure[g];
    const model::Items &items = S.items;
    const int off = groups.parameter_offset[g];
    const Rcpp::List data_g = DATA[g];
    const pairs::GroupData data(data_g);
    std::vector<Eigen::MatrixXd> cell_scores(items.npairs);
    std::vector<std::vector<int>> score_index(items.npairs);
    std::vector<Eigen::VectorXd> pair_gradient(items.npairs);
    workers::PairHessian hw(pars[g], data, cell_scores, score_index, pair_gradient);
    RcppParallel::parallelReduce(0, items.npairs, hw);
    double ll = 0;
    Eigen::MatrixXd univariate_scores = Eigen::MatrixXd::Zero(S.d, data.ufreq.size());
    pairs::add_univariate(pars[g], data, Order::hessian, ll, hw.gradient, hw.hessian, univariate_scores);

    const Eigen::MatrixXd cluster_weights = get<Eigen::MatrixXd>(data_g, "CLUSTER_WEIGHTS");
    const Eigen::VectorXd uweights = get<Eigen::VectorXd>(data_g, "UWEIGHTS"), jweights = get<Eigen::VectorXd>(data_g, "JWEIGHTS");
    Eigen::MatrixXd case_scores = Eigen::MatrixXd::Zero(S.d, data.responses.rows());
    workers::CaseScores cw(pars[g], data, cell_scores, score_index, univariate_scores, uweights, case_scores);
    RcppParallel::parallelFor(0, data.responses.rows(), cw, 256);
    const Eigen::MatrixXd weighted = cluster_weights.size() > 0 ? Eigen::MatrixXd(case_scores * cluster_weights) :
      Eigen::MatrixXd(case_scores * jweights.asDiagonal());
    H.block(off, off, S.d, S.d) = hw.hessian / groups.n;
    J.block(off, off, S.d, S.d) = weighted * weighted.transpose() / groups.n;
    gradient.segment(off, S.d) = hw.gradient;
    if(PAIRVAR){
      Eigen::MatrixXd univariate_items = Eigen::MatrixXd::Zero(S.d, S.p), pair_items = Eigen::MatrixXd::Zero(S.d, S.p);
      Eigen::MatrixXd pv = Eigen::MatrixXd::Zero(S.d, S.d);
      const int ucells = items.cat.sum();
      for(int t = 0; t < data.strata.cols() && data.ufreq.size() > 0; t++){
        for(int k = 0, uk = 0; k < S.p; uk += items.cat(k), k++){
          univariate_items.col(k) += univariate_scores.middleCols(t * ucells + uk, items.cat(k)) * data.ufreq.segment(t * ucells + uk, items.cat(k));
        }
      }
      for(int r = 0; r < items.npairs; r++){
        const std::vector<int> &I = score_index[r];
        const Eigen::VectorXd &s = pair_gradient[r];
        for(std::size_t a = 0; a < I.size(); a++){
          pair_items(I[a], items.pair_k[r]) += s(a);
          pair_items(I[a], items.pair_l[r]) += s(a);
          for(std::size_t b = 0; b < I.size(); b++) pv(I[a], I[b]) += s(a) * s(b);
        }
      }
      const Eigen::MatrixXd cross = pair_items * univariate_items.transpose();
      const Eigen::VectorXd u = univariate_items.rowwise().sum();
      const double m = S.p - 1.0;
      pv += (cross + cross.transpose()) / m + ((S.p - 2.0) * univariate_items * univariate_items.transpose() + u * u.transpose()) / (m * m);
      pairvar.block(off, off, S.d, S.d) = pv;
    }
  }
  return Rcpp::List::create(Rcpp::Named("H") = H, Rcpp::Named("J") = J, Rcpp::Named("gradient") = gradient,
                            Rcpp::Named("pairvar") = pairvar);
}

// [[Rcpp::export]]
Rcpp::List cpp_reduced_form(const Eigen::Map<Eigen::VectorXd> PARAMETER_VECTOR, const Rcpp::List MODEL, const Rcpp::List DATA){
  const model::Groups groups(MODEL, DATA);
  const std::vector<model::Matrices> pars = groups.matrices(PARAMETER_VECTOR);
  int nitem = 0;
  for(int g = 0; g < groups.ngroups; g++) nitem += pairs::ReducedIndex(groups.structure[g].items, groups.structure[g].nx).size;
  Eigen::VectorXd rho(groups.npairs), tau(nitem);
  Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(groups.npairs, groups.d), tau_jacobian = Eigen::MatrixXd::Zero(nitem, groups.d);
  auto put = [](const pairs::Pair &PAIR, const int ROW, const int OFFSET, const int E, Eigen::MatrixXd &OUT){
    const std::vector<int> I = PAIR.parameter_indices();
    const Eigen::VectorXd v = PAIR.to_local(Eigen::VectorXd::Unit(PAIR.nreduced, E));
    for(std::size_t a = 0; a < I.size(); a++) OUT(ROW, OFFSET + I[a]) += v(a);
  };
  for(int g = 0, ot = 0; g < groups.ngroups; g++){
    const model::Structure &S = groups.structure[g];
    const model::Matrices &par = pars[g];
    const pairs::ReducedIndex index(S.items, S.nx);
    const int off = groups.parameter_offset[g];
    tau.segment(ot, S.items.nthr) = par.tau_std;
    for(int r = 0; r < S.items.npairs; r++){
      const pairs::Pair pair(par, r, true);
      rho(groups.pair_offset[g] + r) = pair.rho;
      put(pair, groups.pair_offset[g] + r, off, pair.nloc, jacobian);
    }
    for(int k = 0; k < S.p; k++){
      const pairs::Pair item(par, -1, true, k);
      for(int a = 0; a < item.cat_k - 1; a++) put(item, ot + S.items.tau_offset(k) + a, off, a, tau_jacobian);
      for(int c = 0; c < S.nx; c++){
        const int row = ot + index.slope(c, k);
        tau(row) = par.pi_std(k, c);
        put(item, row, off, item.nloc + 1 + c, tau_jacobian);
      }
      if(item.numeric_k){
        tau(ot + index.location(k)) = par.location(k);
        tau(ot + index.scale(k)) = par.scale(k);
        put(item, ot + index.location(k), off, 0, tau_jacobian);
        put(item, ot + index.scale(k), off, item.nloc + 1 + S.nx, tau_jacobian);
      }
    }
    ot += index.size;
  }
  return Rcpp::List::create(Rcpp::Named("rho") = rho, Rcpp::Named("jacobian") = jacobian, Rcpp::Named("tau") = tau,
                            Rcpp::Named("tau_jacobian") = tau_jacobian);
}

// [[Rcpp::export]]
Rcpp::List cpp_saturated(const Eigen::Map<Eigen::VectorXd> ITEM_PARAMETERS, const Eigen::Map<Eigen::VectorXd> RHO,
                         const Rcpp::List DATA, const int ORDER, const bool SCORES){
  const model::Items items(get<Eigen::VectorXd>(DATA, "CAT"));
  const pairs::GroupData data(DATA);
  const Order order = static_cast<Order>(SCORES ? std::max(ORDER, 1) : ORDER);
  const int nx = data.strata.rows(), nthr = items.nthr, nitem = ITEM_PARAMETERS.size(), ucells = items.cat.sum();
  const pairs::ReducedIndex index(items, nx);
  const Eigen::VectorXd tau = ITEM_PARAMETERS.head(nthr);
  const Eigen::MatrixXd slope = Eigen::Map<const Eigen::MatrixXd>(ITEM_PARAMETERS.data() + nthr, items.p, nx);
  Eigen::VectorXd location = Eigen::VectorXd::Zero(items.p), scale = Eigen::VectorXd::Ones(items.p);
  for(int k = 0; k < items.p; k++){
    if(items.numeric_index(k) < 0) continue;
    location(k) = ITEM_PARAMETERS(index.location(k));
    scale(k) = ITEM_PARAMETERS(index.scale(k));
  }
  const pairs::Reduced reduced(tau, slope, location, scale);

  const bool keep = SCORES && data.tables;
  Eigen::MatrixXd pattern_scores = Eigen::MatrixXd::Zero(keep ? 2 * items.cat.maxCoeff() - 1 : 0, keep ? data.strata.cols() * items.ncell : 0);
  workers::SaturatedLoglik sat(items, data, reduced, RHO, nitem, order, SCORES, pattern_scores);
  RcppParallel::parallelReduce(0, items.npairs, sat);
  Eigen::MatrixXd univariate_scores = Eigen::MatrixXd::Zero(SCORES ? nitem : 0, SCORES ? data.ufreq.size() : 0);
  for(int t = 0; t < data.strata.cols() && data.ufreq.size() > 0; t++){
    const auto xt = data.strata.col(t);
    const Eigen::VectorXd shift = slope * xt;
    Eigen::VectorXd tau_t = tau;
    for(int k = 0; k < items.p; k++) tau_t.segment(items.tau_offset(k), items.cat(k) - 1).array() -= shift(k);
    Eigen::VectorXd gu = Eigen::VectorXd::Zero(nthr);
    Eigen::MatrixXd hu = Eigen::MatrixXd::Zero(order == Order::hessian ? nthr : 0, order == Order::hessian ? nthr : 0);
    Eigen::MatrixXd us = Eigen::MatrixXd::Zero(SCORES ? nthr : 0, SCORES ? ucells : 0);
    pairs::univariate_margins(tau_t, items, data.ufreq.segment(t * ucells, ucells), order, sat.ll, gu, hu, us);
    if(order == Order::value) continue;
    for(int k = 0, uk = 0; k < items.p; uk += items.cat(k), k++){
      const int m = items.cat(k) - 1, tk = items.tau_offset(k);
      std::vector<int> idx(m + nx);
      for(int a = 0; a < m; a++) idx[a] = tk + a;
      for(int c = 0; c < nx; c++) idx[m + c] = index.slope(c, k);
      Eigen::MatrixXd Jk = Eigen::MatrixXd::Zero(m, m + nx);
      Jk.leftCols(m).setIdentity();
      Jk.rightCols(nx).rowwise() = -xt.transpose();
      const Eigen::VectorXd gl = Jk.transpose() * gu.segment(tk, m);
      for(int a = 0; a < m + nx; a++) sat.gradient(idx[a]) += gl(a);
      for(int c = 0; c < items.cat(k) && SCORES; c++){
        const Eigen::VectorXd sl = Jk.transpose() * us.block(tk, uk + c, m, 1);
        for(int a = 0; a < m + nx; a++) univariate_scores(idx[a], t * ucells + uk + c) += sl(a);
      }
      if(order != Order::hessian) continue;
      const Eigen::MatrixXd HL = Jk.transpose() * hu.block(tk, tk, m, m) * Jk;
      for(int a = 0; a < m + nx; a++){
        for(int b = 0; b < m + nx; b++){
          sat.htt(idx[a], idx[b]) += HL(a, b);
          sat.hschur(idx[a], idx[b]) += HL(a, b);
        }
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("ll") = sat.ll, Rcpp::Named("gradient") = sat.gradient,
    Rcpp::Named("htt") = sat.htt, Rcpp::Named("htr") = sat.htr, Rcpp::Named("hrr") = sat.hrr,
    Rcpp::Named("hschur") = sat.hschur
  );
  if(SCORES){
    const Eigen::VectorXd uweights = get<Eigen::VectorXd>(DATA, "UWEIGHTS");
    Eigen::MatrixXd score_tau = Eigen::MatrixXd::Zero(nitem, data.responses.rows());
    Eigen::MatrixXd score_rho = Eigen::MatrixXd::Zero(items.npairs, data.responses.rows());
    workers::SaturatedCaseScores cases(items, data, reduced, RHO, pattern_scores, univariate_scores, uweights, score_tau, score_rho);
    RcppParallel::parallelFor(0, data.responses.rows(), cases, 256);
    out["st"] = score_tau;
    out["sr"] = score_rho;
  }
  return out;
}

// [[Rcpp::export]]
Eigen::VectorXd cpp_saturated_traces(
    const Eigen::Map<Eigen::MatrixXd> HTT,
    const Eigen::Map<Eigen::MatrixXd> HTR,
    const Eigen::Map<Eigen::VectorXd> HRR,
    const Eigen::Map<Eigen::MatrixXd> HSCHUR,
    const Eigen::Map<Eigen::MatrixXd> SCORE_TAU,
    const Eigen::Map<Eigen::MatrixXd> SCORE_RHO,
    const double N,
    const bool FULL = false
){
  const Eigen::VectorXd dinv = HRR.cwiseInverse();
  Eigen::MatrixXd M;
  if(FULL){
    const int rows = SCORE_TAU.rows() + SCORE_RHO.rows();
    Eigen::MatrixXd S(rows, SCORE_TAU.cols()), X(rows, SCORE_TAU.cols());
    S << SCORE_TAU, SCORE_RHO;
    const Eigen::MatrixXd xt = HSCHUR.ldlt().solve(SCORE_TAU - HTR * (dinv.asDiagonal() * SCORE_RHO));
    X << xt, dinv.asDiagonal() * (SCORE_RHO - HTR.transpose() * xt);
    M = (S.cols() <= S.rows()) ? Eigen::MatrixXd(S.transpose() * X / N) : Eigen::MatrixXd(X * S.transpose() / N);
  }else{
    const Eigen::MatrixXd U = SCORE_RHO - HTR.transpose() * HTT.ldlt().solve(SCORE_TAU);
    const Eigen::MatrixXd Ud = dinv.asDiagonal() * U;
    const Eigen::MatrixXd AinvU = Ud + dinv.asDiagonal() * (HTR.transpose() * HSCHUR.ldlt().solve(HTR * Ud));
    M = (U.cols() <= U.rows()) ? Eigen::MatrixXd(U.transpose() * AinvU / N) : Eigen::MatrixXd(AinvU * U.transpose() / N);
  }
  Eigen::VectorXd out(2);
  out(0) = M.trace();
  out(1) = M.cwiseProduct(M.transpose()).sum();
  return out;
}

//' Pairwise frequencies
//'
//' @param Y \eqn{n*p} data matrix, with categories coded from zero and `NA` for missing responses.
//' @param CAT Number of categories of each item.
//' @param WEIGHTS Unit weights, or empty for unit weights.
//' @param STRATA Stratum of each row, from zero, or empty for a single stratum.
//'
//' @return The weighted frequencies of the bivariate response patterns: pairs `(k, l)` with `l < k`, ordered by `k`,
//'   and within each pair the categories of `k` then `l`, stratum after stratum.
//' @seealso [lavaanpl_data()]
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd cpp_pairs_freq(const Eigen::Map<Eigen::MatrixXd> Y, const Eigen::Map<Eigen::VectorXd> CAT,
                           const Rcpp::NumericVector WEIGHTS = Rcpp::NumericVector::create(),
                           const Rcpp::NumericVector STRATA = Rcpp::NumericVector::create()){
  const model::Items items(CAT);
  const Eigen::VectorXd w = WEIGHTS.size() > 0 ? Rcpp::as<Eigen::VectorXd>(WEIGHTS) : Eigen::VectorXd::Ones(Y.rows());
  const int ns = STRATA.size() > 0 ? Rcpp::max(STRATA) + 1 : 1;
  Eigen::VectorXd freq = Eigen::VectorXd::Zero(ns * items.ncell);
  for(int r = 0; r < items.npairs; r++){
    const int k = items.pair_k[r], l = items.pair_l[r];
    if(items.cat(k) == 1 || items.cat(l) == 1) continue;
    for(int i = 0; i < Y.rows(); i++){
      if(std::isnan(Y(i, k)) || std::isnan(Y(i, l))) continue;
      const int o = STRATA.size() > 0 ? STRATA[i] * items.ncell : 0;
      freq(o + items.cell_offset[r] + (int)Y(i, k) * items.cat(l) + (int)Y(i, l)) += w(i);
    }
  }
  return freq;
}

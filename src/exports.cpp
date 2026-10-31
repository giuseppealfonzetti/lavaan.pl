#include "workers.h"

template <typename T>
static T get(const Rcpp::List &X, const char *NAME){
  return Rcpp::as<T>(X[NAME]);
}

// [[Rcpp::export]]
Rcpp::List cpp_pl(const Eigen::Map<Eigen::VectorXd> THETA, const Rcpp::List MODEL, const Rcpp::List DATA, const int ORDER = 0){
  const model::Spec S(MODEL, get<Eigen::VectorXd>(DATA, "C_VEC"));
  const model::Setup ST(THETA, S);
  const Eigen::VectorXd freq = get<Eigen::VectorXd>(DATA, "FREQ"), ufreq = get<Eigen::VectorXd>(DATA, "UFREQ");
  std::vector<int> rs(S.L.P);
  std::iota(rs.begin(), rs.end(), 0);
  workers::PL w(ST, freq, rs, ORDER);
  RcppParallel::parallelReduce(0, rs.size(), w);
  Eigen::MatrixXd hd, ud;
  pairs::univariate_contribution(ST.tau, S.L, ufreq, ORDER, w.ll, w.gradient, hd, ud);
  return Rcpp::List::create(Rcpp::Named("pl") = w.ll, Rcpp::Named("gradient") = w.gradient);
}

// [[Rcpp::export]]
Rcpp::List cpp_HJ(const Eigen::Map<Eigen::VectorXd> THETA, const Rcpp::List MODEL, const Rcpp::List DATA, const bool PAIRVAR = false){
  const model::Spec S(MODEL, get<Eigen::VectorXd>(DATA, "C_VEC"));
  const model::Setup ST(THETA, S);
  const Eigen::VectorXd freq = get<Eigen::VectorXd>(DATA, "FREQ"), ufreq = get<Eigen::VectorXd>(DATA, "UFREQ");
  const double n = get<double>(DATA, "N");

  std::vector<Eigen::MatrixXd> score(S.L.P);
  std::vector<std::vector<int>> sidx(S.L.P);
  workers::Hess hw(ST, freq, score, sidx);
  RcppParallel::parallelReduce(0, S.L.P, hw);
  double ll = 0;
  Eigen::MatrixXd uscore = Eigen::MatrixXd::Zero(S.L.nthr, ufreq.size());
  pairs::univariate_contribution(ST.tau, S.L, ufreq, 2, ll, hw.gradient, hw.hessian, uscore);

  const Eigen::MatrixXd y = get<Eigen::MatrixXd>(DATA, "Y"), wm = get<Eigen::MatrixXd>(DATA, "WM");
  const Eigen::VectorXd uw = get<Eigen::VectorXd>(DATA, "UW"), wv = get<Eigen::VectorXd>(DATA, "WV");
  Eigen::MatrixXd cases(S.d, y.rows());
  workers::Cases cw(S.L, y, score, sidx, uscore, uw, cases);
  RcppParallel::parallelFor(0, y.rows(), cw);
  if(std::isfinite(S.rob)){
    Eigen::VectorXd m = Eigen::VectorXd::Zero(S.d);
    for(int r = 0; r < S.L.P; r++){
      for(int c = 0; c < score[r].cols(); c++){
        const double nc = freq(S.L.off[r] + c);
        for(std::size_t a = 0; a < sidx[r].size(); a++) m(sidx[r][a]) += nc * score[r](a, c);
      }
    }
    cases.colwise() -= m / n;
  }
  const Eigen::MatrixXd sc = wm.size() > 0 ? Eigen::MatrixXd(cases * wm) : Eigen::MatrixXd(cases * wv.asDiagonal());

  Eigen::MatrixXd pairvar;
  if(PAIRVAR){
    Eigen::VectorXd ug = Eigen::VectorXd::Zero(S.d);
    Eigen::MatrixXd hd, ud;
    double l0 = 0;
    pairs::univariate_contribution(ST.tau, S.L, ufreq, 1, l0, ug, hd, ud);
    ug /= S.p - 1;
    pairvar = Eigen::MatrixXd::Zero(S.d, S.d);
    for(int r = 0; r < S.L.P; r++){
      const std::vector<int> &I = sidx[r];
      Eigen::VectorXd g = score[r] * freq.segment(S.L.off[r], score[r].cols());
      for(int a = 0; a < S.L.cat(S.L.pk[r]) + S.L.cat(S.L.pl[r]) - 2; a++) g(a) += ug(I[a]);
      for(std::size_t a = 0; a < I.size(); a++){
        for(std::size_t b = 0; b < I.size(); b++) pairvar(I[a], I[b]) += g(a) * g(b);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("H") = hw.hessian / n,
    Rcpp::Named("J") = Eigen::MatrixXd(sc * sc.transpose() / n),
    Rcpp::Named("gradient") = hw.gradient,
    Rcpp::Named("pairvar") = pairvar
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_rho(const Eigen::Map<Eigen::VectorXd> THETA, const Rcpp::List MODEL, const Rcpp::List DATA){
  const model::Spec S(MODEL, get<Eigen::VectorXd>(DATA, "C_VEC"));
  const model::Setup ST(THETA, S);
  Eigen::VectorXd rho(S.L.P);
  Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(S.L.P, S.d);
  for(int r = 0; r < S.L.P; r++){
    const pairs::Pair PR(ST, r, true);
    rho(r) = PR.rho;
    for(int a = 0; a < PR.s; a++) jac(r, PR.idx[a]) = PR.grho(a);
  }
  return Rcpp::List::create(Rcpp::Named("rho") = rho, Rcpp::Named("jacobian") = jac);
}

// [[Rcpp::export]]
Rcpp::List cpp_saturated(const Eigen::Map<Eigen::VectorXd> TAU, const Eigen::Map<Eigen::VectorXd> RHO, const Rcpp::List DATA,
                         const int ORDER, const bool SCORES){
  const model::Layout L(get<Eigen::VectorXd>(DATA, "C_VEC"));
  const Eigen::VectorXd freq = get<Eigen::VectorXd>(DATA, "FREQ"), ufreq = get<Eigen::VectorXd>(DATA, "UFREQ");
  const int order = SCORES ? std::max(ORDER, 1) : ORDER;

  Eigen::MatrixXd patscore = Eigen::MatrixXd::Zero(SCORES ? 2 * L.cat.maxCoeff() - 1 : 0, SCORES ? L.ncell : 0);
  workers::Saturated sat(L, TAU, RHO, freq, order, SCORES, patscore);
  RcppParallel::parallelReduce(0, L.P, sat);
  Eigen::MatrixXd uscore = Eigen::MatrixXd::Zero(SCORES ? L.nthr : 0, SCORES ? ufreq.size() : 0);
  const Eigen::MatrixXd htt_pairs = sat.htt;
  pairs::univariate_contribution(TAU, L, ufreq, order, sat.ll, sat.gradient, sat.htt, uscore);
  sat.hschur += sat.htt - htt_pairs;

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("ll") = sat.ll, Rcpp::Named("gradient") = sat.gradient,
    Rcpp::Named("htt") = sat.htt, Rcpp::Named("htr") = sat.htr, Rcpp::Named("hrr") = sat.hrr,
    Rcpp::Named("hschur") = sat.hschur
  );
  if(SCORES){
    const Eigen::MatrixXd y = get<Eigen::MatrixXd>(DATA, "Y");
    const Eigen::VectorXd uw = get<Eigen::VectorXd>(DATA, "UW");
    Eigen::MatrixXd st = Eigen::MatrixXd::Zero(L.nthr, y.rows()), sr = Eigen::MatrixXd::Zero(L.P, y.rows());
    workers::SaturatedCases cases(L, y, patscore, uscore, uw, st, sr);
    RcppParallel::parallelFor(0, y.rows(), cases);
    out["st"] = st;
    out["sr"] = sr;
  }
  return out;
}

// [[Rcpp::export]]
Eigen::VectorXd cpp_saturated_traces(
    const Eigen::Map<Eigen::MatrixXd> HTT,
    const Eigen::Map<Eigen::MatrixXd> HTR,
    const Eigen::Map<Eigen::VectorXd> D,
    const Eigen::Map<Eigen::MatrixXd> B,
    const Eigen::Map<Eigen::MatrixXd> ST,
    const Eigen::Map<Eigen::MatrixXd> SR,
    const double N,
    const bool FULL = false
){
  const Eigen::VectorXd dinv = D.cwiseInverse();
  Eigen::MatrixXd M;
  if(FULL){
    Eigen::MatrixXd S(ST.rows() + SR.rows(), ST.cols()), X(ST.rows() + SR.rows(), ST.cols());
    S << ST, SR;
    const Eigen::MatrixXd xt = B.ldlt().solve(ST - HTR * (dinv.asDiagonal() * SR));
    X << xt, dinv.asDiagonal() * (SR - HTR.transpose() * xt);
    M = (S.cols() <= S.rows()) ? Eigen::MatrixXd(S.transpose() * X / N) : Eigen::MatrixXd(X * S.transpose() / N);
  }else{
    const Eigen::MatrixXd U = SR - HTR.transpose() * HTT.ldlt().solve(ST);
    const Eigen::MatrixXd Ud = dinv.asDiagonal() * U;
    const Eigen::MatrixXd AinvU = Ud + dinv.asDiagonal() * (HTR.transpose() * B.ldlt().solve(HTR * Ud));
    M = (U.cols() <= U.rows()) ? Eigen::MatrixXd(U.transpose() * AinvU / N) : Eigen::MatrixXd(AinvU * U.transpose() / N);
  }
  Eigen::VectorXd out(2);
  out(0) = M.trace();
  out(1) = M.cwiseProduct(M.transpose()).sum();
  return out;
}

//' Pairwise frequencies
//'
//' @param Y Data matrix of dimension \eqn{n*p}, with categories coded from zero and `NA` for
//'   missing responses.
//' @param C_VEC Number of categories of each item.
//' @param WEIGHTS Optional vector of unit weights; empty for unit weights.
//'
//' @return The (weighted) frequencies of the bivariate response patterns: pairs \eqn{(k, l)},
//'   \eqn{l < k}, ordered by \eqn{k}, and within each pair the categories of \eqn{k} then \eqn{l}.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd pairs_freq(const Eigen::Map<Eigen::MatrixXd> Y, const Eigen::Map<Eigen::VectorXd> C_VEC,
                           const Rcpp::NumericVector WEIGHTS = Rcpp::NumericVector::create()){
  const model::Layout L(C_VEC);
  const Eigen::VectorXd w = WEIGHTS.size() > 0 ? Rcpp::as<Eigen::VectorXd>(WEIGHTS) : Eigen::VectorXd::Ones(Y.rows());
  Eigen::VectorXd freq = Eigen::VectorXd::Zero(L.ncell);
  for(int r = 0; r < L.P; r++){
    const int k = L.pk[r], l = L.pl[r];
    for(int i = 0; i < Y.rows(); i++){
      if(std::isnan(Y(i, k)) || std::isnan(Y(i, l))) continue;
      freq(L.off[r] + (int)Y(i, k) * L.cat(l) + (int)Y(i, l)) += w(i);
    }
  }
  return freq;
}

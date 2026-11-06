#include "workers.h"

void workers::PairLoglik::operator()(std::size_t BEGIN, std::size_t END){
  for(std::size_t i = BEGIN; i < END; i++) pairs::add_pair_gradient(par, data, pair_indices[i], order, ll, gradient);
}

void workers::PairLoglik::join(const PairLoglik &O){
  ll += O.ll;
  gradient += O.gradient;
}

void workers::PairHessian::operator()(std::size_t BEGIN, std::size_t END){
  for(std::size_t r = BEGIN; r < END; r++){
    const pairs::Pair pair(par, r, true);
    const int nthr = pair.nthr, nlocal = pair.nlocal;
    const std::vector<int> I = pair.parameter_indices();

    const bool keep = !cell_scores.empty() && data.tables && !pair.numeric_k && !pair.numeric_l &&
      data.strata.cols() * par.structure.items.ncell * (nthr + nlocal) <= 5e7;
    Eigen::MatrixXd cells(keep ? nthr + nlocal : 0, keep ? data.strata.cols() * pair.cat_k * pair.cat_l : 0);
    Eigen::MatrixXd info = Eigen::MatrixXd::Zero(pair.nreduced, pair.nreduced);
    Eigen::VectorXd score = Eigen::VectorXd::Zero(pair.nreduced);
    double ll = 0;
    pairs::pair_loglik(par.structure.items, pairs::Reduced(par), data, pair, pairs::Order::hessian, ll, score, info, keep ? &cells : nullptr);
    const Eigen::VectorXd local_score = pair.to_local(score);

    Eigen::MatrixXd local_hessian(nthr + nlocal, nthr + nlocal);
    if(pair.standardised){
      local_hessian.topLeftCorner(nthr, nthr) = info.topLeftCorner(nthr, nthr);
      local_hessian.topRightCorner(nthr, nlocal) = info.col(nthr).head(nthr) * pair.d_cov.transpose();
      local_hessian.bottomLeftCorner(nlocal, nthr) = local_hessian.topRightCorner(nthr, nlocal).transpose();
      local_hessian.bottomRightCorner(nlocal, nlocal) = info(nthr, nthr) * pair.d_cov * pair.d_cov.transpose() -
        score(nthr) * pairs::cov_hessian(par, pair);
    }else{
      local_hessian = pair.jacobian.transpose() * info * pair.jacobian - pairs::reduced_curvature(par, pair, score);
    }
    for(int i = 0; i < nthr + nlocal; i++){
      gradient(I[i]) += local_score(i);
      for(int j = 0; j < nthr + nlocal; j++) hessian(I[i], I[j]) += local_hessian(i, j);
    }
    if(!cell_scores.empty()){
      cell_scores[r] = std::move(cells);
      score_index[r] = I;
      pair_gradient[r] = local_score;
    }
  }
}

void workers::PairHessian::join(const PairHessian &O){
  hessian += O.hessian;
  gradient += O.gradient;
}

static void add_univariate_case_scores(const model::Items &ITEMS, const pairs::GroupData &DATA,
                                       const Eigen::Ref<const Eigen::MatrixXd> UNIVARIATE_SCORES,
                                       const Eigen::Ref<const Eigen::VectorXd> UWEIGHTS, const std::size_t BEGIN,
                                       const std::size_t END, Eigen::MatrixXd &SCORES){
  if(UWEIGHTS.size() == 0) return;
  for(std::size_t u = BEGIN; u < END; u++){
    if(UWEIGHTS(u) == 0) continue;
    const int ou = DATA.stratum(u) * ITEMS.cat.sum();
    for(int k = 0, uk = 0; k < ITEMS.p; uk += ITEMS.cat(k), k++){
      if(!std::isnan(DATA.responses(u, k))) SCORES.col(u) += UWEIGHTS(u) * UNIVARIATE_SCORES.col(ou + uk + DATA.responses(u, k));
    }
  }
}

void workers::CaseScores::operator()(std::size_t BEGIN, std::size_t END){
  const model::Items &items = par.structure.items;
  const pairs::Reduced reduced(par);
  pairs::UnitWork W;
  for(int r = 0; r < items.npairs; r++){
    const int k = items.pair_k[r], l = items.pair_l[r], ck = items.cat(k), cl = items.cat(l);
    const std::vector<int> &I = score_index[r];
    if(cell_scores[r].size() > 0){
      for(std::size_t u = BEGIN; u < END; u++){
        if(std::isnan(data.responses(u, k)) || std::isnan(data.responses(u, l))) continue;
        const int c = data.stratum(u) * ck * cl + data.responses(u, k) * cl + data.responses(u, l);
        for(std::size_t a = 0; a < I.size(); a++) scores(I[a], u) += cell_scores[r](a, c);
      }
      continue;
    }
    const pairs::Pair pair(par, r, true);
    pairs::Grid grid(ck, cl, pair.rho, pairs::Order::gradient);
    Eigen::VectorXd score(pair.nreduced);
    for(std::size_t u = BEGIN; u < END; u++){
      if(std::isnan(data.responses(u, k)) || std::isnan(data.responses(u, l))) continue;
      pairs::unit_loglik(reduced, data, pair, grid, u, pairs::Order::gradient, W);
      score.setZero();
      pair.covariate_weights(data.strata.col(data.stratum(u)), data.responses(u, k), data.responses(u, l), W.weights_k, W.weights_l);
      pair.add_score(1, W.score, W.weights_k, W.weights_l, score);
      pair.add_log_scale(1, reduced.scale, score, nullptr);
      const Eigen::VectorXd v = pair.to_local(score);
      for(std::size_t a = 0; a < I.size(); a++) scores(I[a], u) += v(a);
    }
  }
  add_univariate_case_scores(items, data, univariate_scores, uweights, BEGIN, END, scores);
}

static std::vector<int> saturated_index(const pairs::ReducedIndex &INDEX, const pairs::Pair &PAIR){
  std::vector<int> idx(PAIR.nreduced);
  for(int a = 0; a < PAIR.cat_k - 1; a++) idx[a] = PAIR.tau_offset_k + a;
  for(int a = 0; a < PAIR.cat_l - 1; a++) idx[PAIR.nloc_k + a] = PAIR.tau_offset_l + a;
  if(PAIR.numeric_k) idx[0] = INDEX.location(PAIR.k);
  if(PAIR.numeric_l) idx[PAIR.nloc_k] = INDEX.location(PAIR.l);
  idx[PAIR.nloc] = INDEX.rho(PAIR.index);
  const int slope_k = PAIR.nloc + 1, slope_l = slope_k + PAIR.nslope_k;
  for(int c = 0; c < INDEX.nx; c++){
    idx[slope_k + c] = INDEX.slope(c, PAIR.k);
    idx[slope_l + c] = INDEX.slope(c, PAIR.l);
  }
  if(PAIR.numeric_k) idx[slope_k + INDEX.nx] = INDEX.scale(PAIR.k);
  if(PAIR.numeric_l) idx[slope_l + INDEX.nx] = INDEX.scale(PAIR.l);
  return idx;
}

void workers::SaturatedLoglik::operator()(std::size_t BEGIN, std::size_t END){
  const pairs::ReducedIndex index(items, reduced.slope.cols());
  Eigen::VectorXd g;
  Eigen::MatrixXd cell_hessian;
  for(std::size_t r = BEGIN; r < END; r++){
    const pairs::Pair pair(items, r, rho(r), index.nx);
    const int rr = pair.nloc, m = pair.nreduced;
    const std::vector<int> idx = saturated_index(index, pair);
    Eigen::VectorXd score = Eigen::VectorXd::Zero(m);
    Eigen::MatrixXd info = Eigen::MatrixXd::Zero(m, m);
    pairs::pair_loglik(items, reduced, data, pair, order, ll, score, info);
    if(order != pairs::Order::value){
      for(int i = 0; i < m; i++) gradient(idx[i]) += score(i);
    }
    if(order == pairs::Order::hessian){
      const double hr = info(rr, rr);
      for(int i = 0; i < m; i++){
        if(i == rr) continue;
        for(int j = 0; j < m; j++){
          if(j == rr) continue;
          htt(idx[i], idx[j]) += info(i, j);
          hschur(idx[i], idx[j]) += info(i, j) - info(i, rr) * info(j, rr) / hr;
        }
        htr(idx[i], r) += info(i, rr);
      }
      hrr(r) += hr;
    }
    if(scores && data.tables && !pair.numeric_k && !pair.numeric_l){
      for(int t = 0; t < data.strata.cols(); t++){
        const pairs::Grid grid(reduced.tau, pair.tau_offset_k, pair.tau_offset_l, pair.cat_k, pair.cat_l, pair.rho, pairs::Order::gradient,
                               data.strata.col(t).dot(reduced.slope.row(pair.k).transpose()),
                               data.strata.col(t).dot(reduced.slope.row(pair.l).transpose()));
        for(int sk = 0; sk < pair.cat_k; sk++){
          for(int sl = 0; sl < pair.cat_l; sl++){
            const double pe = grid.cell_probability(sk, sl, g, cell_hessian);
            pattern_scores.col(t * items.ncell + pair.cell_offset + sk * pair.cat_l + sl).head(pair.nthr + 1) = g / pe;
          }
        }
      }
    }
  }
}

void workers::SaturatedLoglik::join(const SaturatedLoglik &O){
  ll += O.ll;
  gradient += O.gradient;
  htt += O.htt;
  htr += O.htr;
  hrr += O.hrr;
  hschur += O.hschur;
}

void workers::SaturatedCaseScores::operator()(std::size_t BEGIN, std::size_t END){
  const pairs::ReducedIndex index(items, reduced.slope.cols());
  const int nitem = score_tau.rows();
  pairs::UnitWork W;
  for(int r = 0; r < items.npairs; r++){
    const pairs::Pair pair(items, r, rho(r), index.nx);
    const int k = pair.k, l = pair.l;
    const bool tables = data.tables && !pair.numeric_k && !pair.numeric_l;
    pairs::Grid grid(pair.cat_k, pair.cat_l, pair.rho, pairs::Order::gradient);
    const std::vector<int> idx = saturated_index(index, pair);
    Eigen::VectorXd score(pair.nreduced);
    for(std::size_t u = BEGIN; u < END; u++){
      if(std::isnan(data.responses(u, k)) || std::isnan(data.responses(u, l))) continue;
      const int t = data.stratum(u);
      if(tables){
        const int sk = data.responses(u, k), sl = data.responses(u, l);
        W.score = pattern_scores.col(t * items.ncell + pair.cell_offset + sk * pair.cat_l + sl).head(pair.nloc + 1);
      }else{
        pairs::unit_loglik(reduced, data, pair, grid, u, pairs::Order::gradient, W);
      }
      score.setZero();
      pair.covariate_weights(data.strata.col(t), data.responses(u, k), data.responses(u, l), W.weights_k, W.weights_l);
      pair.add_score(1, W.score, W.weights_k, W.weights_l, score);
      pair.add_log_scale(1, reduced.scale, score, nullptr);
      for(int i = 0; i < pair.nreduced; i++){
        if(idx[i] < nitem) score_tau(idx[i], u) += score(i);
      }
      score_rho(r, u) = score(pair.nloc);
    }
  }
  add_univariate_case_scores(items, data, univariate_scores, uweights, BEGIN, END, score_tau);
}

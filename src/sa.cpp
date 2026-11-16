#include <random>
#include <chrono>
#include "workers.h"

// [[Rcpp::export]]
Rcpp::List cpp_sa(
    const Eigen::Map<Eigen::VectorXd> START,
    const Eigen::Map<Eigen::MatrixXd> K,
    const Eigen::Map<Eigen::VectorXd> C,
    const Rcpp::List MODEL,
    const Rcpp::List DATA,
    const int PAIRS_PER_ITERATION,
    const double ALPHA,
    const int ADAM_EPOCHS,
    const int PRECONDITION_EPOCHS,
    const int MAX_PRECONDITION_EPOCHS,
    const int AVERAGE_EPOCHS,
    const double TOL,
    const double DECAY,
    const double AVG_FROM,
    const bool HESSIAN_GAIN,
    const int LEVEL_EPOCHS,
    const int HALVINGS,
    const int SEED,
    const bool VERBOSE = false
){
  typedef std::chrono::steady_clock clk;
  auto secs = [](clk::time_point a){ return std::chrono::duration<double>(clk::now() - a).count(); };
  using pairs::Order;

  const model::Groups groups(MODEL, DATA);
  const int G = groups.ngroups;
  std::vector<pairs::GroupData> data;
  for(int g = 0; g < G; g++){
    data.emplace_back(Rcpp::as<Rcpp::List>(DATA[g]));
  }
  const double n = groups.n;
  const int P = groups.npairs, batch = std::min(PAIRS_PER_ITERATION, P), epoch_length = (P + batch - 1) / batch, df = START.size();
  const bool constrained = K.size() > 0;
  const Eigen::SparseMatrix<double> K_sparse = K.sparseView();
  const double beta1 = .9, beta2 = .999;
  const double alpha_start = ALPHA / std::sqrt(static_cast<double>(epoch_length));
  double alpha = alpha_start;

  auto expand = [&](const Eigen::VectorXd &X){ return constrained ? Eigen::VectorXd(K_sparse * X + C) : X; };

  std::vector<std::vector<int>> rescov_pairs(G), all_pairs(G);
  for(int g = 0; g < G; g++){
    const model::Structure &S = groups.structure[g];
    all_pairs[g].resize(S.items.npairs);
    std::iota(all_pairs[g].begin(), all_pairs[g].end(), 0);
    for(int r = 0; r < S.items.npairs && S.theta.size() > 0; r++){
      const int k = S.items.pair_k[r], l = S.items.pair_l[r];
      if(S.rescov_index(k, l) >= 0 || S.theta(k, l) != 0) rescov_pairs[g].push_back(r);
    }
  }

  auto feasible = [&](const std::vector<model::Matrices> &PARS){
    for(int g = 0; g < G; g++){
      const model::Matrices &par = PARS[g];
      const model::Structure &S = groups.structure[g];
      if(!par.lambda_ib_inv.allFinite()) return false;
      for(int k = 0; k < S.p; k++){
        const int tk = S.items.tau_offset(k);
        for(int c = 1; c < S.items.cat(k) - 1; c++) if(par.tau(tk + c) <= par.tau(tk + c - 1)) return false;
      }
      if(S.latcov == model::LatCov::direct && par.psi.llt().info() != Eigen::Success) return false;
      const Eigen::MatrixXd LP = par.lambda_ib_inv * par.psi;
      const Eigen::VectorXd sc = S.standardised ? Eigen::VectorXd::Ones(S.p) : par.scale;
      const Eigen::VectorXd v = par.common_var.cwiseProduct(sc.cwiseAbs2());
      auto ok = [&](const int k, const int l){ return std::abs(sc(k) * sc(l) * (LP.row(k).dot(par.lambda_ib_inv.row(l)) + par.rescov(k, l))) < 1; };
      for(const int r : rescov_pairs[g]) if(!ok(S.items.pair_k[r], S.items.pair_l[r])) return false;
      if(v.maxCoeff() < 1) continue;
      for(int k = 0; k < S.p; k++){
        if(v(k) < 1) continue;
        for(int l = 0; l < S.p; l++) if(l != k && v(k) * v(l) >= 1 && !ok(k, l)) return false;
      }
    }
    return true;
  };

  std::mt19937 rng(SEED);
  std::vector<int> pool(P), drawn(batch);
  std::vector<std::vector<int>> drawn_by_group(G);
  std::iota(pool.begin(), pool.end(), 0);
  bool shuffled = false;
  int cursor = P;
  auto draw_pairs = [&](){
    for(int g = 0; g < G; g++) drawn_by_group[g].clear();
    if(shuffled && cursor + batch > P){
      for(int i = 0; i < P - 1; i++){
        std::uniform_int_distribution<int> u(i, P - 1);
        std::swap(pool[i], pool[u(rng)]);
      }
      cursor = 0;
    }
    for(int i = 0; i < batch; i++){
      if(shuffled){
        drawn[i] = pool[cursor++];
      }else{
        std::uniform_int_distribution<int> u(i, P - 1);
        std::swap(pool[i], pool[u(rng)]);
        drawn[i] = pool[i];
      }
      int g = G - 1;
      while(groups.pair_offset[g] > drawn[i]) g--;
      drawn_by_group[g].push_back(drawn[i] - groups.pair_offset[g]);
    }
  };

  auto stochastic_gradient = [&](const std::vector<model::Matrices> &PARS, Eigen::VectorXd &GRADIENT){
    double ll = 0;
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(groups.d);
    for(int g = 0; g < G; g++){
      const model::Structure &S = groups.structure[g];
      workers::PairLoglik w(PARS[g], data[g], drawn_by_group[g], Order::gradient);
      RcppParallel::parallelReduce(0, drawn_by_group[g].size(), w, 4);
      w.ll *= static_cast<double>(P) / batch;
      w.gradient *= static_cast<double>(P) / batch;
      if(data[g].ufreq.size() > 0){
        Eigen::VectorXd item_weights = Eigen::VectorXd::Zero(S.p);
        for(const int r : drawn_by_group[g]){
          item_weights(S.items.pair_k[r]) += 1;
          item_weights(S.items.pair_l[r]) += 1;
        }
        item_weights *= static_cast<double>(P) / (batch * (S.p - 1));
        Eigen::MatrixXd hd, ud;
        pairs::add_univariate(PARS[g], data[g], Order::gradient, w.ll, w.gradient, hd, ud, item_weights);
      }
      ll += w.ll;
      grad.segment(groups.parameter_offset[g], S.d) = w.gradient;
    }
    if(!std::isfinite(ll) || !grad.allFinite()) return false;
    GRADIENT = constrained ? Eigen::VectorXd(-K_sparse.transpose() * grad / n) : Eigen::VectorXd(-grad / n);
    return true;
  };

  auto objective = [&](const Eigen::VectorXd &X){
    const std::vector<model::Matrices> pars = groups.matrices(expand(X));
    double ll = 0;
    for(int g = 0; g < G; g++){
      if(groups.structure[g].latcov == model::LatCov::direct && pars[g].psi.llt().info() != Eigen::Success) return R_PosInf;
      workers::PairLoglik w(pars[g], data[g], all_pairs[g], Order::value);
      RcppParallel::parallelReduce(0, all_pairs[g].size(), w);
      Eigen::VectorXd gd;
      Eigen::MatrixXd hd, ud;
      pairs::add_univariate(pars[g], data[g], Order::value, w.ll, gd, hd, ud);
      ll += w.ll;
    }
    return std::isfinite(ll) ? -ll / n : R_PosInf;
  };

  bool gain_pd = false;
  auto hessian_gain = [&](const Eigen::VectorXd &X, Eigen::MatrixXd &GAIN, const double DAMP){
    const std::vector<model::Matrices> pars = groups.matrices(expand(X));
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(groups.d, groups.d);
    for(int g = 0; g < G; g++){
      std::vector<Eigen::MatrixXd> cell_scores;
      std::vector<std::vector<int>> score_index;
      std::vector<Eigen::VectorXd> pair_gradient;
      workers::PairHessian hw(pars[g], data[g], cell_scores, score_index, pair_gradient);
      RcppParallel::parallelReduce(0, groups.structure[g].items.npairs, hw);
      double ll = 0;
      Eigen::MatrixXd ud;
      pairs::add_univariate(pars[g], data[g], Order::hessian, ll, hw.gradient, hw.hessian, ud);
      H.block(groups.parameter_offset[g], groups.parameter_offset[g], groups.structure[g].d, groups.structure[g].d) = hw.hessian / n;
    }
    if(constrained) H = K_sparse.transpose() * H * K_sparse;
    if(!H.allFinite()) return false;
    const Eigen::VectorXd s = H.diagonal().cwiseAbs().cwiseMax(1e-12).cwiseSqrt().cwiseInverse();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s.asDiagonal() * (0.5 * (H + H.transpose())) * s.asDiagonal());
    gain_pd = es.eigenvalues().minCoeff() > 1e-6;
    Eigen::VectorXd ev = es.eigenvalues().cwiseAbs();
    ev = ev.array() + DAMP * ev.maxCoeff();
    GAIN = s.asDiagonal() * es.eigenvectors() * ev.cwiseInverse().asDiagonal() * es.eigenvectors().transpose() * s.asDiagonal();
    return GAIN.allFinite();
  };

  auto apply_gain = [&](const Eigen::MatrixXd &GAIN, const Eigen::VectorXd &GRADIENT){
    Eigen::VectorXd out = Eigen::VectorXd::Zero(df);
    for(int j = 0; j < df; j++) if(GRADIENT(j) != 0) out += GRADIENT(j) * GAIN.col(j);
    return out;
  };

  Eigen::VectorXd x = START, g;
  std::unique_ptr<std::vector<model::Matrices>> current;
  auto move_to = [&](const Eigen::VectorXd &X){
    std::unique_ptr<std::vector<model::Matrices>> pars(new std::vector<model::Matrices>(groups.matrices(expand(X))));
    if(!feasible(*pars)) return false;
    x = X;
    current = std::move(pars);
    return true;
  };

  int rejections = 0, skipped = 0, iter = 0;
  auto sa_step = [&](auto &&DIRECTION){
    iter++;
    draw_pairs();
    if(!stochastic_gradient(*current, g)){
      skipped++;
      return;
    }
    Eigen::VectorXd step = DIRECTION(g);
    for(int k = 0; k < 20; k++){
      if(move_to(x - step)) return;
      step /= 2;
      rejections++;
    }
    skipped++;
  };

  std::vector<Eigen::VectorXd> path_vector;
  std::vector<double> path_nll;
  std::vector<int> path_iter, path_phase;
  auto record_path = [&](const Eigen::VectorXd &X, const double F, const int PHASE){
    path_vector.push_back(expand(X));
    if(PHASE < 3) path_nll.push_back(F);
    path_iter.push_back(iter);
    path_phase.push_back(PHASE);
    if(VERBOSE){
      Rcpp::Rcout << "Phase " << PHASE << " | iteration " << iter;
      if(PHASE < 3) Rcpp::Rcout << " | objective " << F;
      Rcpp::Rcout << "\n";
    }
  };

  double t_adam = 0, t_precondition = 0, t_hess = 0, t_avg = 0, t_mon = 0, prev_f = 0, best_f = 0;
  Eigen::VectorXd best, epoch_mean;
  int calm_epochs = 0;
  bool plateau = false;
  auto burn_epoch = [&](auto &&DIRECTION, const int PHASE){
    Eigen::VectorXd epoch_average = Eigen::VectorXd::Zero(df);
    for(int i = 0; i < epoch_length; i++){
      sa_step(DIRECTION);
      epoch_average += (x - epoch_average) / (i + 1);
      Rcpp::checkUserInterrupt();
    }
    clk::time_point tm = clk::now();
    epoch_mean = epoch_average;
    const double f = objective(epoch_average);
    t_mon += secs(tm);
    record_path(epoch_average, f, PHASE);
    if(f < best_f){
      best = epoch_average;
      best_f = f;
    }
    calm_epochs = std::isfinite(f) && std::abs(n * (prev_f - f)) < TOL * df ? calm_epochs + 1 : 0;
    prev_f = f;
    return calm_epochs >= 2;
  };

  Eigen::VectorXd adam_m = Eigen::VectorXd::Zero(df), adam_v = Eigen::VectorXd::Zero(df);
  int adam_t = 0, adam_epochs = 0, precondition_epochs = 0, n_hess = 0;
  auto adam = [&](const Eigen::VectorXd &GRADIENT){
    adam_t++;
    adam_m = beta1 * adam_m + (1 - beta1) * GRADIENT;
    adam_v = beta2 * adam_v + (1 - beta2) * GRADIENT.cwiseAbs2();
    return Eigen::VectorXd(alpha * (adam_m / (1 - std::pow(beta1, adam_t))).array() /
                           ((adam_v / (1 - std::pow(beta2, adam_t))).array().sqrt() + 1e-8));
  };
  clk::time_point t0 = clk::now();
  if(!move_to(START)) Rcpp::stop("The starting values are outside the parameter space.");
  clk::time_point tm = clk::now();
  prev_f = objective(x);
  t_mon += secs(tm);
  if(!std::isfinite(prev_f)) Rcpp::stop("The starting values are outside the parameter space.");
  record_path(x, prev_f, 0);
  best = x;
  best_f = prev_f;
  const bool schedule = LEVEL_EPOCHS > 0;
  Eigen::VectorXd level = Eigen::VectorXd::Zero(df);
  int halvings = 0;
  for(int e = 0; e < ADAM_EPOCHS && (schedule || !plateau); e++){
    plateau = burn_epoch(adam, 1);
    adam_epochs++;
    if(!schedule) continue;
    level += epoch_mean / LEVEL_EPOCHS;
    if(adam_epochs % LEVEL_EPOCHS != 0) continue;
    clk::time_point tc = clk::now();
    const double f = objective(level);
    t_mon += secs(tc);
    if(f < best_f){
      best = level;
      best_f = f;
    }
    move_to(best);
    level.setZero();
    if(halvings < HALVINGS){
      alpha /= 2;
      halvings++;
    }
  }
  if(schedule) plateau = false;
  move_to(best);
  prev_f = best_f;
  t_adam = secs(t0) - t_mon;

  Eigen::MatrixXd gain;
  Eigen::VectorXd diagonal_gain;
  bool hessian = HESSIAN_GAIN;
  double gain_scale = 1.0 / epoch_length;
  t0 = clk::now();
  const double mon0 = t_mon;
  double damp = 1e-3;
  if(hessian){
    calm_epochs = 0;
    plateau = false;
    shuffled = true;
    while(true){
      clk::time_point th = clk::now();
      hessian = hessian_gain(x, gain, damp);
      t_hess += secs(th);
      n_hess++;
      if(!hessian) break;
      const bool enough = plateau || precondition_epochs >= PRECONDITION_EPOCHS;
      if((enough && gain_pd) || precondition_epochs >= std::max(PRECONDITION_EPOCHS, MAX_PRECONDITION_EPOCHS)) break;
      const double before = best_f;
      plateau = burn_epoch([&](const Eigen::VectorXd &GRADIENT){ return Eigen::VectorXd(gain_scale * apply_gain(gain, GRADIENT)); }, 2);
      damp = best_f < before ? std::max(damp / 10, 1e-6) : std::min(damp * 10, 1.0);
      precondition_epochs++;
      move_to(best);
    }
    shuffled = false;
  }
  if(!hessian){
    diagonal_gain = alpha_start * (1 - beta1) * ((adam_v / (1 - std::pow(beta2, std::max(adam_t, 1)))).array().sqrt() + 1e-8).inverse();
    gain_scale = 1;
  }
  t_precondition = secs(t0) - (t_mon - mon0);
  const Eigen::VectorXd x_burn = x;

  const int skipped_burn = skipped, rejections_burn = rejections;
  Eigen::VectorXd avg = x;
  const int T2 = AVERAGE_EPOCHS * epoch_length, T0 = static_cast<int>(AVG_FROM * T2);
  t0 = clk::now();
  for(int s = 1; s <= T2; s++){
    const double step_size = gain_scale * std::pow(1.0 + static_cast<double>(s) / epoch_length, -DECAY);
    if(hessian){
      sa_step([&](const Eigen::VectorXd &GRADIENT){ return Eigen::VectorXd(step_size * apply_gain(gain, GRADIENT)); });
    }else{
      sa_step([&](const Eigen::VectorXd &GRADIENT){ return Eigen::VectorXd(step_size * diagonal_gain.cwiseProduct(GRADIENT)); });
    }
    if(s <= T0){
      avg = x;
    }else{
      avg += (x - avg) / (s - T0 + 1);
    }
    if(s % epoch_length == 0){
      record_path(avg, 0, 3);
      Rcpp::checkUserInterrupt();
    }
  }
  t_avg = secs(t0);

  return Rcpp::List::create(
    Rcpp::Named("parameter_vector") = expand(avg),
    Rcpp::Named("last_iterate") = expand(x),
    Rcpp::Named("burn_in_iterate") = expand(x_burn),
    Rcpp::Named("epoch_length") = epoch_length,
    Rcpp::Named("pairs_per_iteration") = batch,
    Rcpp::Named("adam_epochs") = adam_epochs,
    Rcpp::Named("precondition_epochs") = precondition_epochs,
    Rcpp::Named("plateau") = plateau,
    Rcpp::Named("halvings") = halvings,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("averaged") = T2 - T0,
    Rcpp::Named("gain") = hessian ? "hessian" : "diagonal",
    Rcpp::Named("n_hessian") = n_hess,
    Rcpp::Named("invh") = hessian && T2 == 0 ? gain : Eigen::MatrixXd(),
    Rcpp::Named("rejections") = rejections,
    Rcpp::Named("skipped") = skipped,
    Rcpp::Named("rejections_averaging") = rejections - rejections_burn,
    Rcpp::Named("skipped_averaging") = skipped - skipped_burn,
    Rcpp::Named("path") = Rcpp::List::create(
      Rcpp::Named("iter") = path_iter, Rcpp::Named("phase") = path_phase,
      Rcpp::Named("parameter_vector") = path_vector, Rcpp::Named("nll") = path_nll),
    Rcpp::Named("time") = Rcpp::NumericVector::create(
      Rcpp::Named("adam") = t_adam, Rcpp::Named("precondition") = t_precondition, Rcpp::Named("hessian") = t_hess,
      Rcpp::Named("averaging") = t_avg, Rcpp::Named("monitor") = t_mon)
  );
}

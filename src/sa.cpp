#include <random>
#include <chrono>
#include "workers.h"

// [[Rcpp::export]]
Rcpp::List cpp_sa(
    const Eigen::Map<Eigen::MatrixXd> STARTS,
    const Eigen::Map<Eigen::MatrixXd> K,
    const Eigen::Map<Eigen::VectorXd> C,
    const Rcpp::List MODEL,
    const Rcpp::List DATA,
    const Eigen::Map<Eigen::VectorXd> VALFREQ,
    const Eigen::Map<Eigen::VectorXd> VALUFREQ,
    const int PAIRS_PER_ITERATION,
    const double ALPHA,
    const int ADAM_EPOCHS,
    const int BURN_EPOCHS,
    const int EPOCHS,
    const double TOL,
    const double DECAY,
    const double AVG_FROM,
    const int GAIN,
    const int SEED,
    const bool ADAM_ONLY = false,
    const bool VERBOSE = false
){
  typedef std::chrono::steady_clock clk;
  auto secs = [](clk::time_point a){ return std::chrono::duration<double>(clk::now() - a).count(); };

  const model::Spec S(MODEL, Rcpp::as<Eigen::VectorXd>(DATA["C_VEC"]));
  const Eigen::VectorXd freq = Rcpp::as<Eigen::VectorXd>(DATA["FREQ"]), ufreq = Rcpp::as<Eigen::VectorXd>(DATA["UFREQ"]);
  const double n = Rcpp::as<double>(DATA["N"]);
  const int P = S.L.P, nu = std::min(PAIRS_PER_ITERATION, P), E = (P + nu - 1) / nu, df = STARTS.rows();
  const bool constrained = K.size() > 0;
  const Eigen::SparseMatrix<double> Ks = K.sparseView();
  const double b1 = .9, b2 = .999, alpha = ALPHA / std::sqrt(static_cast<double>(E));

  auto full = [&](const Eigen::VectorXd &X){ return constrained ? Eigen::VectorXd(Ks * X + C) : X; };

  std::vector<int> rpairs;
  for(int r = 0; r < P && S.rescov.size() > 0; r++){
    const int k = S.L.pk[r], l = S.L.pl[r];
    if(S.resid(k, l) >= 0 || S.rescov(k, l) != 0) rpairs.push_back(r);
  }

  auto feasible = [&](const model::Setup &ST){
    if(!ST.Lt.allFinite()) return false;
    for(int k = 0; k < S.p; k++){
      for(int c = 1; c < S.L.cat(k) - 1; c++) if(ST.tau(S.L.tk(k) + c) <= ST.tau(S.L.tk(k) + c - 1)) return false;
    }
    Eigen::MatrixXd Psi;
    if(S.direct){
      if(ST.Sig.llt().info() != Eigen::Success) return false;
      Psi = ST.Sig;
    }else{
      Psi = ST.sd.asDiagonal() * ST.U.transpose() * ST.U * ST.sd.asDiagonal();
    }
    const Eigen::MatrixXd LP = ST.Lt * Psi;
    const Eigen::VectorXd v = LP.cwiseProduct(ST.Lt).rowwise().sum();
    auto ok = [&](const int k, const int l){ return std::abs(LP.row(k).dot(ST.Lt.row(l)) + ST.rc(k, l)) < 1; };
    for(const int r : rpairs) if(!ok(S.L.pk[r], S.L.pl[r])) return false;
    if(v.maxCoeff() < 1) return true;
    for(int k = 0; k < S.p; k++){
      if(v(k) < 1) continue;
      for(int l = 0; l < S.p; l++) if(l != k && v(k) * v(l) >= 1 && !ok(k, l)) return false;
    }
    return true;
  };

  std::mt19937 rng(SEED);
  std::vector<int> pool(P), drawn(nu), all(P);
  std::iota(pool.begin(), pool.end(), 0);
  std::iota(all.begin(), all.end(), 0);
  auto draw = [&](){
    for(int i = 0; i < nu; i++){
      std::uniform_int_distribution<int> u(i, P - 1);
      std::swap(pool[i], pool[u(rng)]);
      drawn[i] = pool[i];
    }
  };

  Eigen::VectorXd itemw = Eigen::VectorXd::Zero(S.p);
  auto stoch = [&](const model::Setup &ST, Eigen::VectorXd &G){
    workers::PL w(ST, freq, drawn, 1);
    RcppParallel::parallelReduce(0, nu, w, 4);
    double ll = w.ll * P / nu;
    Eigen::VectorXd grad = w.gradient * P / nu;
    Eigen::MatrixXd hd, ud;
    if(ufreq.size() > 0){
      itemw.setZero();
      for(const int r : drawn){
        itemw(S.L.pk[r]) += 1;
        itemw(S.L.pl[r]) += 1;
      }
      itemw *= static_cast<double>(P) / (nu * (S.p - 1));
      pairs::univariate_contribution(ST.tau, S.L, ufreq, 1, ll, grad, hd, ud, itemw);
    }
    if(!std::isfinite(ll) || !grad.allFinite()) return false;
    G = constrained ? Eigen::VectorXd(-Ks.transpose() * grad / n) : Eigen::VectorXd(-grad / n);
    return true;
  };

  auto objective = [&](const Eigen::VectorXd &X){
    const model::Setup ST(full(X), S);
    if(S.direct && ST.Sig.llt().info() != Eigen::Success) return R_PosInf;
    workers::PL w(ST, VALFREQ, all, 0);
    RcppParallel::parallelReduce(0, P, w);
    Eigen::VectorXd gd;
    Eigen::MatrixXd hd, ud;
    pairs::univariate_contribution(ST.tau, S.L, VALUFREQ, 0, w.ll, gd, hd, ud);
    return std::isfinite(w.ll) ? -w.ll / n : R_PosInf;
  };

  auto hessian_gain = [&](const Eigen::VectorXd &X, Eigen::MatrixXd &G, const double DAMP){
    const model::Setup ST(full(X), S);
    std::vector<Eigen::MatrixXd> sc;
    std::vector<std::vector<int>> si;
    workers::Hess hw(ST, freq, sc, si);
    RcppParallel::parallelReduce(0, P, hw);
    double ll = 0;
    Eigen::MatrixXd ud;
    pairs::univariate_contribution(ST.tau, S.L, ufreq, 2, ll, hw.gradient, hw.hessian, ud);
    Eigen::MatrixXd H = hw.hessian / n;
    if(constrained) H = Ks.transpose() * H * Ks;
    if(!H.allFinite()) return false;
    const Eigen::VectorXd s = H.diagonal().cwiseAbs().cwiseMax(1e-12).cwiseSqrt().cwiseInverse();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s.asDiagonal() * (0.5 * (H + H.transpose())) * s.asDiagonal());
    Eigen::VectorXd ev = es.eigenvalues().cwiseAbs();
    ev = ev.cwiseMax(DAMP * ev.maxCoeff());
    G = s.asDiagonal() * es.eigenvectors() * ev.cwiseInverse().asDiagonal() * es.eigenvectors().transpose() * s.asDiagonal();
    return G.allFinite();
  };

  auto gain_times = [&](const Eigen::MatrixXd &Gm, const Eigen::VectorXd &G){
    Eigen::VectorXd out = Eigen::VectorXd::Zero(df);
    for(int j = 0; j < df; j++) if(G(j) != 0) out += G(j) * Gm.col(j);
    return out;
  };

  Eigen::VectorXd x = STARTS.col(0), g;
  std::unique_ptr<model::Setup> cur;
  auto setx = [&](const Eigen::VectorXd &X){
    std::unique_ptr<model::Setup> st(new model::Setup(full(X), S));
    if(!feasible(*st)) return false;
    x = X;
    cur = std::move(st);
    return true;
  };

  int rejections = 0, skipped = 0, iter = 0;
  auto update = [&](auto &&DIR){
    iter++;
    draw();
    if(!stoch(*cur, g)){
      skipped++;
      return;
    }
    Eigen::VectorXd step = DIR(g);
    for(int k = 0; k < 20; k++){
      if(setx(x - step)) return;
      step /= 2;
      rejections++;
    }
    skipped++;
  };

  std::vector<Eigen::VectorXd> path_theta;
  std::vector<double> path_nll;
  std::vector<int> path_iter, path_phase;
  auto record = [&](const Eigen::VectorXd &X, const double F, const int PHASE){
    path_theta.push_back(full(X));
    if(PHASE < 3) path_nll.push_back(F);
    path_iter.push_back(iter);
    path_phase.push_back(PHASE);
    if(VERBOSE){
      Rcpp::Rcout << "Phase " << PHASE << " | iteration " << iter;
      if(PHASE < 3) Rcpp::Rcout << " | objective " << F;
      Rcpp::Rcout << "\n";
    }
  };

  double t_adam = 0, t_newton = 0, t_hess = 0, t_avg = 0, t_mon = 0, prev_f = 0, best_f = 0;
  Eigen::VectorXd best;
  int calm = 0;
  bool plateau = false;
  auto burn_epoch = [&](auto &&DIR, const int PHASE){
    Eigen::VectorXd ea = Eigen::VectorXd::Zero(df);
    for(int i = 0; i < E; i++){
      update(DIR);
      ea += (x - ea) / (i + 1);
      Rcpp::checkUserInterrupt();
    }
    clk::time_point tm = clk::now();
    const double f = objective(ea);
    t_mon += secs(tm);
    record(ea, f, PHASE);
    if(f < best_f){
      best = ea;
      best_f = f;
    }
    calm = std::isfinite(f) && std::abs(n * (prev_f - f)) < TOL * df ? calm + 1 : 0;
    prev_f = f;
    return calm >= 2;
  };

  Eigen::VectorXd m = Eigen::VectorXd::Zero(df), v = Eigen::VectorXd::Zero(df);
  int ta = 0, adam_epochs = 0, newton_epochs = 0, n_hess = 0;
  auto adam = [&](const Eigen::VectorXd &G){
    ta++;
    m = b1 * m + (1 - b1) * G;
    v = b2 * v + (1 - b2) * G.cwiseAbs2();
    return Eigen::VectorXd(alpha * (m / (1 - std::pow(b1, ta))).array() / ((v / (1 - std::pow(b2, ta))).array().sqrt() + 1e-8));
  };
  const int max_adam = ADAM_ONLY || GAIN == 1 ? ADAM_EPOCHS : ADAM_EPOCHS + BURN_EPOCHS;
  Eigen::VectorXd xw, vw;
  double fw = R_PosInf;
  int taw = 0;
  bool pw = false;
  clk::time_point t0 = clk::now();
  for(int j = 0; j < STARTS.cols(); j++){
    if(!setx(STARTS.col(j))) Rcpp::stop("The starting values are outside the parameter space.");
    clk::time_point tm = clk::now();
    prev_f = objective(x);
    t_mon += secs(tm);
    if(!std::isfinite(prev_f)) Rcpp::stop("The starting values are outside the parameter space.");
    record(x, prev_f, 0);
    best = x;
    best_f = prev_f;
    m.setZero();
    v.setZero();
    ta = calm = 0;
    plateau = false;
    for(int e = 0; e < max_adam && !plateau; e++){
      plateau = burn_epoch(adam, 1);
      adam_epochs++;
    }
    if(best_f < fw){
      xw = best;
      fw = best_f;
      vw = v;
      taw = ta;
      pw = plateau;
    }
  }
  setx(xw);
  best = xw;
  prev_f = best_f = fw;
  v = vw;
  ta = taw;
  plateau = pw;
  t_adam = secs(t0) - t_mon;

  Eigen::MatrixXd Gm;
  Eigen::VectorXd Gd;
  bool hessian = GAIN == 1 && !ADAM_ONLY;
  double kappa = 1.0 / E;
  t0 = clk::now();
  const double mon0 = t_mon;
  double damp = 1e-3;
  if(hessian){
    calm = 0;
    plateau = false;
    while(newton_epochs < BURN_EPOCHS && !plateau){
      clk::time_point th = clk::now();
      if(!hessian_gain(x, Gm, damp)){
        hessian = false;
        break;
      }
      t_hess += secs(th);
      n_hess++;
      const double before = best_f;
      plateau = burn_epoch([&](const Eigen::VectorXd &G){ return Eigen::VectorXd(kappa * gain_times(Gm, G)); }, 2);
      damp = best_f < before ? std::max(damp / 10, 1e-6) : std::min(damp * 10, 1.0);
      newton_epochs++;
      setx(best);
    }
    if(hessian){
      clk::time_point th = clk::now();
      hessian = hessian_gain(x, Gm, damp);
      t_hess += secs(th);
      n_hess++;
    }
  }
  if(!hessian){
    Gd = alpha * (1 - b1) * ((v / (1 - std::pow(b2, std::max(ta, 1)))).array().sqrt() + 1e-8).inverse();
    kappa = 1;
  }
  t_newton = secs(t0) - (t_mon - mon0);
  const Eigen::VectorXd x_burn = x;

  const int skipped_burn = skipped, rejections_burn = rejections;
  Eigen::VectorXd avg = x;
  const int T2 = ADAM_ONLY ? 0 : EPOCHS * E, T0 = static_cast<int>(AVG_FROM * T2);
  t0 = clk::now();
  for(int s = 1; s <= T2; s++){
    const double gam = kappa * std::pow(1.0 + static_cast<double>(s) / E, -DECAY);
    if(hessian){
      update([&](const Eigen::VectorXd &G){ return Eigen::VectorXd(gam * gain_times(Gm, G)); });
    }else{
      update([&](const Eigen::VectorXd &G){ return Eigen::VectorXd(gam * Gd.cwiseProduct(G)); });
    }
    if(s <= T0){
      avg = x;
    }else{
      avg += (x - avg) / (s - T0 + 1);
    }
    if(s % E == 0){
      record(avg, 0, 3);
      Rcpp::checkUserInterrupt();
    }
  }
  t_avg = secs(t0);

  return Rcpp::List::create(
    Rcpp::Named("theta") = full(ADAM_ONLY ? x_burn : avg),
    Rcpp::Named("theta_last") = full(x),
    Rcpp::Named("theta_burn") = full(x_burn),
    Rcpp::Named("epoch_length") = E,
    Rcpp::Named("adam_epochs") = adam_epochs,
    Rcpp::Named("newton_epochs") = newton_epochs,
    Rcpp::Named("plateau") = plateau,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("averaged") = T2 - T0,
    Rcpp::Named("gain") = hessian ? "hessian" : "diagonal",
    Rcpp::Named("n_hessian") = n_hess,
    Rcpp::Named("rejections") = rejections,
    Rcpp::Named("skipped") = skipped,
    Rcpp::Named("rejections_averaging") = rejections - rejections_burn,
    Rcpp::Named("skipped_averaging") = skipped - skipped_burn,
    Rcpp::Named("path") = Rcpp::List::create(
      Rcpp::Named("iter") = path_iter, Rcpp::Named("phase") = path_phase,
      Rcpp::Named("theta") = path_theta, Rcpp::Named("nll") = path_nll),
    Rcpp::Named("time") = Rcpp::NumericVector::create(
      Rcpp::Named("adam") = t_adam, Rcpp::Named("newton") = t_newton, Rcpp::Named("hessian") = t_hess,
      Rcpp::Named("averaging") = t_avg, Rcpp::Named("monitor") = t_mon)
  );
}

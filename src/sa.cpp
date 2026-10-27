#include <random>
#include "workers.h"

// Rescale the loadings (and free log-sds) of the items whose implied variance of the
// underlying variable exceeds one (Alfonzetti et al., 2025, Appendix B.2)
static void proj2(const model::Spec &S, Eigen::VectorXd &THETA, bool &CHECKEVENT){
  const model::Setup ST(THETA, S);
  Eigen::MatrixXd Lam = ST.Lam;
  Eigen::VectorXd sd = ST.sd;
  const Eigen::MatrixXd Sigma = sd.asDiagonal() * ST.U.transpose() * ST.U * sd.asDiagonal();

  for(int k = 0; k < S.p; k++){
    const Eigen::VectorXd lambdak = Lam.row(k);
    const double vk = lambdak.transpose() * Sigma * lambdak;
    if(vk <= 1) continue;
    CHECKEVENT = true;
    for(int h = 0; h < S.q; h++){
      if(!std::isfinite(S.constrmat(k, h))){
        if(!std::isfinite(S.constrlogsd(h))){
          Lam(k, h) /= std::sqrt(std::sqrt(vk));
          sd(h) /= std::sqrt(std::sqrt(vk));
        }else{
          Lam(k, h) /= std::sqrt(vk);
        }
      }else if(S.constrmat(k, h) != 0.0){
        sd(h) /= std::sqrt(vk);
      }
    }
  }
  if(!CHECKEVENT) return;
  for(int v = 0; v < S.nvar; v++) THETA(S.ov + v) = std::log(sd(S.sfree[v]));
  for(int j = 0; j < S.p; j++){
    for(std::size_t a = 0; a < S.lidx[j].size(); a++) THETA(S.lidx[j][a]) = Lam(j, S.lfac[j][a]);
  }
}

// [[Rcpp::export]]
Eigen::VectorXd cpp_sa_proj(const Eigen::Map<Eigen::VectorXd> THETA, const Rcpp::List MODEL, const Rcpp::List DATA){
  const model::Spec S(MODEL, Rcpp::as<Eigen::VectorXd>(DATA["C_VEC"]));
  Eigen::VectorXd theta = THETA;
  bool flag = false;
  proj2(S, theta, flag);
  return theta;
}

// [[Rcpp::export]]
Rcpp::List cpp_plSA2(
    const Eigen::Map<Eigen::VectorXd> THETA_INIT,
    const Rcpp::List MODEL,
    const Rcpp::List DATA,
    const Eigen::Map<Eigen::VectorXd> VALFREQ,
    const Eigen::Map<Eigen::VectorXd> VALUFREQ,
    const int PAIRS_PER_ITERATION,
    const double STEP0,
    const double STEP1,
    const double STEP2,
    const double STEP3,
    const int BURNE,
    const int MAXE,
    const int EPOCHS_PER_CHECK,
    const int MAX_TOL_COUNTER,
    const double TOL_BURN,
    const double TOL_END,
    const int CHECK_TOL,
    const int SEED,
    const bool VERBOSE,
    bool VERBOSE_ITER = false
){
  const model::Spec S(MODEL, Rcpp::as<Eigen::VectorXd>(DATA["C_VEC"]));
  const Eigen::VectorXd freq = Rcpp::as<Eigen::VectorXd>(DATA["FREQ"]), ufreq = Rcpp::as<Eigen::VectorXd>(DATA["UFREQ"]);
  const int n = Rcpp::as<int>(DATA["N"]);
  Eigen::MatrixXd hdummy, udummy;

  if(!VERBOSE) VERBOSE_ITER = false;

  const int d = THETA_INIT.size();
  const int pairs = S.L.P;
  if(d != S.d) Rcpp::stop("check theta dimensions");
  if(PAIRS_PER_ITERATION > pairs) Rcpp::stop("too many pairs per iterations");
  const int NTHR = S.L.nthr, NLOAD = S.nload, NCORR = S.ncorr, NVAR = S.nvar;

  int burnt              = 0;
  int burne              = BURNE;
  const double prob      = static_cast<double>(PAIRS_PER_ITERATION) / static_cast<double>(pairs);
  int tol_counter        = 0;
  int tol_counter_burn   = 0;
  double stepsize        = STEP0;

  std::vector<int> full_pool(pairs);
  std::iota(std::begin(full_pool), std::end(full_pool), 0);

  Eigen::VectorXd theta      = THETA_INIT;
  Eigen::VectorXd avtheta    = THETA_INIT;
  Eigen::VectorXd thetapdiff = Eigen::VectorXd::Zero(d);

  std::vector<Eigen::VectorXd> path_theta;
  std::vector<Eigen::VectorXd> path_avtheta;
  std::vector<Eigen::VectorXd> path_thetapdiff;
  std::vector<double>          path_nll;
  std::vector<int>             path_iters;
  std::vector<int>             post_index;

  // Full negative pairwise log-likelihood per unit on the validation frequencies
  auto full_nll = [&](const Eigen::VectorXd &TH){
    const model::Setup ST(TH, S);
    workers::PL w(ST, VALFREQ, full_pool, 0);
    RcppParallel::parallelReduce(0, pairs, w);
    double ull = 0;
    Eigen::VectorXd ugrad;
    pairs::univariate_contribution(ST.tau, S.L, VALUFREQ, 0, ull, ugrad, hdummy, udummy);
    return -(w.ll + ull) / n;
  };

  const int maxtpe = pairs / PAIRS_PER_ITERATION;
  double nll = full_nll(theta);
  if(VERBOSE) Rcpp::Rcout << "Initial full npll: " << nll << " | iterations per epoch:" << maxtpe << "\n";

  const double scale = prob / static_cast<double>(n);

  // Diagnostics
  bool neg_pdiff        = false;
  bool proj_after_burn  = false;
  bool finite_gr        = true;
  bool convergence_burn = true;
  int  convergence_full = 0;
  int last_iter = 0;
  double pdiff = 0;

  Eigen::VectorXd prev_theta = theta;
  path_iters.push_back(0);
  path_theta.push_back(theta);
  path_avtheta.push_back(avtheta);
  path_thetapdiff.push_back(thetapdiff);
  path_nll.push_back(nll);

  int t = 1;
  for(int epoch = 0; epoch < MAXE; epoch++){
    Rcpp::checkUserInterrupt();
    Eigen::VectorXd prev_epoch_theta = theta;
    Eigen::VectorXd prev_epoch_avtheta = avtheta;

    // Sampling scheme: one random permutation of the pairs per epoch
    std::mt19937 randomizer(SEED + epoch);
    std::shuffle(full_pool.begin(), full_pool.end(), randomizer);
    int idx_1pair = 0;

    for(int te = 0; te < maxtpe; te++){
      Rcpp::checkUserInterrupt();

      // Stochastic gradient on the sampled pairs
      std::vector<int> iter_chosen_pairs(full_pool.begin() + idx_1pair, full_pool.begin() + idx_1pair + PAIRS_PER_ITERATION);
      idx_1pair += PAIRS_PER_ITERATION;
      const model::Setup ST(theta, S);
      workers::PL w(ST, freq, iter_chosen_pairs, 1);
      RcppParallel::parallelReduce(0, iter_chosen_pairs.size(), w);
      Eigen::VectorXd iter_gradient = -scale * w.gradient;
      if(ufreq.size() > 0){
        double ull = 0;
        Eigen::VectorXd ugrad = Eigen::VectorXd::Zero(d);
        pairs::univariate_contribution(ST.tau, S.L, ufreq, 1, ull, ugrad, hdummy, udummy);
        iter_gradient -= scale * prob * ugrad;
      }

      // Parameters update
      bool update = true;
      int attempt_counter = 0;
      double step_length = stepsize;
      while(update && attempt_counter < 10){
        step_length = stepsize * STEP1 * pow(1 + STEP2 * stepsize * t, -STEP3);
        iter_gradient.segment(NTHR, NLOAD) *= 5;
        if(NCORR > 0) iter_gradient.segment(NTHR + NLOAD, NCORR) *= 100;
        if(NVAR > 0)  iter_gradient.segment(NTHR + NLOAD + NCORR, NVAR) *= 5;
        if(std::isfinite(iter_gradient.norm())){
          update = false;
          prev_theta = theta;
          theta -= step_length * iter_gradient;
        }else{
          stepsize *= .75;
          theta = prev_theta;
          finite_gr = false;
        }
        attempt_counter++;
      }

      // Post-update projection of the loadings, if needed
      bool checkevent = false;
      proj2(S, theta, checkevent);
      if(checkevent) post_index.push_back(t);
      if(checkevent && epoch > burne) proj_after_burn = true;

      thetapdiff = (theta.array() - prev_theta.array()) / prev_theta.array();
      if(VERBOSE_ITER) Rcpp::Rcout << " - Iter: " << t << " | step: " << step_length << " | mean abs theta pdiff: " <<
        thetapdiff.array().abs().mean() << " | gradient norm: " << iter_gradient.norm() << "\n";

      if(epoch <= burne){
        avtheta = theta;
      }else{
        avtheta = ((t - burnt) * avtheta + theta) / (t - burnt + 1);
      }
      t++;
    }

    // Full negative pairwise log-likelihood, path storage and convergence monitoring
    const double prev_nll = nll;
    if(epoch % EPOCHS_PER_CHECK == 0){
      nll = full_nll(theta);
      path_iters.push_back(t);
      path_theta.push_back(theta);
      path_avtheta.push_back(avtheta);
      path_thetapdiff.push_back(thetapdiff);
      path_nll.push_back(nll);
      pdiff = - (nll - prev_nll) / (prev_nll * maxtpe);
      if(pdiff < TOL_END){tol_counter++;}else{tol_counter = 0;}
      if(pdiff < TOL_BURN){tol_counter_burn++;}else{tol_counter_burn = 0;}
      if((pdiff < 0) && (epoch <= burne)) neg_pdiff = true;
    }

    if(!std::isfinite(nll)){
      if(VERBOSE) Rcpp::Rcout << "Exiting because of non-finite nll\n";
      last_iter        = t;
      convergence_full = -1;
      break;
    }

    if(VERBOSE && (epoch % EPOCHS_PER_CHECK == 0)){
      const Eigen::VectorXd epdiff = (theta.array() - prev_epoch_theta.array()) / prev_epoch_theta.array();
      Rcpp::Rcout << "End of Epoch: " << epoch << " | full npll: " << nll << " | av nll pdiff from prev epoch: " << pdiff <<
        " | tc:" << tol_counter << "| mean abs theta pdiff from prev epoch: " << epdiff.array().abs().mean() << "\n";
    }

    // Burn-in: ends at BURNE, or earlier when the tolerance is reached
    if(epoch == burne){
      if(VERBOSE) Rcpp::Rcout << "Burn-in ended\n";
      convergence_burn = false;
      burnt = t;
    }else if(CHECK_TOL && (tol_counter_burn >= MAX_TOL_COUNTER) && (epoch <= burne) && (epoch % EPOCHS_PER_CHECK == 0)){
      if(VERBOSE) Rcpp::Rcout << "Burn-in stopped\n";
      burne = epoch;
      tol_counter_burn = 0;
      convergence_burn = true;
      burnt = t;
    }

    // Exit
    if(CHECK_TOL && (tol_counter >= MAX_TOL_COUNTER) && (epoch > burne) && (epoch % EPOCHS_PER_CHECK == 0)){
      if(VERBOSE) Rcpp::Rcout << "Converged\n";
      convergence_full = 1;
      last_iter = t;
      if(pdiff < 0){
        if(VERBOSE) Rcpp::Rcout << "Reporting parameter values from epoch " << epoch - 1 << "\n";
        theta = prev_epoch_theta;
        avtheta = prev_epoch_avtheta;
      }
      break;
    }else if(epoch == (MAXE - 1)){
      if(VERBOSE) Rcpp::Rcout << "Ended\n";
      last_iter = t;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("burnt") = burnt,
    Rcpp::Named("scale") = scale,
    Rcpp::Named("nll") = nll,
    Rcpp::Named("path_theta") = path_theta,
    Rcpp::Named("path_thetapdiff") = path_thetapdiff,
    Rcpp::Named("path_avtheta") = path_avtheta,
    Rcpp::Named("path_iters") = path_iters,
    Rcpp::Named("path_nll") = path_nll,
    Rcpp::Named("post_index") = post_index,
    Rcpp::Named("last_iter") = last_iter,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("avtheta") = avtheta,
    Rcpp::Named("convergence_full") = convergence_full,
    Rcpp::Named("convergence_burn") = convergence_burn,
    Rcpp::Named("proj_after_burn") = proj_after_burn,
    Rcpp::Named("finite_gr") = finite_gr,
    Rcpp::Named("neg_pdiff") = neg_pdiff
  );
}

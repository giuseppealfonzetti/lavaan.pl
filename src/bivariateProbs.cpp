#include "bivariateProbs.h"




/*  PI (the probability of a bivariate response pattern) */
// Compute specific pi_sksl
double biprobs::compute_pi(
    const Eigen::Ref<const Eigen::VectorXd> c_vec,
    const Eigen::Ref<const Eigen::VectorXd> pi_thresholds,
    double rho_kl,

    const unsigned int k,
    const unsigned int l,
    const unsigned int sk,
    const unsigned int sl
){
  unsigned int ck = c_vec(k);
  unsigned int cl = c_vec(l);

  // read pi related thresholds
  double t_sk = pi_thresholds(0);
  double t_sl = pi_thresholds(1);
  double t_sk_prev = pi_thresholds(2);
  double t_sl_prev = pi_thresholds(3);

  // Phi(t_sk, t_sl; rho_kl)
  double cum1;
  if ((sk == (ck-1)) && (sl == (cl-1))) {
    cum1 = 1;
  } else if(sk == (ck-1)){
    cum1 = R::pnorm(t_sl, 0, 1, 1, 0);
  } else if(sl == (cl-1)){
    cum1 = R::pnorm(t_sk, 0, 1, 1, 0);
  } else {
    cum1 = binorm::pbvnorm( t_sk, t_sl, rho_kl);
  }

  // Phi(t_sk, t_sl-1; rho_kl)
  double cum2;
  if(sl == 0){
    cum2 = 0;
  } else {
    cum2 = binorm::pbvnorm( t_sk, t_sl_prev, rho_kl);
  }
  // Phi(t_sk-1, t_sl; rho_kl)
  double cum3;
  if(sk == 0){
    cum3 = 0;
  } else {
    cum3 = binorm::pbvnorm( t_sk_prev, t_sl, rho_kl);
  }
  // Phi(t_sk-1, t_sl-1; rho_kl)
  double cum4;
  if ((sl == 0) || (sk == 0)) {
    cum4 = 0;
  } else {
    cum4 = binorm::pbvnorm( t_sk_prev, t_sl_prev, rho_kl);
  }

  double pi_sksl = cum1 - cum2 - cum3 + cum4;

  return pi_sksl;
}


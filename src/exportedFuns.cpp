#include "exportedFuns.h"

Rcpp::List cpp_compute_pair_ext(
     Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
     Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
     Eigen::Map<Eigen::VectorXd> C_VEC,
     Eigen::VectorXd THETA,
     const int NTHR,
     const int NLOAD,
     const int NCORR,
     const int NVAR,
     const unsigned int K,
     const unsigned int L,
     Eigen::MatrixXd PAIRS_TABLE,
     const unsigned int GRADFLAG,
     const int CORRFLAG,
     Rcpp::Nullable<Rcpp::NumericMatrix> REGMAT,
     Rcpp::Nullable<Rcpp::NumericMatrix> RESCOV
 ){
   const int d = THETA.size();

   double ll = 0;
   Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
   const pairs::Setup S(THETA, CONSTRMAT, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, CORRFLAG,
                        params::nullable_mat(REGMAT), params::nullable_mat(RESCOV));
   if(d != NTHR + S.m) Rcpp::stop("Check theta size");
   pairs::pair_contribution(S, C_VEC, PAIRS_TABLE, K, L, GRADFLAG, ll, gradient);

   Rcpp::List output =
     Rcpp::List::create(
       Rcpp::Named("nll") = -ll,
       Rcpp::Named("ngradient") = -gradient
     );
   return(output);

 }

Eigen::VectorXd cpp_loadings_theta2vec(
    Eigen::Map<Eigen::VectorXd> THETA,
    const int NTHR,
    const int NLOAD
){
  return params::loadings::theta2vec(THETA, NTHR, NLOAD);
}

Eigen::MatrixXd cpp_loadings_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int NTHR,
    const int NLOAD
){
  return params::loadings::theta2mat(THETA, CONSTRMAT, NTHR, NLOAD);
}


Eigen::MatrixXd cpp_loadings_mat2vec(
    Eigen::Map<Eigen::MatrixXd> LOADINGS,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int  NLOAD){
  return params::loadings::mat2vec(LOADINGS, CONSTRMAT, NLOAD);
}

Eigen::MatrixXd cpp_latvar_vec2cmat(Eigen::Map<Eigen::VectorXd> VEC,
                             const int NCORR,
                             const int Q){
  return params::latvar::vec2cmat(VEC, NCORR, Q);
}

Eigen::MatrixXd cpp_latvar_theta2cmat(
    Eigen::Map<Eigen::VectorXd> THETA,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q){

  Eigen::VectorXd vec = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR);
  return params::latvar::vec2cmat(vec, NCORR, Q);
}

Eigen::MatrixXd cpp_latvar_theta2dmat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q){


  return params::latvar::theta2dmat(THETA,CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q);
}

Eigen::MatrixXd cpp_latvar_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q,
    const int CORRFLAG){
  if(CORRFLAG == 2){
    const Eigen::VectorXd sd = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q).diagonal();
    Eigen::MatrixXd out = sd.array().square().matrix().asDiagonal();
    for(int c = 1, a = 0; c < Q; c++){
      for(int r = 0; r < c; r++, a++){
        out(r, c) = THETA(NTHR + NLOAD + a);
        out(c, r) = THETA(NTHR + NLOAD + a);
      }
    }
    return out;
  }
  return params::latvar::theta2mat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q);
}

Eigen::VectorXd cpp_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S,
                                   Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
                                   const int NCORR,
                                   const int NVAR,
                                   const int CORRFLAG){
  if(CORRFLAG == 2){
    const int q = S.cols();
    Eigen::VectorXd out(NCORR + NVAR);
    for(int c = 1, a = 0; c < q; c++){
      for(int r = 0; r < c; r++, a++) out(a) = S(r, c);
    }
    if(NVAR > 0) out.tail(NVAR) = params::latvar::edvec2dvec(params::latvar::mat2edvec(S), CONSTRLOGSD, q, NVAR);
    return out;
  }
  return params::latvar::mat2vec(S, CONSTRLOGSD, NCORR, NVAR);
}

Eigen::VectorXd cpp_sa_proj(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    const int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR
){
  bool flag = false;
  Eigen::VectorXd theta = THETA;
  sa::proj2(CONSTRMAT,
           CONSTRLOGSD,
           C_VEC,
           CORRFLAG,
           NTHR,
           NLOAD,
           NCORR,
           NVAR,
           theta,
           flag);
  return theta;
}

Rcpp::List cpp_sample_estimators_HJ(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    Eigen::Map<Eigen::MatrixXd> DATA,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    int N,
    int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    Rcpp::NumericVector UFREQ,
    Rcpp::NumericVector UW,
    Rcpp::NumericVector GW,
    Rcpp::NumericVector JW,
    Rcpp::NumericVector CL,
    Rcpp::Nullable<Rcpp::NumericMatrix> REGMAT,
    Rcpp::Nullable<Rcpp::NumericMatrix> RESCOV
){
  // Rcpp::Rcout << "ef|theta:\n";
  // Rcpp::Rcout << THETA.transpose()<<"\n";
  Rcpp::List output = fullPairwise::multiThreadSamplEstHJ(
    THETA,
    FREQ,
    DATA,
    C_VEC,
    CONSTRMAT,
    CONSTRLOGSD,
    N,
    CORRFLAG,
    NTHR,
    NLOAD,
    NCORR,
    NVAR,
    Eigen::Map<Eigen::VectorXd>(UFREQ.begin(), UFREQ.size()),
    Eigen::Map<Eigen::VectorXd>(UW.begin(), UW.size()),
    Eigen::Map<Eigen::VectorXd>(GW.begin(), GW.size()),
    Eigen::Map<Eigen::VectorXd>(JW.begin(), JW.size()),
    Eigen::Map<Eigen::VectorXd>(CL.begin(), CL.size()),
    params::nullable_mat(REGMAT),
    params::nullable_mat(RESCOV)
  );

  return output;


}

Rcpp::List cpp_saturated(
    Eigen::Map<Eigen::VectorXd> TAU,
    Eigen::Map<Eigen::VectorXd> RHO,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    Eigen::Map<Eigen::MatrixXd> DATA,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    const int ORDER,
    const bool SCORES,
    Rcpp::NumericVector UFREQ,
    Rcpp::NumericVector UW
){
  return fullPairwise::multiThreadSaturated(TAU, RHO, FREQ, DATA, C_VEC, ORDER, SCORES,
                                            Eigen::Map<Eigen::VectorXd>(UFREQ.begin(), UFREQ.size()),
                                            Eigen::Map<Eigen::VectorXd>(UW.begin(), UW.size()));
}

Eigen::VectorXd cpp_saturated_traces(
    Eigen::Map<Eigen::MatrixXd> HTT,
    Eigen::Map<Eigen::MatrixXd> HTR,
    Eigen::Map<Eigen::VectorXd> D,
    Eigen::Map<Eigen::MatrixXd> ST,
    Eigen::Map<Eigen::MatrixXd> SR,
    Eigen::Map<Eigen::VectorXd> W,
    const double N,
    const bool FULL
){
  const Eigen::VectorXd dinv = D.cwiseInverse();
  const Eigen::MatrixXd B = HTT - HTR * dinv.asDiagonal() * HTR.transpose();
  Eigen::MatrixXd M;
  if(FULL){
    const Eigen::MatrixXd st = ST * W.cwiseSqrt().asDiagonal(), sr = SR * W.cwiseSqrt().asDiagonal();
    Eigen::MatrixXd S(st.rows() + sr.rows(), st.cols()), X(st.rows() + sr.rows(), st.cols());
    S << st, sr;
    const Eigen::MatrixXd xt = B.ldlt().solve(st - HTR * (dinv.asDiagonal() * sr));
    X << xt, dinv.asDiagonal() * (sr - HTR.transpose() * xt);
    M = (S.cols() <= S.rows()) ? Eigen::MatrixXd(S.transpose() * X / N) : Eigen::MatrixXd(X * S.transpose() / N);
  }else{
    Eigen::MatrixXd U = SR - HTR.transpose() * HTT.ldlt().solve(ST);
    U *= W.cwiseSqrt().asDiagonal();
    const Eigen::MatrixXd Ud = dinv.asDiagonal() * U;
    const Eigen::MatrixXd AinvU = Ud + dinv.asDiagonal() * (HTR.transpose() * B.ldlt().solve(HTR * Ud));
    M = (U.cols() <= U.rows()) ? Eigen::MatrixXd(U.transpose() * AinvU / N) : Eigen::MatrixXd(AinvU * U.transpose() / N);
  }

  Eigen::VectorXd out(2);
  out(0) = M.trace();
  out(1) = M.cwiseProduct(M.transpose()).sum();
  return out;
}

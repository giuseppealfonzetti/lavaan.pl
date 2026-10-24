#include "exportedFuns.h"

Rcpp::List cpp_compute_pair_ext(
     Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
     Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
     const std::vector<std::vector<std::vector<double>>> LLC,
     Eigen::Map<Eigen::VectorXd> C_VEC,
     Eigen::VectorXd THETA,
     const int NTHR,
     const int NLOAD,
     const int NCORR,
     const int NVAR,
     const unsigned int K,
     const unsigned int L,
     Eigen::MatrixXd PAIRS_TABLE,
     const unsigned int GRADFLAG
 ){
   const int d = NTHR+NLOAD+NCORR+NVAR;

   if(THETA.size()!=(NTHR+NLOAD+NCORR+NVAR))Rcpp::stop("Check theta size");

   double ll = 0;
   Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
   const pairs::Setup S(THETA, CONSTRMAT, CONSTRLOGSD, LLC, NTHR, NLOAD, NCORR, NVAR);
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
    const std::vector<std::vector<std::vector<double>>> LLC,
    const int NTHR,
    const int NLOAD
){
  return params::loadings::theta2mat(THETA, CONSTRMAT, LLC, NTHR, NLOAD);
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
    const int Q){
  return params::latvar::theta2mat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q);
}

Eigen::VectorXd cpp_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S,
                                   Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
                                   const int NCORR,
                                   const int NVAR){
  return params::latvar::mat2vec(S, CONSTRLOGSD, NCORR, NVAR);
}

Eigen::VectorXd cpp_sa_proj(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const std::vector<std::vector<std::vector<double>>> LLC,
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
           LLC,
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
    const std::vector<std::vector<std::vector<double>>> LLC,
    int N,
    int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR
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
    LLC,
    N,
    CORRFLAG,
    NTHR,
    NLOAD,
    NCORR,
    NVAR
  );

  return output;


}

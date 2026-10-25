#ifndef exportedFuns_H
#define exportedFuns_H

#include "pairs.h"
#include "optimisationUtils.h"
#include "fullPairwise.h"


// [[Rcpp::export]]
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
 );

// [[Rcpp::export]]
Eigen::VectorXd cpp_loadings_theta2vec(
    Eigen::Map<Eigen::VectorXd> THETA,
    const int NTHR,
    const int NLOAD
);

// [[Rcpp::export]]
Eigen::MatrixXd cpp_loadings_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const std::vector<std::vector<std::vector<double>>> LLC,
    const int NTHR,
    const int NLOAD
);


// [[Rcpp::export]]
Eigen::MatrixXd cpp_loadings_mat2vec(
    Eigen::Map<Eigen::MatrixXd> LOADINGS,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int  NLOAD);

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_vec2cmat(Eigen::Map<Eigen::VectorXd> VEC,
                             const int NCORR,
                             const int Q);

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_theta2cmat(
    Eigen::Map<Eigen::VectorXd> THETA,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q);

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_theta2dmat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q);

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q);

// [[Rcpp::export]]
Eigen::VectorXd cpp_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S,
                                   Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
                                   const int NCORR,
                                   const int NVAR);

// [[Rcpp::export]]
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
);

// [[Rcpp::export]]
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
    const int NVAR,
    Rcpp::NumericVector UFREQ = Rcpp::NumericVector::create(),
    Rcpp::NumericVector UW = Rcpp::NumericVector::create()
);

// [[Rcpp::export]]
Rcpp::List cpp_saturated(
    Eigen::Map<Eigen::VectorXd> TAU,
    Eigen::Map<Eigen::VectorXd> RHO,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    Eigen::Map<Eigen::MatrixXd> DATA,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    const int ORDER,
    const bool SCORES,
    Rcpp::NumericVector UFREQ = Rcpp::NumericVector::create(),
    Rcpp::NumericVector UW = Rcpp::NumericVector::create()
);

// [[Rcpp::export]]
Eigen::VectorXd cpp_saturated_traces(
    Eigen::Map<Eigen::MatrixXd> HTT,
    Eigen::Map<Eigen::MatrixXd> HTR,
    Eigen::Map<Eigen::VectorXd> D,
    Eigen::Map<Eigen::MatrixXd> ST,
    Eigen::Map<Eigen::MatrixXd> SR,
    Eigen::Map<Eigen::VectorXd> W,
    const double N
);


#endif

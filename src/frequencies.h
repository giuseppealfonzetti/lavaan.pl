#ifndef frequencies_H
#define frequencies_H
#include <RcppEigen.h>

//' Compute pairwise frequencies
//'
//' @param Y Integer matrix of dimension \eqn{n*p}, where \eqn{n} is the sample size
//' and \eqn{p} is the number of items considered. Categories must be coded starting from zero.
//' For example, an item with three categories can only accept values contained in
//' \eqn{\{0, 1, 2\}}.
//' @param C_VEC Integer vector indicating how many possible categories are associated to
//' each item in 'Y'.
//' @param WEIGHTS Optional vector of unit weights. Empty for unit weights.
//'
//' @return
//' It returns a 5-rows matrix with each combination of items and categories as columns.
//' Row0: item k, Row1: item l, Row2; category item k, Row3: category item l, Row4: freq
//' It is computed just once, before the optimization of the complete pairwise
//'
//' @export
// [[Rcpp::export]]
 Eigen::MatrixXd pairs_freq(
     Eigen::Map<Eigen::MatrixXd> Y,
     Eigen::Map<Eigen::VectorXd> C_VEC,
     Rcpp::NumericVector WEIGHTS = Rcpp::NumericVector::create()
 );

#endif

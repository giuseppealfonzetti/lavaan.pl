#ifndef bivariate_normal_H
#define bivariate_normal_H
#include <Rcpp.h>


/* BIVARIATE NORMAL */
// Code adapted from pbv package https://github.com/cran/pbv/
// Needed to switch from Rcpp NumericVector to std::vector to avoid
// otherwise unsolved memory problems during code parallelisation.
namespace binorm{


  double dbvnorm(double A, double B, double RHO);

  //pbv_rcpp_pbvnorm0 (Drezner & Wesolowksy, 1990, JCSC)
  double pbvnorm(double A, double B, double RHO);


}

#endif

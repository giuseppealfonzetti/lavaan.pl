#ifndef bivariateNormal_H
#define bivariateNormal_H
const double pi = 3.1415926535897;
#include <Rcpp.h>


/* BIVARIATE NORMAL */
// Code adapted from pbv package https://github.com/cran/pbv/
// Needed to switch from Rcpp NumericVector to std::vector to avoid
// otherwise unsolved memory problems during code parallelisation.
namespace binorm{


  double dbvnorm( double X, double Y, double RHO, bool USE_LOG);

  //pbv_rcpp_pbvnorm0 (Drezner & Wesolowksy, 1990, JCSC)
  double pbvnorm( double H1, double HK, double R);


}






#endif

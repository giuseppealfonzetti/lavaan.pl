#include "bivariate_normal.h"

double binorm::dbvnorm(double A, double B, double RHO){
  double r2 = 1-RHO*RHO;
  double z = A*A - 2*RHO*A*B + B*B;
  return std::exp(- z / r2 / 2.0) / (2*M_PI) / std::sqrt(r2);
}

//pbv_rcpp_pbvnorm0 (Drezner & Wesolowksy, 1990, JCSC)
double binorm::pbvnorm(double A, double B, double RHO){
  constexpr int NX = 5;
  // data
  static constexpr double X[NX] = {.04691008, .23076534, .5, .76923466, .95308992};
  static constexpr double W[NX] = {.018854042, .038088059, .0452707394, .038088059, .018854042};
  // declarations
  double bv = 0;
  double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
  double cor_max = 0.7;
  double bv_fac1 = 0.13298076;
  double bv_fac2 = 0.053051647;
  // computation
  double h2 = B;
  double h12 = (A*A+h2*h2)/2;
  double r_abs = std::abs(RHO);
  if (r_abs > cor_max){
    r2 = 1.0 - RHO*RHO;
    r3 = std::sqrt(r2);
    if (RHO<0){
      h2 = -h2;
    }
    h3 = A*h2;
    h7 = std::exp( -h3 / 2.0);
    if ( r_abs < 1){
      h6 = std::abs(A-h2);
      h5 = h6*h6 / 2.0;
      h6 = h6 / r3;
      aa = 0.5 - h3 / 8.0;
      ab = 3.0 - 2.0 * aa * h5;
      bv = bv_fac1*h6*ab*(1-R::pnorm(h6, 0, 1, 1, 0))-std::exp(-h5/r2)*(ab + aa*r2)*bv_fac2;
      for (int ii=0; ii<NX; ii++){
        r1 = r3*X[ii];
        rr = r1*r1;
        r2 = std::sqrt( 1.0 - rr);
        bv += - W[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);
      }
    }
    h11 = std::min(A,h2);
    bv = bv*r3*h7 + R::pnorm(h11, 0, 1, 1, 0);
    if (RHO < 0){
      bv = R::pnorm(A, 0, 1, 1, 0) - bv;
    }

  } else {
    h3=A*h2;
    for (int ii=0; ii<NX; ii++){
      r1 = RHO*X[ii];
      rr2 = 1.0 - r1*r1;
      bv += W[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
    }
    bv = R::pnorm(A, 0, 1, 1, 0)*R::pnorm(h2, 0, 1, 1, 0) + RHO*bv;
  }
  return bv;
}

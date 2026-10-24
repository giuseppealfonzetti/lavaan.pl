#include "bivariateNormal.h"

double binorm::dbvnorm( double X, double Y, double RHO, bool USE_LOG){
  double pi2 = 2*pi;
  double r2 = 1-RHO*RHO;
  double r3 = std::sqrt(r2);
  double z = X*X - 2*RHO*X*Y + Y*Y;
  z = - z / r2 / 2.0;
  if(!USE_LOG){
    z = std::exp(z) / pi2 / r3;
  } else {
    z += - std::log(r3*pi2);
  }

  return z;
}

//pbv_rcpp_pbvnorm0 (Drezner & Wesolowksy, 1990, JCSC)
double binorm::pbvnorm( double H1, double HK, double R){
  int NX=5;
  std::vector<double> X(NX);
  std::vector<double> W(NX);
  // data
  X[0]=.04691008;
  X[1]=.23076534;
  X[2]=.5;
  X[3]=.76923466;
  X[4]=.95308992;
  W[0]=.018854042;
  W[1]=.038088059;
  W[2]=.0452707394;
  W[3]=.038088059;
  W[4]=.018854042;
  // declarations
  double bv = 0;
  double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
  double cor_max = 0.7;
  double bv_fac1 = 0.13298076;
  double bv_fac2 = 0.053051647;
  // computation
  double h2 = HK;
  double h12 = (H1*H1+h2*h2)/2;
  double r_abs = std::abs(R);
  if (r_abs > cor_max){
    r2 = 1.0 - R*R;
    r3 = std::sqrt(r2);
    if (R<0){
      h2 = -h2;
    }
    h3 = H1*h2;
    h7 = std::exp( -h3 / 2.0);
    if ( r_abs < 1){
      h6 = std::abs(H1-h2);
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
    h11 = std::min(H1,h2);
    bv = bv*r3*h7 + R::pnorm(h11, 0, 1, 1, 0);
    if (R < 0){
      bv = R::pnorm(H1, 0, 1, 1, 0) - bv;
    }

  } else {
    h3=H1*h2;
    for (int ii=0; ii<NX; ii++){
      r1 = R*X[ii];
      rr2 = 1.0 - r1*r1;
      bv += W[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
    }
    bv = R::pnorm(H1, 0, 1, 1, 0)*R::pnorm(h2, 0, 1, 1, 0) + R*bv;
  }
  return bv;
}


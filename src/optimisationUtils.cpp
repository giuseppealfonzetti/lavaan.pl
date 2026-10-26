#include "optimisationUtils.h"
void sa::proj2(
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    const int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    Eigen::VectorXd &THETA,
    bool &CHECKEVENT){
  //dimensions
  const int p = CONSTRMAT.rows(); // number of items
  const int q = CONSTRMAT.cols(); // number of latents

  // rearrange parameters
  Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, CONSTRMAT, NTHR, NLOAD);
  Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
  Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);


  // Normalisation if needed
  for(unsigned int k = 0; k < p; k++){
    const Eigen::VectorXd lambdak = Lam.row(k);
    const double vk = lambdak.transpose()*Sigma_u*lambdak;
    if(vk>1){
      CHECKEVENT = true;

      for(int h=0; h<q; h++){

        if(!std::isfinite(CONSTRMAT(k,h))){
          // Free loading


          if(!std::isfinite(CONSTRLOGSD(h))){
            // Free loading and free variance
            Lam(k,h)/=sqrt(sqrt(vk));
            Du(h,h)/=sqrt(sqrt(vk));

          }else{

            // Free loading but fixed variance
            Lam(k,h)/= sqrt(vk);

          }
        }else{
          // Fixed loading

          if(CONSTRMAT(k,h)!=0.0){
            // Non-zero loading
            Du(h,h)/=sqrt(vk);
          }

        }
      }
    }
  }



  if(CHECKEVENT){

    int idx = 0;
    if(NVAR>0){
      Eigen::VectorXd norm_logsd(NVAR);
      for(unsigned int h = 0; h < q; h++){
        if(!std::isfinite(CONSTRLOGSD(h))){
          norm_logsd(idx) = log(Du(h,h));
          idx++;
        }
      }
      THETA.segment(NTHR+NLOAD+NCORR, NVAR) = norm_logsd;
    }

    idx=0;
    Eigen::VectorXd norm_lambda(NLOAD);
    for(unsigned int h = 0; h < q; h++){
      for(unsigned int j = 0; j < p; j++){
        if(!std::isfinite(CONSTRMAT(j,h)))
        {
          norm_lambda(idx) = Lam(j, h);
          idx ++;
        };
      };
    }
    THETA.segment(NTHR, NLOAD) = norm_lambda;

  }
}

#include "genericUtils.h"

/* EXTRACT PI THRESHOLDS */
// Extract from tau the thresholds related to pi_sksl
Eigen::VectorXd params::extract_thresholds(const Eigen::Ref<const Eigen::VectorXd> TAU,
                                   const Eigen::Ref<const Eigen::VectorXd> C_VEC,
                                   const unsigned int K,
                                   const unsigned int L,
                                   const unsigned int SK,
                                   const unsigned int SL){

  unsigned int ck = C_VEC(K);
  unsigned int cl = C_VEC(L);

  // identify tau_sk, tau_sl, tau_sk-1, tau_sl-1
  unsigned int sk_tau_index = C_VEC.segment(0, K).sum() - (K) + SK;   // index tau_sk in tau vector
  unsigned int sl_tau_index = C_VEC.segment(0, L).sum() - (L) + SL;   // index tau_sl in tau vector
  double t_sk;                                               // tau_sk
  if(SK == ck-1){
    t_sk = 100;
  } else {
    t_sk = TAU(sk_tau_index);
  }
  double t_sk_prev;                                           // tau_sk-1
  if(SK == 0){
    t_sk_prev = -100;
  }else{
    t_sk_prev = TAU(sk_tau_index-1);
  }

  double t_sl;                                               // tau_sl
  if(SL == cl-1){
    t_sl = 100;
  } else {
    t_sl = TAU(sl_tau_index);
  }
  double t_sl_prev;                                           // tau_sl-1
  if(SL == 0){
    t_sl_prev = -100;
  }else{
    t_sl_prev = TAU(sl_tau_index-1);
  }

  Eigen::VectorXd pi_thresholds(4);
  pi_thresholds << t_sk, t_sl, t_sk_prev, t_sl_prev;
  return pi_thresholds;
}


  namespace params::thresholds{
    Eigen::VectorXd theta2vec(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const int NTHR){
      return THETA.segment(0, NTHR);
      }
  }
  namespace params::loadings{
    Eigen::VectorXd theta2vec(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const int NTHR,
        const int NLOAD
    ){
      return THETA.segment(NTHR, NLOAD);
    }

    Eigen::VectorXd mat2vec(const Eigen::Ref<const Eigen::MatrixXd> LOADINGS,
                            const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
                            const int  NLOAD){

      Eigen::VectorXd out(NLOAD);
      const int p = LOADINGS.rows();
      const int q = LOADINGS.cols();
      int idx = 0;

      // construct vector with parameters free to be estimated
      for(int h = 0; h < q; h++){
        for(int j = 0; j < p; j++){
          if(!std::isfinite(CONSTRMAT(j,h))){
            out(idx)=LOADINGS(j,h);
            idx ++;
          };
        };
      }

      return out;


    }

    Eigen::MatrixXd theta2mat(const Eigen::Ref<const Eigen::VectorXd> THETA,
                              const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
                              const int NTHR,
                              const int NLOAD){

      Eigen::VectorXd vec = params::loadings::theta2vec(THETA, NTHR, NLOAD);
      const int p = CONSTRMAT.rows();
      const int q = CONSTRMAT.cols();
      Eigen::MatrixXd out(p, q);



      int idx = 0;

      // construct vector with parameters free to be estimated
      for(int h = 0; h < q; h++){
        for(int j = 0; j < p; j++){
          if(!std::isfinite(CONSTRMAT(j,h))){
            out(j, h)=vec(idx);
            idx ++;
          }else{
            out(j, h)= double(CONSTRMAT(j, h));
          };
        };
      }



      return out;


    }


  }
  namespace params::latvar{
    Eigen::VectorXd theta2vec(const Eigen::Ref<const Eigen::VectorXd> THETA,
                              const int NTHR,
                              const int NLOAD,
                              const int NCORR,
                              const int NVAR){
      return THETA.segment(NTHR+NLOAD, NCORR+NVAR);
    }


    Eigen::MatrixXd vec2cmat(const Eigen::Ref<const Eigen::VectorXd> VEC,
                             const int NCORR,
                             const int Q){

      if(NCORR>0){
        // extract uncostrained parameters related to latent correlations

        // inverse of Fisher's transformation
        Eigen::VectorXd tanh_entries = ((Eigen::VectorXd(2*VEC.segment(0, NCORR))).array().exp() - 1)/
          ((Eigen::VectorXd(2*VEC.segment(0, NCORR))).array().exp() + 1);

        // place tanh_entries in upper triangular matrix
        // with unit diagonal
        Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(Q,Q);
        unsigned int iterator = 0;
        for(unsigned int col = 1; col < Q; col ++){
          for(unsigned int row = 0; row < col; row++){
            Z(row, col) = tanh_entries(iterator);
            iterator ++;
          }
        }

        // construct upper Cholesky from Z
        Eigen::MatrixXd U(Q,Q);

        for( int col = 0; col < Q; col++){
          for( int row = 0; row < Q; row++){
            if(row > col) U(row, col) = 0;
            else if ((row == col) && (row == 0)) U(row, col) = 1;
            else if ((row > 0) && (row == col)) {
              if(Z(row-1,col)==0){
                double prod = 1;
                for(unsigned int row1 = 0; row1 < row; row1++){
                  prod *= pow(1-pow(Z(row1, col), 2),.5);
                }
                U(row, col) = prod;
              }else{
                U(row, col) = U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5)/Z(row-1, col);
              }
            }
            else if ((row == 0) && (row < col)) U(row, col) = Z(row, col);
            else if ((row > 0) && (row < col)) {
              if(Z(row-1,col)==0){
                double prod = Z(row, col);
                for( int row1 = 0; row1 < row; row1++){
                  prod *= pow(1-pow(Z(row1, col), 2),.5);
                }
                U(row, col) = prod;
              }else{
                U(row, col) = Z(row, col)/Z(row-1, col)*U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5);
              }
            }
          }
        }


        // latent variable correlation matrix
        Eigen::MatrixXd R = U.transpose()*U;

        return R;
      }else{
        return Eigen::MatrixXd::Identity(Q,Q);
      }
     }

    Eigen::MatrixXd vec2dmat(const Eigen::Ref<const Eigen::VectorXd> VEC,
                             const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                             const int NCORR,
                             const int NVAR,
                             const int Q){
      Eigen::VectorXd out(Q);

      if(NVAR>0){
        int idx = NCORR;
        for(int j=0; j<Q; j++){
          if(!std::isfinite(CONSTRLOGSD(j))){
            out(j) = exp(VEC(idx));
            idx++;
          }else{
            out(j) = exp(CONSTRLOGSD(j));
          }
        }

      }else{
        out = CONSTRLOGSD.array().exp();
      }

      return out.asDiagonal();
    }


    Eigen::MatrixXd theta2cmat(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const int NTHR,
        const int NLOAD,
        const int NCORR,
        const int NVAR,
        const int Q){

      Eigen::VectorXd vec = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR);
      return params::latvar::vec2cmat(vec, NCORR, Q);
    }

    Eigen::MatrixXd theta2dmat(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
        const int NTHR,
        const int NLOAD,
        const int NCORR,
        const int NVAR,
        const int Q){

        if(NVAR>0){
          Eigen::VectorXd vec = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR);
          return params::latvar::vec2dmat(vec, CONSTRLOGSD, NCORR, NVAR, Q);
        }else{

          return Eigen::VectorXd(CONSTRLOGSD.array().exp()).asDiagonal();
        }
      }

    Eigen::MatrixXd theta2mat(const Eigen::Ref<const Eigen::VectorXd> THETA,
                              const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                              const int NTHR,
                              const int NLOAD,
                              const int NCORR,
                              const int NVAR,
                              const int Q){

      Eigen::VectorXd vec = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR);
      Eigen::MatrixXd R = params::latvar::vec2cmat(vec, NCORR, Q);
      Eigen::MatrixXd D = params::latvar::vec2dmat(vec, CONSTRLOGSD, NCORR, NVAR, Q);
      return D*R*D;
    }

    Eigen::VectorXd mat2edvec(const Eigen::Ref<const Eigen::MatrixXd> S){
      return S.diagonal().array().sqrt().log();
    }
    Eigen::VectorXd edvec2dvec(const Eigen::Ref<const Eigen::VectorXd> EDVEC,
                               const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                               const int Q,
                               const int NVAR){
      Eigen::VectorXd out(NVAR);
      int idx=0;
      for(int j=0; j<Q; j++){
        if(!std::isfinite(CONSTRLOGSD(j))){
          out(idx)=EDVEC(j);
          idx++;
        }
      }

      return(out);


    }

    Eigen::MatrixXd mat2cmat(const Eigen::Ref<const Eigen::MatrixXd> S){
      Eigen::VectorXd invD = S.diagonal().array().sqrt().inverse();
      return invD.asDiagonal() * S * invD.asDiagonal();
    }

    Eigen::VectorXd cmat2cvec(const Eigen::Ref<const Eigen::MatrixXd> R){
      const unsigned int q = R.cols(); // number of latent variables

      Eigen::LLT<Eigen::MatrixXd> llt(R);
      Eigen::MatrixXd U = llt.matrixU();

      Eigen::VectorXd y(q*(q-1)/2);
      Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(q,q);
      unsigned int iterator = 0;
      for(unsigned int col = 1; col < q; col ++){
        for(unsigned int row = 0; row < col; row++){
          if(row == 0) Z(row, col) = U(row, col);
          else if(row > 0){
            double prod = 1;
            for(unsigned int i = 0; i < row; i++){
              prod *= pow(1-pow(Z(i, col), 2),-.5);
              Z(row, col) = U(row, col)*prod;
            }
          }
          y(iterator) = atanh(Z(row, col));
          iterator++;
        }
      }

      return y;
    }

    Eigen::VectorXd mat2vec(const Eigen::Ref<const Eigen::MatrixXd> S,
                            const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                            const int NCORR,
                            const int NVAR){

      Eigen::VectorXd out(NCORR+NVAR);
      int q = S.cols();

      if(NCORR>0){
        Eigen::VectorXd cvec = params::latvar::cmat2cvec(params::latvar::mat2cmat(S));
        if(NVAR>0){
          Eigen::VectorXd dvec = params::latvar::edvec2dvec(params::latvar::mat2edvec(S),
                                                            CONSTRLOGSD, q, NVAR);
          out << cvec, dvec;
        }else{
          out = cvec;
        }
      }else{
        if(NVAR>0){
          Eigen::VectorXd dvec = params::latvar::edvec2dvec(params::latvar::mat2edvec(S),
                                                            CONSTRLOGSD, q, NVAR);
          out = dvec;
        }else{
          out.resize(0);
        }

      }

      return out;

    }






  }

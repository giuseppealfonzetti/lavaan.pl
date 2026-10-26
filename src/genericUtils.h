#ifndef utils_H
#define utils_H
#include <RcppEigen.h>

namespace params{


  /* EXTRACT PI THRESHOLDS */
  // Extract from tau the thresholds related to pi_sksl
  Eigen::VectorXd extract_thresholds(const Eigen::Ref<const Eigen::VectorXd> TAU,
                                     const Eigen::Ref<const Eigen::VectorXd> C_VEC,
                                     const unsigned int K,
                                     const unsigned int L,
                                     const unsigned int SK,
                                     const unsigned int SL);



  namespace thresholds{
    Eigen::VectorXd theta2vec(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const int NTHR);
  }
  namespace loadings{
    Eigen::VectorXd theta2vec(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const int NTHR,
        const int NLOAD
    );

    Eigen::VectorXd mat2vec(const Eigen::Ref<const Eigen::MatrixXd> LOADINGS,
                            const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
                            const int  NLOAD);

    Eigen::MatrixXd theta2mat(const Eigen::Ref<const Eigen::VectorXd> THETA,
                              const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
                              const int NTHR,
                              const int NLOAD);


  }
  namespace latvar{
    Eigen::VectorXd theta2vec(const Eigen::Ref<const Eigen::VectorXd> THETA,
                              const int NTHR,
                              const int NLOAD,
                              const int NCORR,
                              const int NVAR);


    Eigen::MatrixXd vec2cmat(const Eigen::Ref<const Eigen::VectorXd> VEC,
                             const int NCORR,
                             const int Q);

    Eigen::MatrixXd vec2dmat(const Eigen::Ref<const Eigen::VectorXd> VEC,
                             const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                             const int NCORR,
                             const int NVAR,
                             const int Q);


    Eigen::MatrixXd theta2cmat(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const int NTHR,
        const int NLOAD,
        const int NCORR,
        const int NVAR,
        const int Q);

    Eigen::MatrixXd theta2dmat(
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
        const int NTHR,
        const int NLOAD,
        const int NCORR,
        const int NVAR,
        const int Q);

    Eigen::MatrixXd theta2mat(const Eigen::Ref<const Eigen::VectorXd> THETA,
                              const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                              const int NTHR,
                              const int NLOAD,
                              const int NCORR,
                              const int NVAR,
                              const int Q);

    Eigen::VectorXd mat2edvec(const Eigen::Ref<const Eigen::MatrixXd> S);
    Eigen::VectorXd edvec2dvec(const Eigen::Ref<const Eigen::VectorXd> EDVEC,
                               const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                               const int Q,
                               const int NVAR);

    Eigen::MatrixXd mat2cmat(const Eigen::Ref<const Eigen::MatrixXd> S);

    Eigen::VectorXd cmat2cvec(const Eigen::Ref<const Eigen::MatrixXd> R);

    Eigen::VectorXd mat2vec(const Eigen::Ref<const Eigen::MatrixXd> S,
                            const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
                            const int NCORR,
                            const int NVAR);






  }


}


#endif

#ifndef fullPairwise_H
#define fullPairwise_H
#include "parallelWorkers.h"

namespace fullPairwise{

  Rcpp::List multiThreadContribution(
      const int N,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const int GRFLAG,
      const int SILENTFLAG,
      const Eigen::Ref<const Eigen::VectorXd> UFREQ,
      const Eigen::Ref<const Eigen::MatrixXd> REGMAT,
      const Eigen::Ref<const Eigen::MatrixXd> RESCOV
  );




  Rcpp::List multiThreadSamplEstHJ(
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const Eigen::Ref<const Eigen::MatrixXd> DATA,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      int N,
      int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const Eigen::Ref<const Eigen::VectorXd> UFREQ,
      const Eigen::Ref<const Eigen::VectorXd> UW,
      const Eigen::Ref<const Eigen::VectorXd> GW,
      const Eigen::Ref<const Eigen::VectorXd> JW,
      const Eigen::Ref<const Eigen::VectorXd> CL,
      const Eigen::Ref<const Eigen::MatrixXd> REGMAT,
      const Eigen::Ref<const Eigen::MatrixXd> RESCOV
  );

  Rcpp::List multiThreadSaturated(
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const Eigen::Ref<const Eigen::VectorXd> RHO,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const Eigen::Ref<const Eigen::MatrixXd> DATA,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const int ORDER,
      const bool SCORES,
      const Eigen::Ref<const Eigen::VectorXd> UFREQ,
      const Eigen::Ref<const Eigen::VectorXd> UW
  );
}

#endif

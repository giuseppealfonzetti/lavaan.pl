#ifndef fullPairwise_H
#define fullPairwise_H
#include "parallelWorkers.h"

namespace fullPairwise{

  Rcpp::List multiThreadContribution(
      const int N,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const std::vector<std::vector<std::vector<double>>> LLC,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const int GRFLAG,
      const int SILENTFLAG
  );




  Rcpp::List multiThreadSamplEstHJ(
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const Eigen::Ref<const Eigen::MatrixXd> DATA,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const std::vector<std::vector<std::vector<double>>> LLC,
      int N,
      int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR
  );

  Rcpp::List multiThreadSaturated(
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const Eigen::Ref<const Eigen::VectorXd> RHO,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const Eigen::Ref<const Eigen::MatrixXd> DATA,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const int ORDER,
      const bool SCORES
  );
}

#endif

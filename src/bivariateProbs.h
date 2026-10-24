#ifndef bivariateProbs_H
#define bivariateProbs_H

#include <RcppEigen.h>
#include "bivariateNormal.h"

namespace biprobs{
  /*  PI (the probability of a bivariate response pattern) */
  // Compute specific pi_sksl
  double compute_pi(
      const Eigen::Ref<const Eigen::VectorXd> c_vec,
      const Eigen::Ref<const Eigen::VectorXd> pi_thresholds,
      double rho_kl,

      const unsigned int k,
      const unsigned int l,
      const unsigned int sk,
      const unsigned int sl
  );

}




#endif

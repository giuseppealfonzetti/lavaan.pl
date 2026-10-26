#ifndef optimisationUtils_H
#define optimisationUtils_H
#include <random>

#include "genericUtils.h"

namespace sa{


  void proj2(
        Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
        Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
        Eigen::Map<Eigen::VectorXd> C_VEC,
        const int CORRFLAG,
        const int NTHR,
        const int NLOAD,
        const int NCORR,
        const int NVAR,
        Eigen::VectorXd &THETA,
        bool &CHECKEVENT);




}

































#endif



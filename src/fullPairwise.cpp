#include "fullPairwise.h"

Rcpp::List fullPairwise::multiThreadContribution(
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
    const Eigen::Ref<const Eigen::VectorXd> UFREQ
){

  Rcpp::checkUserInterrupt();
  const int d = THETA.size();
  const int p = CONSTRMAT.rows();                                             // number of items
  int R = p*(p-1)/2;                                            // number of pairs of items

  // Copy frequencies, and build pair dictionary
  Eigen::MatrixXd items_pairs(2,R);
  int r = 0;
  for(int k = 1; k < p; k ++){
    for(int l = 0; l < k; l++){
      items_pairs(0, r) = k;
      items_pairs(1, r) = l;
      r++;
    }
  }

  // Rearrange parameters
  Eigen::VectorXd theta(d);
  theta << THETA;                                 // Complete parameters vector

  // Initialize vector of indices for entries in pairs_table
  std::vector<int> vector_pairs(items_pairs.cols()) ;
  std::iota (std::begin(vector_pairs), std::end(vector_pairs), 0);

  double iter_ll = 0;
  Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);

  parallelWorkers::pairsGradient fullpool(CONSTRMAT, CONSTRLOGSD, C_VEC, FREQ, items_pairs, CORRFLAG, NTHR, NLOAD, NCORR, NVAR,
                                          SILENTFLAG, GRFLAG, theta, vector_pairs);
  RcppParallel::parallelReduce(0, R, fullpool);
  iter_ll = fullpool.subset_ll;
  iter_gradient = fullpool.subset_gradient;
  Eigen::MatrixXd hdummy, udummy;
  pairs::univariate_contribution(theta.head(NTHR), C_VEC, UFREQ, GRFLAG, iter_ll, iter_gradient, hdummy, udummy);

  // output list
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("iter_nll") = -iter_ll,
      Rcpp::Named("iter_ngradient") = -iter_gradient
    );
  return(output);
}




Rcpp::List fullPairwise::multiThreadSamplEstHJ(
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
    const Eigen::Ref<const Eigen::VectorXd> CL
){
  const int d  = THETA.size();
  const int pp = FREQ.cols();
  const int n  = DATA.rows();
  const int p  = CONSTRMAT.rows();

  Eigen::MatrixXd gradmat = Eigen::MatrixXd::Zero(d, pp);
  parallelWorkers::patternwiseOuterProd estH(
      CONSTRMAT,
      CONSTRLOGSD,
      C_VEC,
      FREQ,
      CORRFLAG,
      NTHR,
      NLOAD,
      NCORR,
      NVAR,
      THETA,
      gradmat
  );

  RcppParallel::parallelReduce(0, p * (p - 1) / 2, estH);
  Eigen::VectorXd gradH = estH.gradient;
  Eigen::MatrixXd H = estH.hessian;
  double ll = 0;
  Eigen::MatrixXd uscore = Eigen::MatrixXd::Zero(NTHR, UFREQ.size());
  pairs::univariate_contribution(THETA.head(NTHR), C_VEC, UFREQ, 2, ll, gradH, H, uscore);


  parallelWorkers::casewiseOuterProd estJ(
      DATA,
      gradmat,
      C_VEC,
      uscore,
      UW,
      GW,
      JW,
      CL,
      CL.size() > 0 ? static_cast<int>(CL.maxCoeff()) + 1 : 0
  );

  RcppParallel::parallelReduce(0, n, estJ);
  Eigen::VectorXd gradJ = estJ.gradient;
  Eigen::MatrixXd J = estJ.variability;
  if(CL.size() > 0){
    const double G = estJ.ztot.cols();
    const Eigen::MatrixXd zc = estJ.ztot.colwise() - estJ.ztot.rowwise().mean();
    J = G / (G - 1) * zc * zc.transpose();
  }


  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gradmat") = gradmat,
      Rcpp::Named("gradH") = gradH,
      Rcpp::Named("gradJ") = gradJ,
      Rcpp::Named("H") = H/n,
      Rcpp::Named("J") = J/n
    );

  return(output);
}

Rcpp::List fullPairwise::multiThreadSaturated(
    const Eigen::Ref<const Eigen::VectorXd> TAU,
    const Eigen::Ref<const Eigen::VectorXd> RHO,
    const Eigen::Ref<const Eigen::MatrixXd> FREQ,
    const Eigen::Ref<const Eigen::MatrixXd> DATA,
    const Eigen::Ref<const Eigen::VectorXd> C_VEC,
    const int ORDER,
    const bool SCORES,
    const Eigen::Ref<const Eigen::VectorXd> UFREQ,
    const Eigen::Ref<const Eigen::VectorXd> UW
){
  const int nthr = TAU.size();
  const int P    = RHO.size();
  const int n    = DATA.rows();
  const int order = SCORES ? std::max(ORDER, 1) : ORDER;

  Eigen::MatrixXd patscore = Eigen::MatrixXd::Zero(SCORES ? 2 * C_VEC.maxCoeff() - 1 : 0, SCORES ? FREQ.cols() : 0);
  parallelWorkers::saturatedPairs sat(TAU, RHO, C_VEC, FREQ, order, SCORES, patscore);
  RcppParallel::parallelReduce(0, P, sat);
  Eigen::MatrixXd uscore = Eigen::MatrixXd::Zero(SCORES ? nthr : 0, SCORES ? UFREQ.size() : 0);
  pairs::univariate_contribution(TAU, C_VEC, UFREQ, order, sat.ll, sat.gradient, sat.htt, uscore);

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("ll") = sat.ll,
      Rcpp::Named("gradient") = sat.gradient,
      Rcpp::Named("htt") = sat.htt,
      Rcpp::Named("htr") = sat.htr,
      Rcpp::Named("hrr") = sat.hrr
    );

  if(SCORES){
    Eigen::MatrixXd st = Eigen::MatrixXd::Zero(nthr, n);
    Eigen::MatrixXd sr = Eigen::MatrixXd::Zero(P, n);
    parallelWorkers::saturatedScores cases(DATA, patscore, C_VEC, uscore, UW, st, sr);
    RcppParallel::parallelFor(0, n, cases);
    output["st"] = st;
    output["sr"] = sr;
  }

  return(output);
}

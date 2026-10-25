#ifndef paralellWorkers_H
#define paralellWorkers_H
#include <RcppEigen.h>
#include <RcppParallel.h>
#include "pairs.h"

namespace parallelWorkers{

  // RcppParallel Worker to compute pair_contribution
  // on multiple pairs in parallel
  struct pairsGradient : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Ref<const Eigen::MatrixXd> constrmat;
    const Eigen::Ref<const Eigen::VectorXd> constrsd;
    const std::vector<std::vector<std::vector<double>>> llc;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;
    const Eigen::Ref<const Eigen::MatrixXd> pairs_table;
    const Eigen::Ref<const Eigen::MatrixXd> items_pairs;
    const int corrFLAG;
    const int nthr;
    const int nload;
    const int ncorr;
    const int nvar;
    const int silentFLAG;
    const int gradFLAG;

    //// Iteration:
    const Eigen::VectorXd &theta;
    const std::vector<int> &index_vector;

    // Output quantities:
    double subset_ll;
    Eigen::VectorXd subset_gradient;

    // Constructor 1:
    pairsGradient(
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT_,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRSD_,
      const std::vector<std::vector<std::vector<double>>> LLC_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_,
      const Eigen::Ref<const Eigen::MatrixXd> PAIRS_TABLE_,
      const Eigen::Ref<const Eigen::MatrixXd> ITEMS_PAIRS_,
      const int CORRFLAG_,
      const int NTHR_,
      const int NLOAD_,
      const int NCORR_,
      const int NVAR_,
      const int SILENTFLAG_,
      const int GRADFLAG_,
      const Eigen::VectorXd &THETA_,
      const std::vector<int> &INDEX_VECTOR_
    ):
      constrmat(CONSTRMAT_), constrsd(CONSTRSD_), llc(LLC_), c_vec(C_VEC_), pairs_table(PAIRS_TABLE_), items_pairs(ITEMS_PAIRS_),
      corrFLAG(CORRFLAG_), nthr(NTHR_), nload(NLOAD_), ncorr(NCORR_), nvar(NVAR_), silentFLAG(SILENTFLAG_), gradFLAG(GRADFLAG_),
      theta(THETA_), index_vector(INDEX_VECTOR_), subset_ll(0.0),
      subset_gradient(Eigen::VectorXd::Zero(THETA_.size())){}

    // Constructor 2:
    pairsGradient(const pairsGradient &OBJ_, RcppParallel::Split):
      constrmat(OBJ_.constrmat), constrsd(OBJ_.constrsd), llc(OBJ_.llc), c_vec(OBJ_.c_vec), pairs_table(OBJ_.pairs_table), items_pairs(OBJ_.items_pairs),
      corrFLAG(OBJ_.corrFLAG), nthr(OBJ_.nthr), nload(OBJ_.nload), ncorr(OBJ_.ncorr), nvar(OBJ_.nvar),
      silentFLAG(OBJ_.silentFLAG), gradFLAG(OBJ_.gradFLAG),
      theta(OBJ_.theta), index_vector(OBJ_.index_vector), subset_ll(0.0),
      subset_gradient(Eigen::VectorXd::Zero(theta.size())){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const pairsGradient &RHS);

  };

  // RcppParallel Worker to compute gradient outer product of
  // bivariate response patterns. Used to evaluate the
  // sample estimator of H
  struct patternwiseOuterProd : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Ref<const Eigen::MatrixXd> constrmat;
    const Eigen::Ref<const Eigen::VectorXd> constrsd;
    const std::vector<std::vector<std::vector<double>>> llc;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;
    const Eigen::Ref<const Eigen::MatrixXd> freq;
    const int corrFLAG;
    const int nthr;
    const int nload;
    const int ncorr;
    const int nvar;

    //// Iteration:
    const Eigen::Ref<const Eigen::VectorXd> theta;
    Eigen::MatrixXd &gradmat;

    // Output quantities:
    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;


    // Constructor 1:
    patternwiseOuterProd(
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT_,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRSD_,
      const std::vector<std::vector<std::vector<double>>> LLC_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ_,
      const int CORRFLAG_,
      const int NTHR_,
      const int NLOAD_,
      const int NCORR_,
      const int NVAR_,
      const Eigen::Ref<const Eigen::VectorXd> THETA_,
      Eigen::MatrixXd &GRADMAT_

    ):
      constrmat(CONSTRMAT_), constrsd(CONSTRSD_), llc(LLC_), c_vec(C_VEC_), freq(FREQ_),
      corrFLAG(CORRFLAG_), nthr(NTHR_), nload(NLOAD_), ncorr(NCORR_), nvar(NVAR_),
      theta(THETA_),
      gradmat(GRADMAT_),
      gradient(Eigen::VectorXd::Zero(THETA_.size())),
      hessian(Eigen::MatrixXd::Zero(THETA_.size(), THETA_.size()))
      {}

    // Constructor 2:
    patternwiseOuterProd(const patternwiseOuterProd &OBJ_, RcppParallel::Split):
      constrmat(OBJ_.constrmat), constrsd(OBJ_.constrsd), llc(OBJ_.llc), c_vec(OBJ_.c_vec), freq(OBJ_.freq),
      corrFLAG(OBJ_.corrFLAG), nthr(OBJ_.nthr), nload(OBJ_.nload), ncorr(OBJ_.ncorr), nvar(OBJ_.nvar),
      theta(OBJ_.theta),
      gradmat(OBJ_.gradmat),
      gradient(Eigen::VectorXd::Zero(theta.size())),
      hessian(Eigen::MatrixXd::Zero(theta.size(), theta.size())){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const patternwiseOuterProd &RHS);

  };

  // RcppParallel Worker to estimate H and J
  struct casewiseOuterProd : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Ref<const Eigen::MatrixXd> data;
    const Eigen::Ref<const Eigen::MatrixXd> gradmat;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;
    const Eigen::Ref<const Eigen::MatrixXd> uscore;
    const Eigen::Ref<const Eigen::VectorXd> uw;
    // Output quantities:
    Eigen::VectorXd gradient;
    Eigen::MatrixXd variability;



    // Constructor 1:
    casewiseOuterProd(
      const Eigen::Ref<const Eigen::MatrixXd> DATA_,
      const Eigen::Ref<const Eigen::MatrixXd> GRADMAT_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_,
      const Eigen::Ref<const Eigen::MatrixXd> USCORE_,
      const Eigen::Ref<const Eigen::VectorXd> UW_
    ):
      data(DATA_), gradmat(GRADMAT_), c_vec(C_VEC_), uscore(USCORE_), uw(UW_),
      gradient(Eigen::VectorXd::Zero(GRADMAT_.rows())),
      variability(Eigen::MatrixXd::Zero(GRADMAT_.rows(), GRADMAT_.rows())){}

    // Constructor 2:
    casewiseOuterProd(const casewiseOuterProd &OBJ_, RcppParallel::Split):
      data(OBJ_.data), gradmat(OBJ_.gradmat), c_vec(OBJ_.c_vec), uscore(OBJ_.uscore), uw(OBJ_.uw),
      gradient(Eigen::VectorXd::Zero(gradmat.rows())),
      variability(Eigen::MatrixXd::Zero(gradmat.rows(), gradmat.rows())){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const casewiseOuterProd &RHS);

  };

  // RcppParallel Worker for the saturated model, parametrised by
  // thresholds and one polychoric correlation per pair.
  // Returns log-likelihood, gradient and the blocks of the observed
  // information: thresholds-thresholds, thresholds-correlations and
  // the (diagonal) correlations-correlations block
  struct saturatedPairs : public RcppParallel::Worker{
    const Eigen::Ref<const Eigen::VectorXd> tau;
    const Eigen::Ref<const Eigen::VectorXd> rho;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;
    const Eigen::Ref<const Eigen::MatrixXd> freq;
    const int order;
    const bool scores;
    Eigen::MatrixXd &patscore;

    double ll;
    Eigen::VectorXd gradient;
    Eigen::MatrixXd htt;
    Eigen::MatrixXd htr;
    Eigen::VectorXd hrr;

    saturatedPairs(
      const Eigen::Ref<const Eigen::VectorXd> TAU_,
      const Eigen::Ref<const Eigen::VectorXd> RHO_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ_,
      const int ORDER_,
      const bool SCORES_,
      Eigen::MatrixXd &PATSCORE_
    ):
      tau(TAU_), rho(RHO_), c_vec(C_VEC_), freq(FREQ_), order(ORDER_), scores(SCORES_), patscore(PATSCORE_),
      ll(0.0),
      gradient(Eigen::VectorXd::Zero(TAU_.size() + RHO_.size())),
      htt(Eigen::MatrixXd::Zero(TAU_.size(), TAU_.size())),
      htr(Eigen::MatrixXd::Zero(TAU_.size(), RHO_.size())),
      hrr(Eigen::VectorXd::Zero(RHO_.size())){}

    saturatedPairs(const saturatedPairs &OBJ_, RcppParallel::Split):
      tau(OBJ_.tau), rho(OBJ_.rho), c_vec(OBJ_.c_vec), freq(OBJ_.freq), order(OBJ_.order), scores(OBJ_.scores), patscore(OBJ_.patscore),
      ll(0.0),
      gradient(Eigen::VectorXd::Zero(tau.size() + rho.size())),
      htt(Eigen::MatrixXd::Zero(tau.size(), tau.size())),
      htr(Eigen::MatrixXd::Zero(tau.size(), rho.size())),
      hrr(Eigen::VectorXd::Zero(rho.size())){}

    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const saturatedPairs &RHS);
  };

  // RcppParallel Worker collecting the casewise scores of the
  // saturated model from the patternwise scores
  struct saturatedScores : public RcppParallel::Worker{
    const Eigen::Ref<const Eigen::MatrixXd> data;
    const Eigen::Ref<const Eigen::MatrixXd> patscore;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;
    const Eigen::Ref<const Eigen::MatrixXd> uscore;
    const Eigen::Ref<const Eigen::VectorXd> uw;
    Eigen::MatrixXd &st;
    Eigen::MatrixXd &sr;

    saturatedScores(
      const Eigen::Ref<const Eigen::MatrixXd> DATA_,
      const Eigen::Ref<const Eigen::MatrixXd> PATSCORE_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_,
      const Eigen::Ref<const Eigen::MatrixXd> USCORE_,
      const Eigen::Ref<const Eigen::VectorXd> UW_,
      Eigen::MatrixXd &ST_,
      Eigen::MatrixXd &SR_
    ):
      data(DATA_), patscore(PATSCORE_), c_vec(C_VEC_), uscore(USCORE_), uw(UW_), st(ST_), sr(SR_){}

    void operator()(std::size_t BEGIN, std::size_t END);
  };

}

#endif

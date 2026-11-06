#ifndef workers_H
#define workers_H
#include <RcppEigen.h>
#include <RcppParallel.h>
#include "pairs.h"

namespace workers{

  struct PairLoglik : public RcppParallel::Worker{
    const model::Matrices &par;
    const pairs::GroupData &data;
    const std::vector<int> &pair_indices;
    const pairs::Order order;
    double ll;
    Eigen::VectorXd gradient;

    PairLoglik(const model::Matrices &PAR, const pairs::GroupData &DATA, const std::vector<int> &PAIR_INDICES, const pairs::Order ORDER):
      par(PAR), data(DATA), pair_indices(PAIR_INDICES), order(ORDER), ll(0), gradient(Eigen::VectorXd::Zero(PAR.structure.d)){}
    PairLoglik(const PairLoglik &O, RcppParallel::Split):
      par(O.par), data(O.data), pair_indices(O.pair_indices), order(O.order), ll(0), gradient(Eigen::VectorXd::Zero(O.par.structure.d)){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const PairLoglik &O);
  };

  struct PairHessian : public RcppParallel::Worker{
    const model::Matrices &par;
    const pairs::GroupData &data;
    std::vector<Eigen::MatrixXd> &cell_scores;
    std::vector<std::vector<int>> &score_index;
    std::vector<Eigen::VectorXd> &pair_gradient;
    Eigen::MatrixXd hessian;
    Eigen::VectorXd gradient;

    PairHessian(const model::Matrices &PAR, const pairs::GroupData &DATA, std::vector<Eigen::MatrixXd> &CELL_SCORES,
                std::vector<std::vector<int>> &SCORE_INDEX, std::vector<Eigen::VectorXd> &PAIR_GRADIENT):
      par(PAR), data(DATA), cell_scores(CELL_SCORES), score_index(SCORE_INDEX), pair_gradient(PAIR_GRADIENT),
      hessian(Eigen::MatrixXd::Zero(PAR.structure.d, PAR.structure.d)), gradient(Eigen::VectorXd::Zero(PAR.structure.d)){}
    PairHessian(const PairHessian &O, RcppParallel::Split):
      par(O.par), data(O.data), cell_scores(O.cell_scores), score_index(O.score_index), pair_gradient(O.pair_gradient),
      hessian(Eigen::MatrixXd::Zero(O.par.structure.d, O.par.structure.d)), gradient(Eigen::VectorXd::Zero(O.par.structure.d)){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const PairHessian &O);
  };

  struct CaseScores : public RcppParallel::Worker{
    const model::Matrices &par;
    const pairs::GroupData &data;
    const std::vector<Eigen::MatrixXd> &cell_scores;
    const std::vector<std::vector<int>> &score_index;
    const Eigen::Ref<const Eigen::MatrixXd> univariate_scores;
    const Eigen::Ref<const Eigen::VectorXd> uweights;
    Eigen::MatrixXd &scores;

    CaseScores(const model::Matrices &PAR, const pairs::GroupData &DATA, const std::vector<Eigen::MatrixXd> &CELL_SCORES,
               const std::vector<std::vector<int>> &SCORE_INDEX, const Eigen::Ref<const Eigen::MatrixXd> UNIVARIATE_SCORES,
               const Eigen::Ref<const Eigen::VectorXd> UWEIGHTS, Eigen::MatrixXd &SCORES):
      par(PAR), data(DATA), cell_scores(CELL_SCORES), score_index(SCORE_INDEX), univariate_scores(UNIVARIATE_SCORES),
      uweights(UWEIGHTS), scores(SCORES){}
    void operator()(std::size_t BEGIN, std::size_t END);
  };

  struct SaturatedLoglik : public RcppParallel::Worker{
    const model::Items &items;
    const pairs::GroupData &data;
    const pairs::Reduced &reduced;
    const Eigen::Ref<const Eigen::VectorXd> rho;
    const pairs::Order order;
    const bool scores;
    Eigen::MatrixXd &pattern_scores;
    double ll;
    Eigen::VectorXd gradient, hrr;
    Eigen::MatrixXd htt, htr, hschur;

    SaturatedLoglik(const model::Items &ITEMS, const pairs::GroupData &DATA, const pairs::Reduced &REDUCED,
                    const Eigen::Ref<const Eigen::VectorXd> RHO, const int NITEM, const pairs::Order ORDER, const bool SCORES,
                    Eigen::MatrixXd &PATTERN_SCORES):
      items(ITEMS), data(DATA), reduced(REDUCED), rho(RHO), order(ORDER), scores(SCORES), pattern_scores(PATTERN_SCORES), ll(0),
      gradient(Eigen::VectorXd::Zero(NITEM + RHO.size())), hrr(Eigen::VectorXd::Zero(RHO.size())),
      htt(Eigen::MatrixXd::Zero(NITEM, NITEM)), htr(Eigen::MatrixXd::Zero(NITEM, RHO.size())), hschur(Eigen::MatrixXd::Zero(NITEM, NITEM)){}
    SaturatedLoglik(const SaturatedLoglik &O, RcppParallel::Split):
      items(O.items), data(O.data), reduced(O.reduced), rho(O.rho), order(O.order), scores(O.scores),
      pattern_scores(O.pattern_scores), ll(0), gradient(Eigen::VectorXd::Zero(O.gradient.size())),
      hrr(Eigen::VectorXd::Zero(O.hrr.size())), htt(Eigen::MatrixXd::Zero(O.htt.rows(), O.htt.cols())),
      htr(Eigen::MatrixXd::Zero(O.htr.rows(), O.htr.cols())), hschur(Eigen::MatrixXd::Zero(O.htt.rows(), O.htt.cols())){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const SaturatedLoglik &O);
  };

  struct SaturatedCaseScores : public RcppParallel::Worker{
    const model::Items &items;
    const pairs::GroupData &data;
    const pairs::Reduced &reduced;
    const Eigen::Ref<const Eigen::VectorXd> rho;
    const Eigen::Ref<const Eigen::MatrixXd> pattern_scores, univariate_scores;
    const Eigen::Ref<const Eigen::VectorXd> uweights;
    Eigen::MatrixXd &score_tau, &score_rho;

    SaturatedCaseScores(const model::Items &ITEMS, const pairs::GroupData &DATA, const pairs::Reduced &REDUCED,
                        const Eigen::Ref<const Eigen::VectorXd> RHO, const Eigen::Ref<const Eigen::MatrixXd> PATTERN_SCORES,
                        const Eigen::Ref<const Eigen::MatrixXd> UNIVARIATE_SCORES, const Eigen::Ref<const Eigen::VectorXd> UWEIGHTS,
                        Eigen::MatrixXd &SCORE_TAU, Eigen::MatrixXd &SCORE_RHO):
      items(ITEMS), data(DATA), reduced(REDUCED), rho(RHO), pattern_scores(PATTERN_SCORES), univariate_scores(UNIVARIATE_SCORES),
      uweights(UWEIGHTS), score_tau(SCORE_TAU), score_rho(SCORE_RHO){}
    void operator()(std::size_t BEGIN, std::size_t END);
  };

}

#endif

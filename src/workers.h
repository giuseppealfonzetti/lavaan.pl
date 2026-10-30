#ifndef workers_H
#define workers_H
#include <RcppEigen.h>
#include <RcppParallel.h>
#include "pairs.h"

namespace workers{

  struct PL : public RcppParallel::Worker{
    const model::Setup &st;
    const Eigen::Ref<const Eigen::VectorXd> freq;
    const std::vector<int> &rs;
    const int order;
    double ll;
    Eigen::VectorXd gradient;

    PL(const model::Setup &ST, const Eigen::Ref<const Eigen::VectorXd> FREQ, const std::vector<int> &RS, const int ORDER):
      st(ST), freq(FREQ), rs(RS), order(ORDER), ll(0), gradient(Eigen::VectorXd::Zero(ST.S.d)){}
    PL(const PL &O, RcppParallel::Split):
      st(O.st), freq(O.freq), rs(O.rs), order(O.order), ll(0), gradient(Eigen::VectorXd::Zero(O.st.S.d)){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const PL &O);
  };

  struct Hess : public RcppParallel::Worker{
    const model::Setup &st;
    const Eigen::Ref<const Eigen::VectorXd> freq;
    std::vector<Eigen::MatrixXd> &score;
    std::vector<std::vector<int>> &sidx;
    Eigen::MatrixXd hessian;
    Eigen::VectorXd gradient;

    Hess(const model::Setup &ST, const Eigen::Ref<const Eigen::VectorXd> FREQ, std::vector<Eigen::MatrixXd> &SCORE,
         std::vector<std::vector<int>> &SIDX):
      st(ST), freq(FREQ), score(SCORE), sidx(SIDX),
      hessian(Eigen::MatrixXd::Zero(ST.S.d, ST.S.d)), gradient(Eigen::VectorXd::Zero(ST.S.d)){}
    Hess(const Hess &O, RcppParallel::Split):
      st(O.st), freq(O.freq), score(O.score), sidx(O.sidx),
      hessian(Eigen::MatrixXd::Zero(O.st.S.d, O.st.S.d)), gradient(Eigen::VectorXd::Zero(O.st.S.d)){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const Hess &O);
  };

  struct Cases : public RcppParallel::Worker{
    const model::Layout &L;
    const Eigen::Ref<const Eigen::MatrixXd> y;
    const std::vector<Eigen::MatrixXd> &score;
    const std::vector<std::vector<int>> &sidx;
    const Eigen::Ref<const Eigen::MatrixXd> uscore;
    const Eigen::Ref<const Eigen::VectorXd> uw;
    Eigen::MatrixXd &out;

    Cases(const model::Layout &LAY, const Eigen::Ref<const Eigen::MatrixXd> Y, const std::vector<Eigen::MatrixXd> &SCORE,
          const std::vector<std::vector<int>> &SIDX, const Eigen::Ref<const Eigen::MatrixXd> USCORE,
          const Eigen::Ref<const Eigen::VectorXd> UW, Eigen::MatrixXd &CASES):
      L(LAY), y(Y), score(SCORE), sidx(SIDX), uscore(USCORE), uw(UW), out(CASES){}
    void operator()(std::size_t BEGIN, std::size_t END);
  };

  struct Saturated : public RcppParallel::Worker{
    const model::Layout &L;
    const Eigen::Ref<const Eigen::VectorXd> tau, rho, freq;
    const int order;
    const bool scores;
    Eigen::MatrixXd &patscore;
    double ll;
    Eigen::VectorXd gradient, hrr;
    Eigen::MatrixXd htt, htr, hschur;

    Saturated(const model::Layout &LAY, const Eigen::Ref<const Eigen::VectorXd> TAU, const Eigen::Ref<const Eigen::VectorXd> RHO,
              const Eigen::Ref<const Eigen::VectorXd> FREQ, const int ORDER, const bool SCORES, Eigen::MatrixXd &PATSCORE):
      L(LAY), tau(TAU), rho(RHO), freq(FREQ), order(ORDER), scores(SCORES), patscore(PATSCORE), ll(0),
      gradient(Eigen::VectorXd::Zero(TAU.size() + RHO.size())), hrr(Eigen::VectorXd::Zero(RHO.size())),
      htt(Eigen::MatrixXd::Zero(TAU.size(), TAU.size())), htr(Eigen::MatrixXd::Zero(TAU.size(), RHO.size())),
      hschur(Eigen::MatrixXd::Zero(TAU.size(), TAU.size())){}
    Saturated(const Saturated &O, RcppParallel::Split):
      L(O.L), tau(O.tau), rho(O.rho), freq(O.freq), order(O.order), scores(O.scores), patscore(O.patscore), ll(0),
      gradient(Eigen::VectorXd::Zero(O.gradient.size())), hrr(Eigen::VectorXd::Zero(O.hrr.size())),
      htt(Eigen::MatrixXd::Zero(O.htt.rows(), O.htt.cols())), htr(Eigen::MatrixXd::Zero(O.htr.rows(), O.htr.cols())),
      hschur(Eigen::MatrixXd::Zero(O.htt.rows(), O.htt.cols())){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const Saturated &O);
  };

  struct SaturatedCases : public RcppParallel::Worker{
    const model::Layout &L;
    const Eigen::Ref<const Eigen::MatrixXd> y, patscore, uscore;
    const Eigen::Ref<const Eigen::VectorXd> uw;
    Eigen::MatrixXd &st, &sr;

    SaturatedCases(const model::Layout &LAY, const Eigen::Ref<const Eigen::MatrixXd> Y, const Eigen::Ref<const Eigen::MatrixXd> PATSCORE,
                   const Eigen::Ref<const Eigen::MatrixXd> USCORE, const Eigen::Ref<const Eigen::VectorXd> UW,
                   Eigen::MatrixXd &ST, Eigen::MatrixXd &SR):
      L(LAY), y(Y), patscore(PATSCORE), uscore(USCORE), uw(UW), st(ST), sr(SR){}
    void operator()(std::size_t BEGIN, std::size_t END);
  };

}

#endif

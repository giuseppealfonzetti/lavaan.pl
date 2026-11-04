#ifndef workers_H
#define workers_H
#include <RcppEigen.h>
#include <RcppParallel.h>
#include "pairs.h"

namespace workers{

  struct PL : public RcppParallel::Worker{
    const model::Setup &st;
    const pairs::Data &data;
    const std::vector<int> &rs;
    const int order;
    double ll;
    Eigen::VectorXd gradient;

    PL(const model::Setup &ST, const pairs::Data &DATA, const std::vector<int> &RS, const int ORDER):
      st(ST), data(DATA), rs(RS), order(ORDER), ll(0), gradient(Eigen::VectorXd::Zero(ST.S.d)){}
    PL(const PL &O, RcppParallel::Split):
      st(O.st), data(O.data), rs(O.rs), order(O.order), ll(0), gradient(Eigen::VectorXd::Zero(O.st.S.d)){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const PL &O);
  };

  struct Hess : public RcppParallel::Worker{
    const model::Setup &st;
    const pairs::Data &data;
    std::vector<Eigen::MatrixXd> &score;
    std::vector<std::vector<int>> &sidx;
    std::vector<Eigen::VectorXd> &pgrad;
    Eigen::MatrixXd hessian;
    Eigen::VectorXd gradient;

    Hess(const model::Setup &ST, const pairs::Data &DATA, std::vector<Eigen::MatrixXd> &SCORE,
         std::vector<std::vector<int>> &SIDX, std::vector<Eigen::VectorXd> &PGRAD):
      st(ST), data(DATA), score(SCORE), sidx(SIDX), pgrad(PGRAD),
      hessian(Eigen::MatrixXd::Zero(ST.S.d, ST.S.d)), gradient(Eigen::VectorXd::Zero(ST.S.d)){}
    Hess(const Hess &O, RcppParallel::Split):
      st(O.st), data(O.data), score(O.score), sidx(O.sidx), pgrad(O.pgrad),
      hessian(Eigen::MatrixXd::Zero(O.st.S.d, O.st.S.d)), gradient(Eigen::VectorXd::Zero(O.st.S.d)){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const Hess &O);
  };

  struct Cases : public RcppParallel::Worker{
    const model::Setup &st;
    const pairs::Data &data;
    const std::vector<Eigen::MatrixXd> &score;
    const std::vector<std::vector<int>> &sidx;
    const Eigen::Ref<const Eigen::MatrixXd> uscore;
    const Eigen::Ref<const Eigen::VectorXd> uw;
    Eigen::MatrixXd &out;

    Cases(const model::Setup &ST, const pairs::Data &DATA, const std::vector<Eigen::MatrixXd> &SCORE,
          const std::vector<std::vector<int>> &SIDX, const Eigen::Ref<const Eigen::MatrixXd> USCORE,
          const Eigen::Ref<const Eigen::VectorXd> UW, Eigen::MatrixXd &CASES):
      st(ST), data(DATA), score(SCORE), sidx(SIDX), uscore(USCORE), uw(UW), out(CASES){}
    void operator()(std::size_t BEGIN, std::size_t END);
  };

  struct Saturated : public RcppParallel::Worker{
    const model::Layout &L;
    const pairs::Data &data;
    const pairs::Reduced &reduced;
    const Eigen::Ref<const Eigen::VectorXd> rho;
    const int order;
    const bool scores;
    Eigen::MatrixXd &patscore;
    double ll;
    Eigen::VectorXd gradient, hrr;
    Eigen::MatrixXd htt, htr, hschur;

    Saturated(const model::Layout &LAY, const pairs::Data &DATA, const pairs::Reduced &REDUCED, const Eigen::Ref<const Eigen::VectorXd> RHO,
              const int NI, const int ORDER, const bool SCORES, Eigen::MatrixXd &PATSCORE):
      L(LAY), data(DATA), reduced(REDUCED), rho(RHO), order(ORDER), scores(SCORES), patscore(PATSCORE), ll(0),
      gradient(Eigen::VectorXd::Zero(NI + RHO.size())), hrr(Eigen::VectorXd::Zero(RHO.size())),
      htt(Eigen::MatrixXd::Zero(NI, NI)), htr(Eigen::MatrixXd::Zero(NI, RHO.size())), hschur(Eigen::MatrixXd::Zero(NI, NI)){}
    Saturated(const Saturated &O, RcppParallel::Split):
      L(O.L), data(O.data), reduced(O.reduced), rho(O.rho), order(O.order), scores(O.scores), patscore(O.patscore), ll(0),
      gradient(Eigen::VectorXd::Zero(O.gradient.size())), hrr(Eigen::VectorXd::Zero(O.hrr.size())),
      htt(Eigen::MatrixXd::Zero(O.htt.rows(), O.htt.cols())), htr(Eigen::MatrixXd::Zero(O.htr.rows(), O.htr.cols())),
      hschur(Eigen::MatrixXd::Zero(O.htt.rows(), O.htt.cols())){}
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const Saturated &O);
  };

  struct SaturatedCases : public RcppParallel::Worker{
    const model::Layout &L;
    const pairs::Data &data;
    const pairs::Reduced &reduced;
    const Eigen::Ref<const Eigen::VectorXd> rho;
    const Eigen::Ref<const Eigen::MatrixXd> patscore, uscore;
    const Eigen::Ref<const Eigen::VectorXd> uw;
    Eigen::MatrixXd &st, &sr;

    SaturatedCases(const model::Layout &LAY, const pairs::Data &DATA, const pairs::Reduced &REDUCED,
                   const Eigen::Ref<const Eigen::VectorXd> RHO, const Eigen::Ref<const Eigen::MatrixXd> PATSCORE,
                   const Eigen::Ref<const Eigen::MatrixXd> USCORE, const Eigen::Ref<const Eigen::VectorXd> UW, Eigen::MatrixXd &ST,
                   Eigen::MatrixXd &SR):
      L(LAY), data(DATA), reduced(REDUCED), rho(RHO), patscore(PATSCORE), uscore(USCORE), uw(UW), st(ST), sr(SR){}
    void operator()(std::size_t BEGIN, std::size_t END);
  };

}

#endif

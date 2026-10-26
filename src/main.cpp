#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <random>
#include <math.h>

// [[Rcpp::depends(RcppEigen, RcppParallel)]]
#include "genericUtils.h"
#include "frequencies.h"
#include "pairs.h"
#include "optimisationUtils.h"
#include "fullPairwise.h"
#include "exportedFuns.h"
#include "parallelWorkers.h"


//' Full pairwise iteration
//'
//' @description
//' Evaluate negative log pairwise likelihood or gradient of the complete pool of pairs.
//' Used by external optimisers. Multithreading options via `RcppParallel`.
//'
//' @param N Number of observations
//' @param C_VEC Vector containing the number of categories for each item
//' @param CONSTRMAT \eqn{p*q}-dimensional matrix. Elements set to `NA` refers to free loading parameters. Elements set to numerical values denote fixed values constraints.
//' @param CONSTRLOGSD \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent log standard deviations parameters. Elements set to numerical values denote fixed values constraints.
//' @param THETA Parameter vector
//' @param FREQ Frequency table
//' @param CORRFLAG TRUE to estimate latent correlations. 0 for orthogonal latent factors.
//' @param NTHR Number of thresholds parameters.
//' @param NLOAD Number of free loadings parameters
//' @param NCORR Number of free latent correlations parameters.
//' @param NVAR Number of free latent variance parameters.
//' @param GRFLAG 0 to only compute the likelihood. 1 to also compute the gradient.
//' @param SILENTFLAG optional for verbose output
//' @param UFREQ Univariate frequencies of the available-cases term. Empty for none.
//' @param REGMAT \eqn{q*q}-dimensional matrix of structural regressions. Elements set to `NA` refer to free coefficients. Empty for none.
//' @param RESCOV \eqn{p*p}-dimensional matrix of residual covariances. Off-diagonal elements set to `NA` refer to free covariances. Empty for none.
//'
// [[Rcpp::export]]
Rcpp::List cpp_multiThread_completePairwise(
    const int N,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    const int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int GRFLAG = 0,
    const int SILENTFLAG = 1,
    Rcpp::NumericVector UFREQ = Rcpp::NumericVector::create(),
    Rcpp::Nullable<Rcpp::NumericMatrix> REGMAT = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> RESCOV = R_NilValue
){
  const Eigen::Map<Eigen::VectorXd> ufreq(UFREQ.begin(), UFREQ.size());


  // output list
  Rcpp::List output = fullPairwise::multiThreadContribution(
    N,
    C_VEC,
    CONSTRMAT,
    CONSTRLOGSD,
    THETA,
    FREQ,
    CORRFLAG,
    NTHR,
    NLOAD,
    NCORR,
    NVAR,
    GRFLAG,
    SILENTFLAG,
    ufreq,
    params::nullable_mat(REGMAT),
    params::nullable_mat(RESCOV)
  );

  return(output);
}


// [[Rcpp::export]]
Rcpp::List cpp_plSA2(
    Eigen::Map<Eigen::MatrixXd> FREQ,
    Eigen::Map<Eigen::MatrixXd> VALFREQ,
    const int N,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    Eigen::Map<Eigen::VectorXd> THETA_INIT,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int PAIRS_PER_ITERATION,
    const double STEP0,
    const double STEP1,
    const double STEP2,
    const double STEP3,
    const int BURNE,
    const int MAXE,
    const int EPOCHS_PER_CHECK,
    const int MAX_TOL_COUNTER,
    const double TOL_BURN,
    const double TOL_END,
    const int CHECK_TOL,
    const int SEED,
    const bool VERBOSE,
    bool VERBOSE_ITER=false,
    Rcpp::NumericVector UFREQ = Rcpp::NumericVector::create(),
    Rcpp::NumericVector VALUFREQ = Rcpp::NumericVector::create()
){
  const Eigen::Map<Eigen::VectorXd> ufreq(UFREQ.begin(), UFREQ.size());
  const Eigen::Map<Eigen::VectorXd> valufreq(VALUFREQ.begin(), VALUFREQ.size());
  Eigen::MatrixXd hdummy, udummy;

  if(!VERBOSE) VERBOSE_ITER=false;

  const int d = THETA_INIT.size();
  const int n = N;
  const int p = CONSTRMAT.rows();
  const int pairs = p*(p-1)/2;
  const int corrflag = NCORR>0;
  if(d!=NTHR+NLOAD+NCORR+NVAR) Rcpp::stop("check theta dimensions");
  if(PAIRS_PER_ITERATION>pairs)Rcpp::stop("too many pairs per iterations");

  int burnt              = 0;
  int burne              = BURNE;
  const double prob      = static_cast<double>(PAIRS_PER_ITERATION)/static_cast<double>(pairs);
  int tol_counter        = 0;
  int tol_counter_burn   = 0;

  double stepsize   = STEP0;

  // Read frequencies, and initialise items_pairs
  // Eigen::MatrixXd pairs_table = FREQ;
  Eigen::MatrixXd items_pairs(2,pairs);
  const Eigen::MatrixXd nostruct;
  int idx = 0;
  for(int k = 1; k < p; k++){
    for(int l = 0; l < k; l++){
      items_pairs(0, idx) = k;
      items_pairs(1, idx) = l;
      idx++;
    }
  }

  // Initialize vector of indexes for entries in pairs_table
  std::vector<int> full_pool(pairs) ;
  std::iota (std::begin(full_pool), std::end(full_pool), 0);

  // Initialize storage for iterations quantities
  Eigen::VectorXd theta      = THETA_INIT;
  Eigen::VectorXd avtheta    = THETA_INIT;
  Eigen::VectorXd thetapdiff = Eigen::VectorXd::Zero(d);


  // double          nll   = std::numeric_limits<double>::infinity();
  std::vector<Eigen::VectorXd> path_theta;
  std::vector<Eigen::VectorXd> path_avtheta;
  std::vector<Eigen::VectorXd> path_thetapdiff;
  std::vector<double>          path_nll;
  std::vector<int>             path_iters;
  std::vector<int>             post_index;



  ///////////////////////
  // COMPUTE START NLL  //
  ///////////////////////
  int maxtpe = pairs/PAIRS_PER_ITERATION;
  double nll = fullPairwise::multiThreadContribution(
    N,
    C_VEC,
    CONSTRMAT,
    CONSTRLOGSD,
    theta,
    VALFREQ,
    corrflag,
    NTHR,
    NLOAD,
    NCORR,
    NVAR,
    0,
    1,
    valufreq,
    nostruct,
    nostruct
  )["iter_nll"];
  nll/=N;
  if(VERBOSE) Rcpp::Rcout << "Initial full npll: "<< nll << " | iterations per epoch:"<<maxtpe << "\n";




  // Compute scaling constant
  double scale = prob/static_cast<double>(n);


  // Diagnostics
  bool neg_pdiff        = false;
  bool proj_after_burn  = false;
  bool finite_gr        = true;
  bool convergence_burn = true;
  int  convergence_full = 0;


  int last_iter = 0;
  double pdiff = 0;





  //////////////////////////////
  // STORE INITIAL QUANTITIES //
  /////////////////////////////
  Eigen::VectorXd prev_theta = theta;
  path_iters.push_back(0);
  path_theta.push_back(theta);
  path_avtheta.push_back(avtheta);
  path_thetapdiff.push_back(thetapdiff);
  path_nll.push_back(nll);


  //////////////////
  // OPTIMISATION //
  //////////////////
  int t=1;
  for(int epoch=0; epoch < MAXE; epoch++){
    Rcpp::checkUserInterrupt();

    Eigen::VectorXd prev_epoch_theta = theta;
    Eigen::VectorXd prev_epoch_avtheta = avtheta;



    /////////////////////
    // SAMPLING SCHEME //
    /////////////////////
    // Set-up the randomizer
    std::mt19937 randomizer(SEED + epoch);
    std::shuffle(full_pool.begin(), full_pool.end(), randomizer);
    int idx_1pair=0;

    /////////////////////
    // EPOCH UPDATES.  //
    /////////////////////
    for(int te = 0; te < maxtpe; te++){
      Rcpp::checkUserInterrupt();

      // Initialize empty contributions for iteration
      Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);

      ///////////////////////////
      // GRADIENT COMPUTATION  //
      ///////////////////////////
      std::vector<int> iter_chosen_pairs;
      for(int draw = 0; draw < PAIRS_PER_ITERATION; draw++){
        int pair_index = full_pool.at(idx_1pair+draw);
        iter_chosen_pairs.push_back(pair_index);
      }
      idx_1pair += PAIRS_PER_ITERATION;

      parallelWorkers::pairsGradient iteration_subset(CONSTRMAT, CONSTRLOGSD, C_VEC, FREQ, items_pairs, nostruct, nostruct, corrflag,
                                                      NTHR, NLOAD, NCORR, NVAR,
                                                      1, 1, theta, iter_chosen_pairs);
      RcppParallel::parallelReduce(0, iter_chosen_pairs.size(), iteration_subset);

      iter_gradient = -iteration_subset.subset_gradient;
      iter_gradient *= scale;
      if(ufreq.size() > 0){
        double ull = 0;
        Eigen::VectorXd ugrad = Eigen::VectorXd::Zero(d);
        pairs::univariate_contribution(theta.head(NTHR), C_VEC, ufreq, 1, ull, ugrad, hdummy, udummy);
        iter_gradient -= scale * prob * ugrad;
      }


      //////////////////////////
      //   PARAMETERS UPDATE  //
      //////////////////////////
      bool update = true;
      int attempt_counter = 0;
      double step_length;
      while(update && attempt_counter <10){
        step_length = stepsize;
        step_length *= STEP1 * pow(1 + STEP2*stepsize*(t), -STEP3);



        // iter_gradient.segment(0, NTHR) /=10;
        iter_gradient.segment(NTHR, NLOAD) *= 5;
        if(NCORR>0) iter_gradient.segment(NTHR+NLOAD, NCORR) *= 100;
        if(NVAR>0)  iter_gradient.segment(NTHR+NLOAD+NCORR, NVAR) *= 5;



        if(std::isfinite(iter_gradient.norm())){
          update=false;
          prev_theta = theta;
          theta -= step_length * iter_gradient;
          }else{
          stepsize*=.75;
          theta = prev_theta;
          finite_gr = false;
          }
        attempt_counter++;
      }



      ////////////////////////////////////////////////
      // Post update loadings projection, if needed //
      ///////////////////////////////////////////////

      bool checkevent = false;
      sa::proj2(CONSTRMAT, CONSTRLOGSD, C_VEC, corrflag, NTHR, NLOAD, NCORR, NVAR, theta, checkevent);
      if(checkevent) post_index.push_back(t);
      if(checkevent && epoch>burne) proj_after_burn=true;

      thetapdiff = (theta.array()-prev_theta.array())/prev_theta.array();
      if(VERBOSE_ITER) Rcpp::Rcout << " - Iter: " << t << " | step: " << step_length << " | mean abs theta pdiff: "<< thetapdiff.array().abs().mean()<<
        " | gradient norm: "<< iter_gradient.norm()<< "\n";



      if(epoch <= burne){
        avtheta = theta;
      }else{
        avtheta = ( (t - burnt) * avtheta + theta ) / (t - burnt + 1);
      }


      // theta= theta_prev-ste*grad
      // theta-theta_prev= theta_prev-ste*grad-theta_prev=-ste*grad











    t++;}





    ///////////////////////
    // COMPUTE FULL NLL  //
    ///////////////////////
    double prev_nll = nll;
    if(epoch%EPOCHS_PER_CHECK==0){
      parallelWorkers::pairsGradient fullpool_worker(CONSTRMAT, CONSTRLOGSD, C_VEC, VALFREQ, items_pairs, nostruct, nostruct, corrflag, NTHR, NLOAD, NCORR, NVAR,
                                                    1, 0, theta, full_pool);
      RcppParallel::parallelReduce(0, pairs, fullpool_worker);
      nll = -fullpool_worker.subset_ll/n;
      double ull = 0;
      Eigen::VectorXd ugrad;
      pairs::univariate_contribution(theta.head(NTHR), C_VEC, valufreq, 0, ull, ugrad, hdummy, udummy);
      nll -= ull/n;
    }



    ///////////////////
    // STORE EPOCH   //
    ///////////////////
    if(epoch%EPOCHS_PER_CHECK==0){
      path_iters.push_back(t);
      path_theta.push_back(theta);
      path_avtheta.push_back(avtheta);
      path_thetapdiff.push_back(thetapdiff);
      path_nll.push_back(nll);
    }



    ///////////////////////
    // MONITOR FULL NLL //
    ///////////////////////
    // Check convergence
    if(epoch%EPOCHS_PER_CHECK==0){
      pdiff = - (nll - prev_nll)/(prev_nll*maxtpe);
      if(pdiff<TOL_END){tol_counter++;}else{tol_counter=0;}
      if(pdiff<(TOL_BURN)){tol_counter_burn++;}else{tol_counter_burn=0;}
      if((pdiff<0) && (epoch<=burne)) neg_pdiff=true;
    }


    if((t>0) && !std::isfinite(nll)){
      if(VERBOSE) Rcpp::Rcout << "Exiting because of non-finite nll\n";
      last_iter       = t;
      convergence_full= -1;
      break;
    }

    ////////////
    // REPORT //
    ////////////
    Eigen::VectorXd prev_epoch_thetapdiff = (theta.array()-prev_epoch_theta.array())/prev_epoch_theta.array();
    if(VERBOSE && (epoch%EPOCHS_PER_CHECK==0)) Rcpp::Rcout << "End of Epoch: "<<epoch<< " | full npll: "<< nll << " | av nll pdiff from prev epoch: " << pdiff << " | tc:" << tol_counter<<
      "| mean abs theta pdiff from prev epoch: "<< prev_epoch_thetapdiff.array().abs().mean()<<"\n";

    /////////////////
    // AUTO BURNIN //
    /////////////////
    if(epoch==burne){
      if(VERBOSE) Rcpp::Rcout << "Burn-in ended\n";
      convergence_burn=false;
      burnt=t;
    }else if((CHECK_TOL && (tol_counter_burn >= MAX_TOL_COUNTER) && (epoch <= burne) && (epoch%EPOCHS_PER_CHECK==0))){
      if(VERBOSE) Rcpp::Rcout << "Burn-in stopped\n";
        burne = epoch;
        tol_counter_burn=0;
        convergence_burn=true;
        burnt =t;

    }


    /////////////
    // EXITING //
    /////////////
    if(CHECK_TOL && (tol_counter>=MAX_TOL_COUNTER) && ( epoch>(burne) && (epoch%EPOCHS_PER_CHECK==0))){
      if(VERBOSE) Rcpp::Rcout << "Converged\n";
        convergence_full = 1;
        last_iter = t;
        if(pdiff<0){
          if(VERBOSE) Rcpp::Rcout << "Reporting parameter values from epoch " << epoch-1<< "\n";
          theta=prev_epoch_theta;
          avtheta=prev_epoch_avtheta;
        }
        break;
    }else if(epoch==(MAXE-1)){
      if(VERBOSE) Rcpp::Rcout << "Ended\n";
      last_iter = t;
    }



  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("burnt") = burnt,
      Rcpp::Named("scale") = scale,
      Rcpp::Named("nll")   = nll,
      Rcpp::Named("path_theta") = path_theta,
      Rcpp::Named("path_thetapdiff") = path_thetapdiff,
      Rcpp::Named("path_avtheta") = path_avtheta,
      Rcpp::Named("path_iters") = path_iters,
      Rcpp::Named("path_nll") = path_nll,
      Rcpp::Named("post_index") = post_index,
      Rcpp::Named("last_iter") = last_iter,
      Rcpp::Named("theta") = theta,
      Rcpp::Named("avtheta") = avtheta,
      Rcpp::Named("convergence_full") = convergence_full,
      Rcpp::Named("convergence_burn") = convergence_burn,
      Rcpp::Named("proj_after_burn") = proj_after_burn,
      Rcpp::Named("finite_gr") = finite_gr,
      Rcpp::Named("neg_pdiff") = neg_pdiff
    );
  return(output);


}


























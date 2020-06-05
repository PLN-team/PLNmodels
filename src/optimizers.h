#ifndef _optimizers_H
#define _optimizers_H

#include "nlopt_utils.h"
#include "gradients.h"

// ---------------------------------------------------------------------------
// ABSTRACT CLASS OPTIMIZER_PLN
//
// COMMON TO PLN WITH SPHERICAL, DIAGONAL, FULLY PARAMETRIZED COVARIANCES
// MOTHER OF PLN WITH COVARIANCES HAVING RANK AND SPARSITY CONSTRAINTS
class optimizer_PLN {
  protected:

    // problem dimension
    int n;
    int p;
    int d;

    // structure use to stock data (count, covariates, etc.) used during the optimization process
    optim_data data  ;

    // value of vector of parameters at optimum
    stdvec parameter ;

    // the function that computes the objective and thge gradient vector
    // double (*fn_optim) (const stdvec& , stdvec &, void *) ;
    double (*fn_optim) (unsigned , const double* , double* , void*) ;

    // the function that computes the objective and thge gradient vector
    // double (*fn_optim_VE) (const stdvec& , stdvec &, void *) ;
    double (*fn_VEstep) (unsigned , const double* , double* , void*) ;

    // nlopt optimizer
    nlopt_opt optimizer ;

    // result of the optimization process
    nlopt_result status ;

    // element-wise log-likelihood
    arma::vec loglik ;

    // matrices of parameters
    arma::mat Theta ;
    arma::mat Sigma ;
    arma::mat Omega ;
    arma::mat M     ;
    arma::mat S     ;
    arma::mat Z     ;
    arma::mat A     ;

  public:

    // Constructors
    optimizer_PLN() {} ;

    optimizer_PLN(
          arma::vec par,
          const arma::mat & Y,
          const arma::mat & X,
          const arma::mat & O,
          const arma::vec & w,
          Rcpp::List options
    ) ;

    // list encompasing the results of the optimization process sent back to R
    Rcpp::List output ;

    void optimize() ;

    void VEstep(const arma::mat &, const arma::mat &) ;

    // prepare/compute output according to problem dimension
    // will be defined in the child classes
    virtual void export_output() =0 ;

    // export the output an Rcpp::List understandable by R
    virtual Rcpp::List get_output() ;

    // export the output an Rcpp::List understandable by R
    virtual Rcpp::List get_var_par() ;

};

// ---------------------------------------------------------------------------
// CHILD CLASSES of OPTIMIZER_PLN
//

// SPHERICAL COVARIANCE
class optimizer_PLN_spherical: public optimizer_PLN {
  public:
    optimizer_PLN_spherical(
      arma::vec par,
      const arma::mat & Y,
      const arma::mat & X,
      const arma::mat & O,
      const arma::vec & w,
      Rcpp::List options
    ) ;

    void export_output()  ;
    void export_var_par() ;
};

// DIAGONAL COVARIANCE
class optimizer_PLN_diagonal: public optimizer_PLN {
  public:
    optimizer_PLN_diagonal(
      arma::vec par,
      const arma::mat & Y,
      const arma::mat & X,
      const arma::mat & O,
      const arma::vec & w,
      Rcpp::List options
    ) ;

    void export_output()  ;
    void export_var_par() ;
};

// FULLY PARAMETRIZED COVARIANCE
class optimizer_PLN_full: public optimizer_PLN {
  public:
    optimizer_PLN_full(
      arma::vec par,
      const arma::mat & Y,
      const arma::mat & X,
      const arma::mat & O,
      const arma::vec & w,
      Rcpp::List options
    ) ;

    void export_output()  ;
    void export_var_par() ;
};

// RANK-CONSTRAINED COVARIANCE (PCA)
class optimizer_PLN_rank: public optimizer_PLN {
  public:
    optimizer_PLN_rank(
      arma::vec par,
      const arma::mat & Y,
      const arma::mat & X,
      const arma::mat & O,
      const arma::vec & w,
      const int rank,
      Rcpp::List options
    ) ;

    void export_output() ;

    // override Mother's method definition (one additional parameter to send back)
    Rcpp::List get_output() ;

  protected:

    // Rank
    int q ;

    // matrix of scores
    arma::mat B ;
};

// SPARSE INVERSE COVARIANCE (PCA)
class optimizer_PLN_sparse: public optimizer_PLN {
  public:
    optimizer_PLN_sparse(
      arma::vec par,
      const arma::mat & Y,
      const arma::mat & X,
      const arma::mat & O,
      const arma::vec & w,
      const arma::mat & Omega,
      Rcpp::List options
    ) ;

    void export_output() ;
};

#endif

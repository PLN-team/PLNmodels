#ifndef _gradients_H
#define _gradients_H

#include "utils.h"
#include "RcppArmadillo.h"

double fn_optim_PLN(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_weighted(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_spherical(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_weighted_spherical(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_diagonal(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_weighted_diagonal(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_rank(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_weighted_rank(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_sparse(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_weighted_sparse(const stdvec &x, stdvec &grad, void *data) ;

#endif

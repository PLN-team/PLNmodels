#ifndef _gradients_PLN_H
#define _gradients_PLN_H

#include "utils_optim.h"
#include "RcppArmadillo.h"

double fn_optim_PLN(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_weighted(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_spherical(const stdvec &x, stdvec &grad, void *data) ;

double fn_optim_PLN_weighted_spherical(const stdvec &x, stdvec &grad, void *data) ;

#endif

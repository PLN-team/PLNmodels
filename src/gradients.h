#ifndef _gradients_H
#define _gradients_H

#include "data_struct.h"

double fn_optim_PLN(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_weighted(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_spherical(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_weighted_spherical(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_diagonal(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_weighted_diagonal(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_rank(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_weighted_rank(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_sparse(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_PLN_weighted_sparse(unsigned N, const double *x, double *grad, void *data) ;

double fn_optim_VEstep_PLN(unsigned N, const double *x, double *grad, void *data) ;

#endif

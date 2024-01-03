// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_test_nlopt
bool cpp_test_nlopt();
RcppExport SEXP _PLNmodels_cpp_test_nlopt() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cpp_test_nlopt());
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_diagonal
Rcpp::List nlopt_optimize_diagonal(const Rcpp::List& data, const Rcpp::List& params, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_diagonal(SEXP dataSEXP, SEXP paramsSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_diagonal(data, params, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_vestep_diagonal
Rcpp::List nlopt_optimize_vestep_diagonal(const Rcpp::List& data, const Rcpp::List& params, const arma::mat& B, const arma::mat& Omega, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_vestep_diagonal(SEXP dataSEXP, SEXP paramsSEXP, SEXP BSEXP, SEXP OmegaSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_vestep_diagonal(data, params, B, Omega, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_fixed
Rcpp::List nlopt_optimize_fixed(const Rcpp::List& data, const Rcpp::List& params, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_fixed(SEXP dataSEXP, SEXP paramsSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_fixed(data, params, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize
Rcpp::List nlopt_optimize(const Rcpp::List& data, const Rcpp::List& params, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize(SEXP dataSEXP, SEXP paramsSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize(data, params, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_vestep
Rcpp::List nlopt_optimize_vestep(const Rcpp::List& data, const Rcpp::List& params, const arma::mat& B, const arma::mat& Omega, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_vestep(SEXP dataSEXP, SEXP paramsSEXP, SEXP BSEXP, SEXP OmegaSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_vestep(data, params, B, Omega, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_genetic_modeling
Rcpp::List nlopt_optimize_genetic_modeling(const Rcpp::List& init_parameters, const arma::mat& Y, const arma::mat& X, const arma::mat& O, const arma::vec& w, const arma::mat& C, const Rcpp::List& configuration);
RcppExport SEXP _PLNmodels_nlopt_optimize_genetic_modeling(SEXP init_parametersSEXP, SEXP YSEXP, SEXP XSEXP, SEXP OSEXP, SEXP wSEXP, SEXP CSEXP, SEXP configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type init_parameters(init_parametersSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type O(OSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type configuration(configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_genetic_modeling(init_parameters, Y, X, O, w, C, configuration));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_rank
Rcpp::List nlopt_optimize_rank(const Rcpp::List& data, const Rcpp::List& params, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_rank(SEXP dataSEXP, SEXP paramsSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_rank(data, params, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_vestep_rank
Rcpp::List nlopt_optimize_vestep_rank(const Rcpp::List& data, const Rcpp::List& params, const arma::mat& B, const arma::mat& C, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_vestep_rank(SEXP dataSEXP, SEXP paramsSEXP, SEXP BSEXP, SEXP CSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_vestep_rank(data, params, B, C, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_spherical
Rcpp::List nlopt_optimize_spherical(const Rcpp::List& data, const Rcpp::List& params, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_spherical(SEXP dataSEXP, SEXP paramsSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_spherical(data, params, config));
    return rcpp_result_gen;
END_RCPP
}
// nlopt_optimize_vestep_spherical
Rcpp::List nlopt_optimize_vestep_spherical(const Rcpp::List& data, const Rcpp::List& params, const arma::mat& B, const arma::mat& Omega, const Rcpp::List& config);
RcppExport SEXP _PLNmodels_nlopt_optimize_vestep_spherical(SEXP dataSEXP, SEXP paramsSEXP, SEXP BSEXP, SEXP OmegaSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(nlopt_optimize_vestep_spherical(data, params, B, Omega, config));
    return rcpp_result_gen;
END_RCPP
}
// zipln_vloglik
arma::vec zipln_vloglik(const arma::mat& Y, const arma::mat& X, const arma::mat& O, const arma::mat& Pi, const arma::mat& Omega, const arma::mat& B, const arma::mat& R, const arma::mat& M, const arma::mat& S);
RcppExport SEXP _PLNmodels_zipln_vloglik(SEXP YSEXP, SEXP XSEXP, SEXP OSEXP, SEXP PiSEXP, SEXP OmegaSEXP, SEXP BSEXP, SEXP RSEXP, SEXP MSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type O(OSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(zipln_vloglik(Y, X, O, Pi, Omega, B, R, M, S));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_Omega_full
arma::mat optim_zipln_Omega_full(const arma::mat& M, const arma::mat& X, const arma::mat& B, const arma::mat& S);
RcppExport SEXP _PLNmodels_optim_zipln_Omega_full(SEXP MSEXP, SEXP XSEXP, SEXP BSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_Omega_full(M, X, B, S));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_Omega_spherical
arma::mat optim_zipln_Omega_spherical(const arma::mat& M, const arma::mat& X, const arma::mat& B, const arma::mat& S);
RcppExport SEXP _PLNmodels_optim_zipln_Omega_spherical(SEXP MSEXP, SEXP XSEXP, SEXP BSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_Omega_spherical(M, X, B, S));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_Omega_diagonal
arma::mat optim_zipln_Omega_diagonal(const arma::mat& M, const arma::mat& X, const arma::mat& B, const arma::mat& S);
RcppExport SEXP _PLNmodels_optim_zipln_Omega_diagonal(SEXP MSEXP, SEXP XSEXP, SEXP BSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_Omega_diagonal(M, X, B, S));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_B_dense
arma::mat optim_zipln_B_dense(const arma::mat& M, const arma::mat& X);
RcppExport SEXP _PLNmodels_optim_zipln_B_dense(SEXP MSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_B_dense(M, X));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_zipar_covar
Rcpp::List optim_zipln_zipar_covar(const arma::mat& init_B0, const arma::mat& X, const arma::mat& R, const Rcpp::List& configuration);
RcppExport SEXP _PLNmodels_optim_zipln_zipar_covar(SEXP init_B0SEXP, SEXP XSEXP, SEXP RSEXP, SEXP configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type init_B0(init_B0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type configuration(configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_zipar_covar(init_B0, X, R, configuration));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_R
arma::mat optim_zipln_R(const arma::mat& Y, const arma::mat& X, const arma::mat& O, const arma::mat& M, const arma::mat& S, const arma::mat& Pi);
RcppExport SEXP _PLNmodels_optim_zipln_R(SEXP YSEXP, SEXP XSEXP, SEXP OSEXP, SEXP MSEXP, SEXP SSEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type O(OSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_R(Y, X, O, M, S, Pi));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_M
Rcpp::List optim_zipln_M(const arma::mat& init_M, const arma::mat& Y, const arma::mat& X, const arma::mat& O, const arma::mat& R, const arma::mat& S, const arma::mat& B, const arma::mat& Omega, const Rcpp::List& configuration);
RcppExport SEXP _PLNmodels_optim_zipln_M(SEXP init_MSEXP, SEXP YSEXP, SEXP XSEXP, SEXP OSEXP, SEXP RSEXP, SEXP SSEXP, SEXP BSEXP, SEXP OmegaSEXP, SEXP configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type init_M(init_MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type O(OSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type configuration(configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_M(init_M, Y, X, O, R, S, B, Omega, configuration));
    return rcpp_result_gen;
END_RCPP
}
// optim_zipln_S
Rcpp::List optim_zipln_S(const arma::mat& init_S, const arma::mat& O, const arma::mat& M, const arma::mat& R, const arma::mat& B, const arma::vec& diag_Omega, const Rcpp::List& configuration);
RcppExport SEXP _PLNmodels_optim_zipln_S(SEXP init_SSEXP, SEXP OSEXP, SEXP MSEXP, SEXP RSEXP, SEXP BSEXP, SEXP diag_OmegaSEXP, SEXP configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type init_S(init_SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type O(OSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type diag_Omega(diag_OmegaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type configuration(configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_zipln_S(init_S, O, M, R, B, diag_Omega, configuration));
    return rcpp_result_gen;
END_RCPP
}
// cpp_test_packing
bool cpp_test_packing();
RcppExport SEXP _PLNmodels_cpp_test_packing() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cpp_test_packing());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PLNmodels_cpp_test_nlopt", (DL_FUNC) &_PLNmodels_cpp_test_nlopt, 0},
    {"_PLNmodels_nlopt_optimize_diagonal", (DL_FUNC) &_PLNmodels_nlopt_optimize_diagonal, 3},
    {"_PLNmodels_nlopt_optimize_vestep_diagonal", (DL_FUNC) &_PLNmodels_nlopt_optimize_vestep_diagonal, 5},
    {"_PLNmodels_nlopt_optimize_fixed", (DL_FUNC) &_PLNmodels_nlopt_optimize_fixed, 3},
    {"_PLNmodels_nlopt_optimize", (DL_FUNC) &_PLNmodels_nlopt_optimize, 3},
    {"_PLNmodels_nlopt_optimize_vestep", (DL_FUNC) &_PLNmodels_nlopt_optimize_vestep, 5},
    {"_PLNmodels_nlopt_optimize_genetic_modeling", (DL_FUNC) &_PLNmodels_nlopt_optimize_genetic_modeling, 7},
    {"_PLNmodels_nlopt_optimize_rank", (DL_FUNC) &_PLNmodels_nlopt_optimize_rank, 3},
    {"_PLNmodels_nlopt_optimize_vestep_rank", (DL_FUNC) &_PLNmodels_nlopt_optimize_vestep_rank, 5},
    {"_PLNmodels_nlopt_optimize_spherical", (DL_FUNC) &_PLNmodels_nlopt_optimize_spherical, 3},
    {"_PLNmodels_nlopt_optimize_vestep_spherical", (DL_FUNC) &_PLNmodels_nlopt_optimize_vestep_spherical, 5},
    {"_PLNmodels_zipln_vloglik", (DL_FUNC) &_PLNmodels_zipln_vloglik, 9},
    {"_PLNmodels_optim_zipln_Omega_full", (DL_FUNC) &_PLNmodels_optim_zipln_Omega_full, 4},
    {"_PLNmodels_optim_zipln_Omega_spherical", (DL_FUNC) &_PLNmodels_optim_zipln_Omega_spherical, 4},
    {"_PLNmodels_optim_zipln_Omega_diagonal", (DL_FUNC) &_PLNmodels_optim_zipln_Omega_diagonal, 4},
    {"_PLNmodels_optim_zipln_B_dense", (DL_FUNC) &_PLNmodels_optim_zipln_B_dense, 2},
    {"_PLNmodels_optim_zipln_zipar_covar", (DL_FUNC) &_PLNmodels_optim_zipln_zipar_covar, 4},
    {"_PLNmodels_optim_zipln_R", (DL_FUNC) &_PLNmodels_optim_zipln_R, 6},
    {"_PLNmodels_optim_zipln_M", (DL_FUNC) &_PLNmodels_optim_zipln_M, 9},
    {"_PLNmodels_optim_zipln_S", (DL_FUNC) &_PLNmodels_optim_zipln_S, 7},
    {"_PLNmodels_cpp_test_packing", (DL_FUNC) &_PLNmodels_cpp_test_packing, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_PLNmodels(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

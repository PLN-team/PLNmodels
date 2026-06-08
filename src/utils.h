#pragma once

#include <RcppArmadillo.h>

inline arma::vec logfact(arma::mat y) {
    y.replace(0., 1.);
    return sum(y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2., 1);
}

inline arma::mat logfact_mat(arma::mat y) {
    y.replace(0., 1.);
    return y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2.;
}

inline arma::vec ki(arma::mat y) {
    arma::uword p = y.n_cols;
    return -logfact(std::move(y)) + 0.5 * double(p) ;
}

inline arma::mat logistic(arma::mat M) {
  arma::mat e = arma::trunc_exp(M);
  return e / (1. + e);
}

inline arma::mat logit(arma::mat M) {
  return arma::trunc_log(M) - arma::trunc_log(1 - M) ;
}

// ---- Newton step for B: diagonal Hessian approximation with Armijo line search ----
// Updates B, Z, A in-place (Z = O + X*B + M and A = exp(Z + S²/2) after update).
// Requires Xw = X.*w and Xw2 = X².*w precomputed outside the loop.
inline void newton_step_B(
    const arma::mat & Xw, const arma::mat & Xw2,
    const arma::mat & X,  const arma::mat & Y,
    const arma::mat & O,  const arma::vec & w,
    const arma::mat & M,  const arma::mat & S2,
    arma::mat & B, arma::mat & Z, arma::mat & A
) {
    constexpr double c1 = 1e-4;
    arma::mat grad_B = Xw.t() * (A - Y);
    arma::mat hess_B = Xw2.t() * A;
    hess_B.clamp(1e-10, arma::datum::inf);
    const arma::mat step_B  = grad_B / hess_B;
    const arma::mat XstepB  = X * step_B;
    const double f0_B    = arma::accu(w.t() * (A - Y % Z));
    const double slope_B = -arma::accu(grad_B % step_B);
    double alpha_B = 1.0;
    for (int ls = 0; ls < 20; ls++) {
        const arma::mat Zt = Z - alpha_B * XstepB;
        if (arma::accu(w.t() * (arma::exp(Zt + 0.5 * S2) - Y % Zt))
            <= f0_B + c1 * alpha_B * slope_B) break;
        alpha_B *= 0.5;
    }
    B -= alpha_B * step_B;
    Z  = O + X * B + M;
    A  = arma::exp(Z + 0.5 * S2);
}

// ---- Fixed-point update for logS (overflow-safe) ----
// cov_diag: diagonal of Omega broadcast to (n,p) — arma::mat or double (scalar broadcast).
// Updates A, logS, S, S2 in-place using logS = clamp(min(-½log(A+cov_diag), ½log(700-Z))).
template<typename CovDiagType>
inline void fixed_point_logS(
    arma::mat & logS, arma::mat & S, arma::mat & S2,
    const arma::mat & Z, arma::mat & A,
    const CovDiagType & cov_diag
) {
    A = arma::exp(Z + 0.5 * S2);
    const arma::mat logS_cand = -0.5 * arma::log(A + cov_diag);
    const arma::mat logS_ub   = 0.5 * arma::log(arma::clamp(700. - Z, 1., arma::datum::inf));
    logS = arma::clamp(arma::min(logS_cand, logS_ub), -20., arma::datum::inf);
    S    = arma::exp(logS);
    S2   = S % S;
}

// ---- Relative convergence test: |val - prev| < tol * (1 + |prev|) ----
inline bool converged(double val, double prev, double tol) {
    return std::abs(val - prev) < tol * (1.0 + std::abs(prev));
}

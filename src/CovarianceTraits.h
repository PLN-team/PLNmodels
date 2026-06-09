#pragma once
#include <RcppArmadillo.h>
#include "utils.h"

// ─────────────────────────────────────────────────────────────────────────────
// Full (dense p×p) covariance
// ─────────────────────────────────────────────────────────────────────────────
struct FullCovTraits {
    struct State {
        arma::mat Omega;
        arma::mat Sigma;
        arma::vec diag_Omega;

        State(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            Sigma      = (1./w_bar) * (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2));
            Omega      = arma::inv_sympd(Sigma);
            diag_Omega = arma::diagvec(Omega);
        }
        // vestep: Omega known and fixed
        explicit State(const arma::mat & omega)
            : Omega(omega), diag_Omega(arma::diagvec(omega)) {}
    };

    static arma::mat cov_diag(const State & s, const arma::mat & ones_row) {
        return ones_row * s.diag_Omega.t();
    }

    static void grad_hess_M(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & Y, const arma::vec & w, const arma::mat & ones_row,
        arma::mat & grad_M, arma::mat & hess_M)
    {
        arma::mat MO = M * s.Omega;
        grad_M = MO + A - Y; grad_M.each_col() %= w;
        hess_M = A + ones_row * s.diag_Omega.t(); hess_M.each_col() %= w;
    }

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return M * s.Omega; }

    // Takes precomputed MO = M * Omega to avoid redundant matrix multiplies in Armijo
    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    static double objective_cov(const arma::mat & M, const arma::mat & S2, const State & s, const arma::vec & w) {
        arma::mat MO = M * s.Omega;
        return 0.5 * (arma::as_scalar(w.t() * arma::sum(MO % M, 1))
                    + arma::dot(s.diag_Omega, (w.t() * S2).t()));
    }

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.Sigma      = (1./w_bar) * (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2));
        s.Omega      = arma::inv_sympd(s.Sigma);
        s.diag_Omega = arma::diagvec(s.Omega);
    }

    static double elbo_cov(const State & s, double w_bar, arma::uword /*p*/) {
        return -0.5 * w_bar * std::real(arma::log_det(s.Sigma));
    }

    static arma::vec final_loglik(const arma::mat & Y, const arma::mat & Z, const arma::mat & A,
                                   const arma::mat & M, const arma::mat & psi, const State & s) {
        const arma::mat S2 = arma::exp(psi);
        return arma::sum(Y % Z - A + 0.5 * psi
                       - 0.5 * ((M * s.Omega) % M + S2 * arma::diagmat(s.Omega)), 1)
             + 0.5 * std::real(arma::log_det(s.Omega)) + ki(Y);
    }

    static Rcpp::List output_cov(const arma::mat & /*M*/, const arma::mat & /*S2*/,
                                  const arma::vec & /*w*/, double /*w_bar*/, const State & s) {
        return Rcpp::List::create(Rcpp::Named("Sigma", s.Sigma), Rcpp::Named("Omega", s.Omega));
    }

    static constexpr bool has_em = true;
};

// ─────────────────────────────────────────────────────────────────────────────
// Diagonal covariance
// ─────────────────────────────────────────────────────────────────────────────
struct DiagonalCovTraits {
    struct State {
        arma::rowvec omega2;
        arma::rowvec sigma2;

        State(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            sigma2 = (w.t() * (M % M + S2)) / w_bar;
            omega2 = arma::pow(sigma2, -1);
        }
        // vestep: Omega known and fixed (diagonal matrix passed as dense mat)
        explicit State(const arma::mat & omega_mat) {
            omega2 = arma::diagvec(omega_mat).t();
            sigma2 = arma::pow(omega2, -1);
        }
    };

    static arma::mat cov_diag(const State & s, const arma::mat & ones_row) {
        return ones_row * s.omega2;
    }

    static void grad_hess_M(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & Y, const arma::vec & w, const arma::mat & ones_row,
        arma::mat & grad_M, arma::mat & hess_M)
    {
        grad_M = M.each_row() % s.omega2 + A - Y; grad_M.each_col() %= w;
        hess_M = ones_row * s.omega2 + A;          hess_M.each_col() %= w;
    }

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return M.each_row() % s.omega2; }

    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    static double objective_cov(const arma::mat & M, const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * arma::as_scalar((w.t() * (M % M + S2)) * s.omega2.t());
    }

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.sigma2 = (w.t() * (M % M + S2)) / w_bar;
        s.omega2 = arma::pow(s.sigma2, -1);
    }

    static double elbo_cov(const State & s, double w_bar, arma::uword /*p*/) {
        return -0.5 * w_bar * arma::accu(arma::log(s.sigma2));
    }

    static arma::vec final_loglik(const arma::mat & Y, const arma::mat & Z, const arma::mat & A,
                                   const arma::mat & M, const arma::mat & psi, const State & s) {
        const arma::mat S2    = arma::exp(psi);
        const arma::vec omega2_v = s.omega2.t();
        return arma::sum(Y % Z - A + 0.5 * psi, 1)
             - 0.5 * (M % M + S2) * omega2_v
             + 0.5 * arma::accu(arma::log(omega2_v)) + ki(Y);
    }

    static Rcpp::List output_cov(const arma::mat & /*M*/, const arma::mat & /*S2*/,
                                  const arma::vec & /*w*/, double /*w_bar*/, const State & s) {
        arma::uword p = s.omega2.n_elem;
        arma::sp_mat Sigma_out(p, p); Sigma_out.diag() = s.sigma2.t();
        arma::sp_mat Omega_out(p, p); Omega_out.diag() = s.omega2.t();
        return Rcpp::List::create(Rcpp::Named("Sigma", Sigma_out), Rcpp::Named("Omega", Omega_out));
    }

    static constexpr bool has_em = true;
};

// ─────────────────────────────────────────────────────────────────────────────
// Spherical covariance (scalar sigma²)
// ─────────────────────────────────────────────────────────────────────────────
struct SphericalCovTraits {
    struct State {
        double omega2;
        double sigma2;

        State(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            arma::uword p = M.n_cols;
            sigma2 = arma::accu(arma::diagmat(w) * (M % M + S2)) / (double(p) * w_bar);
            omega2 = 1.0 / sigma2;
        }
        // vestep: Omega known and fixed (scalar diagonal matrix)
        explicit State(const arma::mat & omega_mat)
            : omega2(omega_mat(0, 0)), sigma2(1.0 / omega_mat(0, 0)) {}
    };

    // returns double: fixed_point_psi<double> handles scalar broadcast
    static double cov_diag(const State & s, const arma::mat & /*ones_row*/) {
        return s.omega2;
    }

    static void grad_hess_M(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & Y, const arma::vec & w, const arma::mat & /*ones_row*/,
        arma::mat & grad_M, arma::mat & hess_M)
    {
        grad_M = s.omega2 * M + A - Y; grad_M.each_col() %= w;
        hess_M = s.omega2 + A;          hess_M.each_col() %= w;
    }

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return s.omega2 * M; }

    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    static double objective_cov(const arma::mat & M, const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * s.omega2 * arma::accu(arma::diagmat(w) * (M % M + S2));
    }

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword p) {
        s.sigma2 = arma::accu(arma::diagmat(w) * (M % M + S2)) / (double(p) * w_bar);
        s.omega2 = 1.0 / s.sigma2;
    }

    static double elbo_cov(const State & s, double w_bar, arma::uword p) {
        return -0.5 * w_bar * double(p) * std::log(s.sigma2);
    }

    static arma::vec final_loglik(const arma::mat & Y, const arma::mat & Z, const arma::mat & A,
                                   const arma::mat & M, const arma::mat & psi, const State & s) {
        const arma::mat S2 = arma::exp(psi);
        return arma::sum(Y % Z - A - 0.5 * (M % M + S2) / s.sigma2 + 0.5 * (psi - std::log(s.sigma2)), 1)
             + ki(Y);
    }

    static Rcpp::List output_cov(const arma::mat & M, const arma::mat & /*S2*/,
                                  const arma::vec & /*w*/, double /*w_bar*/, const State & s) {
        arma::uword p = M.n_cols;
        arma::sp_mat Sigma_out(p, p); Sigma_out.diag() = arma::ones<arma::vec>(p) * s.sigma2;
        arma::sp_mat Omega_out(p, p); Omega_out.diag() = arma::ones<arma::vec>(p) * s.omega2;
        return Rcpp::List::create(Rcpp::Named("Sigma", Sigma_out), Rcpp::Named("Omega", Omega_out));
    }

    static constexpr bool has_em = true;
};

// ─────────────────────────────────────────────────────────────────────────────
// Fixed covariance (Omega provided externally, not estimated)
// ─────────────────────────────────────────────────────────────────────────────
struct FixedCovTraits {
    struct State {
        arma::mat Omega;
        arma::vec diag_Omega;

        explicit State(const arma::mat & omega)
            : Omega(omega), diag_Omega(arma::diagvec(omega)) {}
    };

    static arma::mat cov_diag(const State & s, const arma::mat & ones_row) {
        return ones_row * s.diag_Omega.t();
    }

    static void grad_hess_M(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & Y, const arma::vec & w, const arma::mat & ones_row,
        arma::mat & grad_M, arma::mat & hess_M)
    {
        arma::mat MO = M * s.Omega;
        grad_M = MO + A - Y; grad_M.each_col() %= w;
        hess_M = A + ones_row * s.diag_Omega.t(); hess_M.each_col() %= w;
    }

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return M * s.Omega; }

    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    static double objective_cov(const arma::mat & M, const arma::mat & S2, const State & s, const arma::vec & w) {
        arma::mat MO = M * s.Omega;
        return 0.5 * (arma::as_scalar(w.t() * arma::sum(MO % M, 1))
                    + arma::dot(s.diag_Omega, (w.t() * S2).t()));
    }

    static void mstep(State & /*s*/, const arma::mat & /*M*/, const arma::mat & /*S2*/,
                      const arma::vec & /*w*/, double /*w_bar*/, arma::uword /*p*/) {}

    static double elbo_cov(const State & /*s*/, double /*w_bar*/, arma::uword /*p*/) {
        return 0.0;
    }

    static arma::vec final_loglik(const arma::mat & Y, const arma::mat & Z, const arma::mat & A,
                                   const arma::mat & M, const arma::mat & psi, const State & s) {
        const arma::mat S2 = arma::exp(psi);
        return arma::sum(Y % Z - A + 0.5 * psi
                       - 0.5 * ((M * s.Omega) % M + S2 * arma::diagmat(s.Omega)), 1)
             + 0.5 * std::real(arma::log_det(s.Omega)) + ki(Y);
    }

    static Rcpp::List output_cov(const arma::mat & M, const arma::mat & S2,
                                  const arma::vec & w, double w_bar, const State & s) {
        arma::mat Sigma = (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2)) / w_bar;
        return Rcpp::List::create(Rcpp::Named("Sigma", Sigma), Rcpp::Named("Omega", s.Omega));
    }

    static constexpr bool has_em = false;
};

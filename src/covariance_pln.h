#pragma once
#include <RcppArmadillo.h>
#include "utils.h"

// ─────────────────────────────────────────────────────────────────────────────
// Shared base for dense (full p×p) Omega variants (FullCovTraits and FixedCovTraits).
// Contains the static methods that are identical in both. Derived traits use
// struct inheritance to expose these methods without repetition.
// ─────────────────────────────────────────────────────────────────────────────
struct DenseOmegaImpl {
    struct State {
        arma::mat Omega;
        arma::vec diag_Omega;
    };

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return M * s.Omega; }

    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    static double penalty_S(const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * arma::dot(s.diag_Omega, (w.t() * S2).t());
    }

    // Joint Newton step for (M, ψ) where ψ = log(S²): diagonal 2×2 per (i,j).
    // MO (output) = M * Omega — returned so the caller can reuse it for penalty/Armijo
    // without an extra O(n p²) matrix product.
    static void compute_joint_step_MS(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & S2,
        const arma::mat & Y, const arma::vec & w, const arma::mat & ones_row,
        arma::mat & grad_M, arma::mat & step_M,
        arma::mat & grad_psi, arma::mat & step_psi,
        arma::mat & MO)
    {
        MO = M * s.Omega;
        const arma::mat omega_d = ones_row * s.diag_Omega.t();
        const arma::mat AS2     = A % S2;

        grad_M   = MO + A - Y;                        grad_M.each_col()   %= w;
        grad_psi = 0.5 * (AS2 + omega_d % S2 - 1.0); grad_psi.each_col() %= w;

        arma::mat h_pp = 0.5 * (S2 % (A % (1.0 + 0.5*S2) + omega_d)); h_pp.each_col() %= w;
        arma::mat h_mp = 0.5 * AS2;                                     h_mp.each_col() %= w;
        arma::mat h_mm = A + omega_d;                                    h_mm.each_col() %= w;

        arma::mat det = h_mm % h_pp - h_mp % h_mp;
        det.clamp(1e-20, arma::datum::inf);
        step_M   = (h_pp % grad_M   - h_mp % grad_psi) / det;
        step_psi = (h_mm % grad_psi - h_mp % grad_M  ) / det;
    }

    // VE-step objective + gradient: B and Omega fixed, only (M, S²) optimized.
    static double vestep_obj_grad(
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const State & s,
        const arma::mat & Y,    const arma::vec & w,
        arma::mat & gM, arma::mat & gPS)
    {
        const arma::mat A  = arma::exp(Z + 0.5 * S2);
        const arma::mat MO = M_res * s.Omega;
        gM  = MO + A - Y;                                              gM.each_col()  %= w;
        gPS = 0.5*(S2.each_row() % s.diag_Omega.t() + S2%A - 1.0);   gPS.each_col() %= w;
        return arma::accu(w.t() * (A - Y%Z - 0.5*logS2))
             + 0.5*(arma::accu(MO % (M_res.each_col() % w))
                  + arma::dot(s.diag_Omega, (w.t() * S2).t()));
    }

    static arma::vec final_loglik(const arma::mat & Y, const arma::mat & Z, const arma::mat & A,
                                   const arma::mat & M, const arma::mat & psi, const State & s) {
        const arma::mat S2 = arma::exp(psi);
        return arma::sum(Y % Z - A + 0.5 * psi
                       - 0.5 * ((M * s.Omega) % M + S2 * arma::diagmat(s.Omega)), 1)
             + 0.5 * std::real(arma::log_det(s.Omega)) + ki(Y);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Full (dense p×p) covariance
// ─────────────────────────────────────────────────────────────────────────────
struct FullCovTraits : DenseOmegaImpl {
    struct State : DenseOmegaImpl::State {
        arma::mat Sigma;
        double log_det_Sigma = 0.0;  // cached by update() — avoids a second Cholesky in elbo_cov

        State(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            update(M, S2, w, w_bar);
        }
        explicit State(const arma::mat & omega) {
            Omega      = omega;
            diag_Omega = arma::diagvec(omega);
        }
        void update(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            Sigma      = (1./w_bar) * (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2));
            Omega      = arma::inv_sympd(Sigma);
            diag_Omega = arma::diagvec(Omega);
            double sign;
            arma::log_det(log_det_Sigma, sign, Sigma);
        }
    };

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.update(M, S2, w, w_bar);
    }

    static double elbo_cov(const State & s, double w_bar, arma::uword /*p*/) {
        return -0.5 * w_bar * s.log_det_Sigma;
    }

    // Profiled joint objective + gradient: profiles Omega from (M_res, S2) at each eval.
    // Returns data_term + (w_bar/2)*log|Sigma_hat|  (envelope theorem: same gradient as
    // fixed-Omega objective).  Only called for small p since it costs O(np² + p³) per eval.
    static double profile_and_grad(
        State & s,
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const arma::mat & Y,    const arma::vec & w,
        double w_bar, arma::uword /*p*/,
        arma::mat & gM, arma::mat & gPS)
    {
        s.update(M_res, S2, w, w_bar);  // O(np²) + O(p³); caches log_det_Sigma
        const arma::mat A  = arma::exp(Z + 0.5 * S2);
        const arma::mat MO = M_res * s.Omega;
        gM  = MO + A - Y;                                              gM.each_col()  %= w;
        gPS = 0.5*(S2.each_row() % s.diag_Omega.t() + S2%A - 1.0);   gPS.each_col() %= w;
        return arma::accu(w.t() * (A - Y%Z - 0.5*logS2)) - elbo_cov(s, w_bar, 0);
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
            update(M, S2, w, w_bar);
        }
        explicit State(const arma::mat & omega_mat) {
            omega2 = arma::diagvec(omega_mat).t();
            sigma2 = arma::pow(omega2, -1);
        }
        void update(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            sigma2 = arma::clamp(w.t() * (M % M + S2) / w_bar, 1e-20, arma::datum::inf);
            omega2 = arma::pow(sigma2, -1);
        }
    };

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return M.each_row() % s.omega2; }

    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    static double penalty_S(const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * arma::as_scalar((w.t() * S2) * s.omega2.t());
    }

    // MO (output) = M.each_row() % omega2 — returned for caller reuse.
    static void compute_joint_step_MS(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & S2,
        const arma::mat & Y, const arma::vec & w, const arma::mat & ones_row,
        arma::mat & grad_M, arma::mat & step_M,
        arma::mat & grad_psi, arma::mat & step_psi,
        arma::mat & MO)
    {
        const arma::mat omega_d = ones_row * s.omega2;
        const arma::mat AS2     = A % S2;
        MO       = M.each_row() % s.omega2;
        grad_M   = MO + A - Y;                            grad_M.each_col()   %= w;
        grad_psi = 0.5 * (AS2 + omega_d % S2 - 1.0);     grad_psi.each_col() %= w;

        arma::mat h_pp = 0.5 * (S2 % (A % (1.0 + 0.5*S2) + omega_d)); h_pp.each_col() %= w;
        arma::mat h_mp = 0.5 * AS2;                                     h_mp.each_col() %= w;
        arma::mat h_mm = A + omega_d;                                    h_mm.each_col() %= w;

        arma::mat det = h_mm % h_pp - h_mp % h_mp;
        det.clamp(1e-20, arma::datum::inf);
        step_M   = (h_pp % grad_M   - h_mp % grad_psi) / det;
        step_psi = (h_mm % grad_psi - h_mp % grad_M  ) / det;
    }

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.update(M, S2, w, w_bar);
    }

    static double elbo_cov(const State & s, double w_bar, arma::uword /*p*/) {
        return -0.5 * w_bar * arma::accu(arma::log(s.sigma2));
    }

    // VE-step objective + gradient: B and Omega fixed, only (M, S²) optimized.
    static double vestep_obj_grad(
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const State & s,
        const arma::mat & Y,    const arma::vec & w,
        arma::mat & gM, arma::mat & gPS)
    {
        const arma::mat A = arma::exp(Z + 0.5 * S2);
        gM  = M_res.each_row() % s.omega2 + A - Y;              gM.each_col()  %= w;
        gPS = 0.5*(S2.each_row() % s.omega2 + S2 % A - 1.0);   gPS.each_col() %= w;
        return arma::accu(w.t() * (A - Y%Z - 0.5*logS2))
             + 0.5 * arma::as_scalar(w.t() * (arma::pow(M_res, 2) + S2) * s.omega2.t());
    }

    // Profiled joint objective + gradient: profiles sigma2 from (M_res, S2) at each eval.
    static double profile_and_grad(
        State & s,
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const arma::mat & Y,    const arma::vec & w,
        double w_bar, arma::uword p,
        arma::mat & gM, arma::mat & gPS)
    {
        s.update(M_res, S2, w, w_bar);
        const arma::mat A = arma::exp(Z + 0.5 * S2);
        gM  = M_res.each_row() % s.omega2 + A - Y;              gM.each_col()  %= w;
        gPS = 0.5*(S2.each_row() % s.omega2 + S2 % A - 1.0);   gPS.each_col() %= w;
        return arma::accu(w.t() * (A - Y%Z - 0.5*logS2)) - elbo_cov(s, w_bar, p);
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
            update(M, S2, w, w_bar);
        }
        explicit State(const arma::mat & omega_mat)
            : omega2(omega_mat(0, 0)), sigma2(1.0 / omega_mat(0, 0)) {}
        void update(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            arma::uword p = M.n_cols;
            sigma2 = arma::accu(w.t() * (M % M + S2)) / (double(p) * w_bar);
            omega2 = 1.0 / sigma2;
        }
    };

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return s.omega2 * M; }

    static double penalty_S(const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * s.omega2 * arma::dot(w, arma::sum(S2, 1));
    }

    // MO (output) = omega2 * M — returned for caller reuse.
    static void compute_joint_step_MS(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & S2,
        const arma::mat & Y, const arma::vec & w, const arma::mat & /*ones_row*/,
        arma::mat & grad_M, arma::mat & step_M,
        arma::mat & grad_psi, arma::mat & step_psi,
        arma::mat & MO)
    {
        const arma::mat AS2 = A % S2;
        MO       = s.omega2 * M;
        grad_M   = MO + A - Y;                          grad_M.each_col()   %= w;
        grad_psi = 0.5 * (AS2 + s.omega2 * S2 - 1.0);  grad_psi.each_col() %= w;

        arma::mat h_pp = 0.5 * (S2 % (A % (1.0 + 0.5*S2) + s.omega2)); h_pp.each_col() %= w;
        arma::mat h_mp = 0.5 * AS2;                                      h_mp.each_col() %= w;
        arma::mat h_mm = A + s.omega2;                                    h_mm.each_col() %= w;

        arma::mat det = h_mm % h_pp - h_mp % h_mp;
        det.clamp(1e-20, arma::datum::inf);
        step_M   = (h_pp % grad_M   - h_mp % grad_psi) / det;
        step_psi = (h_mm % grad_psi - h_mp % grad_M  ) / det;
    }

    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.update(M, S2, w, w_bar);
    }

    static double elbo_cov(const State & s, double w_bar, arma::uword p) {
        return -0.5 * w_bar * double(p) * std::log(s.sigma2);
    }

    // VE-step objective + gradient: B and Omega fixed, only (M, S²) optimized.
    static double vestep_obj_grad(
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const State & s,
        const arma::mat & Y,    const arma::vec & w,
        arma::mat & gM, arma::mat & gPS)
    {
        const arma::mat A = arma::exp(Z + 0.5 * S2);
        gM  = s.omega2 * M_res + A - Y;               gM.each_col()  %= w;
        gPS = 0.5*(s.omega2 * S2 + S2 % A - 1.0);    gPS.each_col() %= w;
        return arma::accu(w.t() * (A - Y%Z - 0.5*logS2))
             + 0.5 * s.omega2 * arma::accu(w.t() * (arma::pow(M_res, 2) + S2));
    }

    // Profiled joint objective + gradient: profiles sigma2 from (M_res, S2) at each eval.
    static double profile_and_grad(
        State & s,
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const arma::mat & Y,    const arma::vec & w,
        double w_bar, arma::uword p,
        arma::mat & gM, arma::mat & gPS)
    {
        s.update(M_res, S2, w, w_bar);
        const arma::mat A = arma::exp(Z + 0.5 * S2);
        gM  = s.omega2 * M_res + A - Y;               gM.each_col()  %= w;
        gPS = 0.5*(s.omega2 * S2 + S2 % A - 1.0);    gPS.each_col() %= w;
        return arma::accu(w.t() * (A - Y%Z - 0.5*logS2)) - elbo_cov(s, w_bar, p);
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
        arma::sp_mat Sigma_out(p, p); Sigma_out.diag().fill(s.sigma2);
        arma::sp_mat Omega_out(p, p); Omega_out.diag().fill(s.omega2);
        return Rcpp::List::create(Rcpp::Named("Sigma", Sigma_out), Rcpp::Named("Omega", Omega_out));
    }

    static constexpr bool has_em = true;
};

// ─────────────────────────────────────────────────────────────────────────────
// Fixed covariance (Omega provided externally, not estimated)
// ─────────────────────────────────────────────────────────────────────────────
struct FixedCovTraits : DenseOmegaImpl {
    struct State : DenseOmegaImpl::State {
        explicit State(const arma::mat & omega) {
            Omega      = omega;
            diag_Omega = arma::diagvec(omega);
        }
    };

    static void mstep(State & /*s*/, const arma::mat & /*M*/, const arma::mat & /*S2*/,
                      const arma::vec & /*w*/, double /*w_bar*/, arma::uword /*p*/) {}

    static double elbo_cov(const State & /*s*/, double /*w_bar*/, arma::uword /*p*/) {
        return 0.0;
    }

    static Rcpp::List output_cov(const arma::mat & M, const arma::mat & S2,
                                  const arma::vec & w, double w_bar, const State & s) {
        arma::mat Sigma = (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2)) / w_bar;
        return Rcpp::List::create(Rcpp::Named("Sigma", Sigma), Rcpp::Named("Omega", s.Omega));
    }

    static constexpr bool has_em = false;
};

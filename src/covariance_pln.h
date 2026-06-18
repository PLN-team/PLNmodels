#pragma once
#include <RcppArmadillo.h>
#include "utils.h"

// ─────────────────────────────────────────────────────────────────────────────
// CRTP base shared by all covariance traits (Full, Fixed, Diagonal, Spherical).
// Implements the optimization machinery that is structurally identical across
// covariance structures, in terms of a handful of primitives each concrete
// trait (Derived) supplies on its own State:
//   - times_Omega(M, s)     = M * Omega                              (n×p)
//   - diag_scale(S2, s)     = S2 ⊙ diag(Omega), broadcast over rows  (n×p)
//   - add_diag(X, s)        = X + diag(Omega), broadcast over rows   (n×p)
//   - penalty_S(S2, s, w)   = ½ Σ w ⊙ S2 ⊙ diag(Omega)               (scalar)
//   - elbo_cov(s, w_bar, p) = -½ w_bar log|Sigma|                    (scalar)
//
// Each method below is itself templated on State (deduced from the call site)
// rather than using `typename Derived::State`: Derived is still incomplete at
// the point FullCovTraits/etc. inherit from CovTraitsBase<Derived>, so State
// can only be referenced inside (lazily-instantiated) method bodies, never in
// a declaration that the compiler must resolve immediately.
// ─────────────────────────────────────────────────────────────────────────────
template <typename Derived>
struct CovTraitsBase {
    static double penalty_M(const arma::mat & MO, const arma::mat & M, const arma::vec & w) {
        return 0.5 * arma::as_scalar(w.t() * arma::sum(MO % M, 1));
    }

    // Objective value only (no gradient) at a given point — used by the builtin Newton
    // solver's Armijo line search, which evaluates several trial points per step without
    // needing their gradient. A, MO are passed in (already available/cheaply updated at
    // the call site) rather than recomputed here.
    template <typename State>
    static double objective(
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const arma::mat & MO,
        const State & s,
        const arma::mat & Y,    const arma::vec & w,
        const arma::mat & A)
    {
        return arma::accu(w.t() * (A - Y % Z - 0.5 * logS2)) + penalty_M(MO, M_res, w) + Derived::penalty_S(S2, s, w);
    }

    // VE-step objective + gradient: B and Omega fixed, only (M, S²) optimized.
    template <typename State>
    static double vestep_obj_grad(
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const State & s,
        const arma::mat & Y,    const arma::vec & w,
        arma::mat & gM, arma::mat & gPS)
    {
        arma::mat MO;
        const double data_term = vestep_core(M_res, Z, S2, logS2, s, Y, w, gM, gPS, MO);
        return data_term + penalty_M(MO, M_res, w) + Derived::penalty_S(S2, s, w);
    }

    // Profiled joint objective + gradient: profiles Omega/sigma² from (M_res, S2) at
    // every eval (envelope theorem: same gradient as the fixed-Omega objective above).
    // Only the returned value differs: the quadratic penalty collapses to elbo_cov
    // because Omega is the MLE of (M_res, S2) at every evaluation.
    template <typename State>
    static double profile_and_grad(
        State & s,
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const arma::mat & Y,    const arma::vec & w,
        double w_bar, arma::uword p,
        arma::mat & gM, arma::mat & gPS)
    {
        s.update(M_res, S2, w, w_bar);
        arma::mat MO;
        const double data_term = vestep_core(M_res, Z, S2, logS2, s, Y, w, gM, gPS, MO);
        return data_term - Derived::elbo_cov(s, w_bar, p);
    }

    // Joint Newton step for (M, ψ) where ψ = log(S²): diagonal 2×2 per (i,j).
    // MO (output) = M * Omega — returned so the caller can reuse it for penalty/Armijo
    // without an extra O(n p²) matrix product.
    template <typename State>
    static void compute_joint_step_MS(
        const arma::mat & M, const State & s,
        const arma::mat & A, const arma::mat & S2,
        const arma::mat & Y, const arma::vec & w,
        arma::mat & grad_M, arma::mat & step_M,
        arma::mat & grad_psi, arma::mat & step_psi,
        arma::mat & MO)
    {
        MO = Derived::times_Omega(M, s);
        const arma::mat AS2      = A % S2;
        const arma::mat diag_AS2 = Derived::diag_scale(S2, s);

        grad_M   = MO + A - Y;                   grad_M.each_col()   %= w;
        grad_psi = 0.5 * (AS2 + diag_AS2 - 1.0); grad_psi.each_col() %= w;

        arma::mat h_pp = 0.5 * (S2 % Derived::add_diag(A % (1.0 + 0.5*S2), s)); h_pp.each_col() %= w;
        arma::mat h_mp = 0.5 * AS2;                                             h_mp.each_col() %= w;
        arma::mat h_mm = Derived::add_diag(A, s);                               h_mm.each_col() %= w;

        arma::mat det = h_mm % h_pp - h_mp % h_mp;
        det.clamp(1e-20, arma::datum::inf);
        step_M   = (h_pp % grad_M   - h_mp % grad_psi) / det;
        step_psi = (h_mm % grad_psi - h_mp % grad_M  ) / det;
    }

private:
    // Shared core: A, MO, gradients and the data-fit term (no covariance penalty —
    // callers above add either the explicit penalty_M/penalty_S or -elbo_cov).
    template <typename State>
    static double vestep_core(
        const arma::mat & M_res, const arma::mat & Z,
        const arma::mat & S2,   const arma::mat & logS2,
        const State & s,
        const arma::mat & Y,    const arma::vec & w,
        arma::mat & gM, arma::mat & gPS, arma::mat & MO)
    {
        const arma::mat A = arma::exp(Z + 0.5 * S2);
        MO  = Derived::times_Omega(M_res, s);
        gM  = MO + A - Y;                                            gM.each_col()  %= w;
        gPS = 0.5 * (Derived::diag_scale(S2, s) + S2 % A - 1.0);   gPS.each_col() %= w;
        return arma::accu(w.t() * (A - Y % Z - 0.5 * logS2));
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Shared base for dense (full p×p) Omega variants (FullCovTraits and FixedCovTraits).
// Provides the primitives required by CovTraitsBase, plus final_loglik (which
// differs too much in form across structures to be worth factoring further).
// ─────────────────────────────────────────────────────────────────────────────
struct DenseOmegaImpl {
    struct State {
        arma::mat Omega;
        arma::vec diag_Omega;
    };

    static arma::mat times_Omega(const arma::mat & M, const State & s) { return M * s.Omega; }

    static arma::mat diag_scale(const arma::mat & S2, const State & s) {
        return S2.each_row() % s.diag_Omega.t();
    }

    static arma::mat add_diag(arma::mat X, const State & s) {
        X.each_row() += s.diag_Omega.t();
        return X;
    }

    static double penalty_S(const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * arma::dot(s.diag_Omega, (w.t() * S2).t());
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
struct FullCovTraits : DenseOmegaImpl, CovTraitsBase<FullCovTraits> {
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

    static Rcpp::List output_cov(const arma::mat & /*M*/, const arma::mat & /*S2*/,
                                  const arma::vec & /*w*/, double /*w_bar*/, const State & s) {
        return Rcpp::List::create(Rcpp::Named("Sigma", s.Sigma), Rcpp::Named("Omega", s.Omega));
    }

    static constexpr bool has_em = true;
};

// ─────────────────────────────────────────────────────────────────────────────
// Diagonal covariance
// ─────────────────────────────────────────────────────────────────────────────
struct DiagonalCovTraits : CovTraitsBase<DiagonalCovTraits> {
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

    static arma::mat diag_scale(const arma::mat & S2, const State & s) { return S2.each_row() % s.omega2; }

    static arma::mat add_diag(arma::mat X, const State & s) {
        X.each_row() += s.omega2;
        return X;
    }

    static double penalty_S(const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * arma::as_scalar((w.t() * S2) * s.omega2.t());
    }

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.update(M, S2, w, w_bar);
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
struct SphericalCovTraits : CovTraitsBase<SphericalCovTraits> {
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

    static arma::mat diag_scale(const arma::mat & S2, const State & s) { return s.omega2 * S2; }

    static arma::mat add_diag(arma::mat X, const State & s) { return X + s.omega2; }

    static double penalty_S(const arma::mat & S2, const State & s, const arma::vec & w) {
        return 0.5 * s.omega2 * arma::dot(w, arma::sum(S2, 1));
    }

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.update(M, S2, w, w_bar);
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
        arma::sp_mat Sigma_out(p, p); Sigma_out.diag().fill(s.sigma2);
        arma::sp_mat Omega_out(p, p); Omega_out.diag().fill(s.omega2);
        return Rcpp::List::create(Rcpp::Named("Sigma", Sigma_out), Rcpp::Named("Omega", Omega_out));
    }

    static constexpr bool has_em = true;
};

// ─────────────────────────────────────────────────────────────────────────────
// Fixed covariance (Omega provided externally, not estimated)
// ─────────────────────────────────────────────────────────────────────────────
struct FixedCovTraits : DenseOmegaImpl, CovTraitsBase<FixedCovTraits> {
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

// ─────────────────────────────────────────────────────────────────────────────
// Genetic / heritability covariance: Sigma = sigma2 * (rho * C + (1-rho) * I_p),
// C a fixed p×p correlation matrix supplied by the user (e.g. a genetic
// relationship matrix). With C = V * diag(Lambda) * V' (eigendecomposition,
// computed once and cached in State — never recomputed by update()):
//   Sigma  = sigma2 * V * diag(u) * V',  u_j = rho * Lambda_j + (1 - rho)
//   Omega  =          V * diag(1 / (sigma2 * u)) * V'
// sigma2 has a closed form given rho (envelope theorem, same trick as the
// other profiled traits); rho itself has no closed form and is found by 1-D
// bisection on its profiled gradient — see VEM-PLN-struct_var.tex, section
// "A covariance modeling heritability", for the derivation (math verified
// term-by-term against that document and against the pre-existing, separate
// nlopt_optimize_genetic_modeling in optim_genet_cov.cpp before this port).
// Because Omega ends up dense, Genetic reuses DenseOmegaImpl's primitives
// (times_Omega, diag_scale, add_diag, penalty_S, final_loglik) unchanged —
// the only genuinely new code is State::update() below.
// ─────────────────────────────────────────────────────────────────────────────
struct GeneticCovTraits : DenseOmegaImpl, CovTraitsBase<GeneticCovTraits> {
    struct State : DenseOmegaImpl::State {
        arma::mat V;
        arma::vec Lambda;
        double sigma2 = 1.0;
        double rho    = 0.5;

        State() = default;
        // VE-step use: Omega given externally (already the result of a previous
        // M-step) — V/Lambda/sigma2/rho are not needed by the VE-step itself.
        explicit State(const arma::mat & omega) {
            Omega      = omega;
            diag_Omega = arma::diagvec(omega);
        }
        // One-time setup from the fixed correlation matrix C: eigendecomposition
        // only (cost O(p^3), same order as Full's inv_sympd, paid once per EM
        // iteration — not on every nlopt eval). Call update() right after to get
        // a usable (rho, sigma2, Omega) fit.
        static State from_correlation(const arma::mat & C) {
            State s;
            arma::eig_sym(s.Lambda, s.V, C);
            s.Lambda = arma::clamp(s.Lambda, 0.0, arma::datum::inf);  // guard tiny negative eigenvalues
            return s;
        }
        void update(const arma::mat & M, const arma::mat & S2, const arma::vec & w, double w_bar) {
            const arma::uword p = M.n_cols;
            const arma::mat R = V.t() * (M.t() * (M.each_col() % w) + arma::diagmat(w.t() * S2)) * V;
            const arma::vec Rdiag = arma::diagvec(R);

            rho = solve_rho(Rdiag, w_bar);
            const arma::vec u = arma::clamp(rho * Lambda + (1.0 - rho), 1e-12, arma::datum::inf);
            sigma2 = arma::accu(Rdiag / u) / (double(p) * w_bar);

            Omega      = V * arma::diagmat(1.0 / (sigma2 * u)) * V.t();
            diag_Omega = arma::diagvec(Omega);
        }

    private:
        // Profiled dJ/drho (sigma2 substituted at its closed-form optimum via the
        // envelope theorem) — see proposition "Gradient of J" (heritability) in
        // VEM-PLN-struct_var.tex.
        double dJ_drho(double r, const arma::vec & Rdiag, double w_bar) const {
            const arma::vec u      = arma::clamp(r * Lambda + (1.0 - r), 1e-12, arma::datum::inf);
            const arma::vec lam_m1 = Lambda - 1.0;
            const double s2        = arma::accu(Rdiag / u) / (double(Rdiag.n_elem) * w_bar);
            return -0.5 * w_bar * arma::accu(lam_m1 / u)
                  + 0.5 / s2 * arma::accu(Rdiag % lam_m1 / arma::square(u));
        }
        // J(rho) is a single variance-ratio profile, unimodal in practice over
        // (0,1); bisection on the sign of dJ/drho is robust and cheap (O(p) per
        // eval) — falls back to the nearest boundary if the gradient never
        // changes sign (optimum at rho = 0 or rho = 1).
        double solve_rho(const arma::vec & Rdiag, double w_bar) const {
            constexpr double eps = 1e-8;
            double lo = eps, hi = 1.0 - eps;
            const double g_lo = dJ_drho(lo, Rdiag, w_bar);
            const double g_hi = dJ_drho(hi, Rdiag, w_bar);
            if (g_lo <= 0.0) return lo;
            if (g_hi >= 0.0) return hi;
            for (int it = 0; it < 60; it++) {
                const double mid = 0.5 * (lo + hi);
                if (dJ_drho(mid, Rdiag, w_bar) > 0.0) lo = mid; else hi = mid;
            }
            return 0.5 * (lo + hi);
        }
    };

    static void mstep(State & s, const arma::mat & M, const arma::mat & S2,
                      const arma::vec & w, double w_bar, arma::uword /*p*/) {
        s.update(M, S2, w, w_bar);
    }

    static double elbo_cov(const State & s, double w_bar, arma::uword p) {
        const arma::vec u = arma::clamp(s.rho * s.Lambda + (1.0 - s.rho), 1e-12, arma::datum::inf);
        return -0.5 * w_bar * (double(p) * std::log(s.sigma2) + arma::accu(arma::log(u)));
    }

    static Rcpp::List output_cov(const arma::mat & /*M*/, const arma::mat & /*S2*/,
                                  const arma::vec & /*w*/, double /*w_bar*/, const State & s) {
        const arma::vec u = arma::clamp(s.rho * s.Lambda + (1.0 - s.rho), 1e-12, arma::datum::inf);
        const arma::mat Sigma = s.V * arma::diagmat(s.sigma2 * u) * s.V.t();
        return Rcpp::List::create(
            Rcpp::Named("Sigma",  Sigma),    Rcpp::Named("Omega", s.Omega),
            Rcpp::Named("sigma2", s.sigma2), Rcpp::Named("rho",   s.rho));
    }

    static constexpr bool has_em = true;
};

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"
#include "covariance_plnpca.h"
#include "trnewton_plnpca.h"

using arma::mat; using arma::vec; using arma::cube;

// A point in θ-space = (B, C). Small helpers for the CG on θ.
namespace {
struct Th { mat B, C; };
inline double dot(const Th & a, const Th & b) { return arma::accu(a.B % b.B) + arma::accu(a.C % b.C); }
inline double nrm(const Th & a) { return std::sqrt(dot(a, a)); }
inline Th axpy(double s, const Th & x, const Th & y) { return { s * x.B + y.B, s * x.C + y.C }; } // s*x + y
inline Th scal(double s, const Th & x) { return { s * x.B, s * x.C }; }
// preconditioner helpers (P = diag, stored as a Th of positive diagonal entries)
inline double dotP(const Th & a, const Th & b, const Th & p) {
    return arma::accu(p.B % a.B % b.B) + arma::accu(p.C % a.C % b.C); }
inline Th Pinv(const Th & r, const Th & p) { return { r.B / p.B, r.C / p.C }; }
} // namespace

// [[Rcpp::export]]
Rcpp::List trnewton_optimize_rank(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params, // List(B, C, M, S2)
    const Rcpp::List & config  // List(maxit_out, ftol_out, cg_maxit, delta0, trace)
) {
    const PlnData d(data);
    mat B   = Rcpp::as<mat>(params["B"]);
    mat C   = Rcpp::as<mat>(params["C"]);
    mat M   = Rcpp::as<mat>(params["M"]);
    mat psi = arma::log(Rcpp::as<mat>(params["S2"]));

    const int    maxit_out = config.containsElementNamed("maxit_out") ? Rcpp::as<int>(config["maxit_out"]) : 50;
    const double ftol_out  = config.containsElementNamed("ftol_out")  ? Rcpp::as<double>(config["ftol_out"]) : 1e-8;
    const int    cg_maxit  = config.containsElementNamed("cg_maxit")  ? Rcpp::as<int>(config["cg_maxit"]) : 25;
    const double gtol      = config.containsElementNamed("gtol")      ? Rcpp::as<double>(config["gtol"]) : 1e-4;
    double       Delta     = config.containsElementNamed("delta0")    ? Rcpp::as<double>(config["delta0"]) : 1.0;
    const int    trace     = config.containsElementNamed("trace")     ? Rcpp::as<int>(config["trace"]) : 0;

    const arma::uword n = d.Y.n_rows, q = C.n_cols;
    const mat Xw = d.X.each_col() % d.w;

    std::vector<double> objective_vec;
    int status = 5;

    mat XB = d.X * B, C2 = C % C;   // kept in sync with B, C across accepted steps
    const mat XX = d.X % d.X;       // for the Jacobi preconditioner (diag L_BB)

    // base gradients at (B,C,M,psi); envelope grad = (gB0,gC0)
    mat gB0, gC0, gM0, gPS0;
    const mat zB = arma::zeros(B.n_rows, B.n_cols), zC = arma::zeros(C.n_rows, C.n_cols);
    const mat zMq = arma::zeros(n, q);

    // Schur reduced Hessian-vector product H_red·v = L_θθ v − L_θφ L_φφ⁻¹ L_φθ v.
    // All blocks computed ANALYTICALLY via `hess_dir`: one call with the θ-direction v
    // yields L_θθ v (θ part) and L_φθ v (φ part); one call with the φ-direction w=L_φφ⁻¹(…)
    // yields L_θφ w (θ part). A, S2 are current values (recomputed once per outer iter).
    auto hessvec = [&](const Th & v, const cube & Hinv, const mat & A, const mat & S2) -> Th {
        mat Ltt_B, Ltt_C, Lft_M, Lft_P;
        trnewton::hess_dir(d, Xw, C, C2, M, A, S2, v.B, v.C, zMq, zMq, Ltt_B, Ltt_C, Lft_M, Lft_P);
        mat wM(n, q), wP(n, q);
        for (arma::uword i = 0; i < n; ++i) {          // w = L_φφ⁻¹ (L_φθ v), per observation
            vec x = Hinv.slice(i) * arma::join_cols(Lft_M.row(i).t(), Lft_P.row(i).t());
            wM.row(i) = x.head(q).t();
            wP.row(i) = x.tail(q).t();
        }
        mat Ltf_B, Ltf_C, dummyM, dummyP;
        trnewton::hess_dir(d, Xw, C, C2, M, A, S2, zB, zC, wM, wP, Ltf_B, Ltf_C, dummyM, dummyP);
        return { Ltt_B - Ltf_B, Ltt_C - Ltf_C };
    };

    // Preconditioned Steihaug-CG: approx-minimise m(s)=g·s+½ sᵀH s over the P-norm
    // ball ‖s‖_P ≤ Δ (P = Jacobi preconditioner). Handles negative curvature; the
    // P-metric absorbs the B-vs-C scale disparity, cutting the CG iteration count.
    auto steihaug = [&](const Th & g, const cube & Hinv, const mat & A, const mat & S2, const Th & p) -> Th {
        Th z{ arma::zeros(B.n_rows, B.n_cols), arma::zeros(C.n_rows, C.n_cols) };
        Th r = g;                       // ∇m at z=0
        Th y = Pinv(r, p);
        Th dvec = scal(-1.0, y);        // preconditioned steepest descent
        double ry = dot(r, y);
        const double tol = 1e-8 * std::sqrt(std::abs(ry));
        auto to_boundary = [&](const Th & zz, const Th & dd) -> Th {  // ‖zz+τ dd‖_P = Δ, τ>0
            double a = dotP(dd, dd, p), b = 2 * dotP(zz, dd, p), c = dotP(zz, zz, p) - Delta * Delta;
            double tau = (-b + std::sqrt(std::max(b * b - 4 * a * c, 0.0))) / (2 * a);
            return axpy(tau, dd, zz);
        };
        for (int j = 0; j < cg_maxit; ++j) {
            if (std::sqrt(std::abs(ry)) <= tol) break;
            Th Hd = hessvec(dvec, Hinv, A, S2);
            double dHd = dot(dvec, Hd);
            if (dHd <= 0) return to_boundary(z, dvec);                 // negative curvature
            double alpha = ry / dHd;
            Th z_new = axpy(alpha, dvec, z);
            if (dotP(z_new, z_new, p) >= Delta * Delta) return to_boundary(z, dvec);
            z = z_new;
            Th r_new = axpy(alpha, Hd, r);
            Th y_new = Pinv(r_new, p);
            double ry_new = dot(r_new, y_new);
            dvec = axpy(ry_new / ry, dvec, scal(-1.0, y_new));        // d = -y_new + β d
            r = r_new; y = y_new; ry = ry_new;
        }
        return z;
    };

    double f = trnewton::ve_solve(d, XB, C, C2, M, psi);   // profile out (M,ψ) at start
    objective_vec.push_back(f);

    int it = 0;
    double g0norm = 0.0;   // ‖∇g‖ at the start, for a scale-invariant stopping rule
    for (; it < maxit_out; ++it) {
        rank_obj_grad(d, Xw, B, C, M, psi, gB0, gC0, gM0, gPS0);  // envelope grad = gB0,gC0
        Th g{ gB0, gC0 };
        double gnorm = nrm(g);
        if (it == 0) g0norm = gnorm;
        // relative gradient tolerance: absolute ‖g‖ scales with n and the counts, so an
        // absolute threshold never triggers on large datasets.
        if (gnorm <= gtol * std::max(g0norm, 1.0)) { status = 3; break; }

        const mat S2cur = arma::exp(psi);
        const mat Acur  = arma::exp(arma::clamp(d.O + XB + M * C.t() + 0.5 * S2cur * C2.t(),
                                                -arma::datum::inf, 30.0));
        cube Hinv = trnewton::inner_blocks_inv(d, XB, C, C2, M, psi);
        Th pcond;
        trnewton::precond_diag(d, XX, C, C2, M, Acur, S2cur, pcond.B, pcond.C);
        // radius lives in the P-norm now; seed it from the preconditioned gradient scale
        if (it == 0) Delta = std::max(1e-6, std::sqrt(dot(g, Pinv(g, pcond))));
        Th s = steihaug(g, Hinv, Acur, S2cur, pcond);

        // predicted reduction  -m(s) = -(g·s + ½ sᵀH s)
        Th Hs = hessvec(s, Hinv, Acur, S2cur);
        double pred = -(dot(g, s) + 0.5 * dot(s, Hs));

        // actual reduction: profile at trial (B+sB, C+sC)
        mat Bt = B + s.B, Ct = C + s.C, C2t = Ct % Ct, XBt = d.X * Bt;
        mat Mt = M, psit = psi;
        double f_new = trnewton::ve_solve(d, XBt, Ct, C2t, Mt, psit);
        double ared = f - f_new;
        double rho = (std::abs(pred) > 1e-14 && std::isfinite(f_new)) ? ared / pred : -1.0;

        if (rho > 0.1 && std::isfinite(f_new)) {       // accept
            B = Bt; C = Ct; C2 = C2t; XB = XBt; M = Mt; psi = psit;
            double df = f - f_new; f = f_new;
            objective_vec.push_back(f);
            if (std::abs(df) <= ftol_out * std::abs(f)) { status = 3; break; }
        }
        // trust-region radius update (radius is in the P-norm)
        double snorm = std::sqrt(dotP(s, s, pcond));
        if (rho < 0.25)                    Delta = 0.25 * Delta;
        else if (rho > 0.75 && snorm > 0.9 * Delta) Delta = std::min(2.0 * Delta, 1e6);
        if (Delta < 1e-10) { status = 4; break; }
        if (trace > 1) Rcpp::Rcout << "  TR it " << it << " f=" << f << " |g|=" << gnorm
                                   << " rho=" << rho << " Delta=" << Delta << "\n";
    }

    // outputs
    mat S2 = arma::exp(psi);
    mat Z  = d.O + d.X * B + M * C.t();
    mat A  = arma::exp(Z + 0.5 * S2 * (C % C).t());
    const double w_bar = arma::accu(d.w);
    Rcpp::List cov_out = rank_output_cov(M, C, S2, d.w, w_bar);
    vec loglik = rank_final_loglik(d.Y, Z, A, M, S2, psi);

    return make_plnpca_result(B, C, M, S2, Z, A, cov_out, loglik, status, "trnewton", objective_vec, it);
}

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

// ---------------------------------------------------------------------------------------
// Rank-constrained PLN — Joint L-BFGS with strong Wolfe line search
//
// All parameters [vec(B); vec(C); vec(M); vec(ψ)] are optimised simultaneously.
// Strong Wolfe line search guarantees s^T y > 0 at every accepted step, so the
// L-BFGS history always accumulates valid curvature pairs including the bilinear
// M·Cᵀ cross-curvature that block-coordinate methods miss.

// ---------------------------------------------------------------------------------------
// L-BFGS two-loop recursion: returns search direction p = -H_k · g
static arma::vec lbfgs_direction(
    const arma::vec & g,
    const std::deque<arma::vec> & sv,
    const std::deque<arma::vec> & yv
) {
    const int m = (int)sv.size();
    arma::vec q = g, alpha(m, arma::fill::zeros);
    for (int i = m-1; i >= 0; i--) {
        double rho = 1.0 / arma::dot(yv[i], sv[i]);
        alpha(i) = rho * arma::dot(sv[i], q);
        q -= alpha(i) * yv[i];
    }
    arma::vec r = q;
    if (m > 0) {
        double sy = arma::dot(sv.back(), yv.back());
        double yy = arma::dot(yv.back(), yv.back());
        if (sy > 0 && yy > 1e-20) r *= (sy / yy);
    }
    for (int i = 0; i < m; i++) {
        double rho  = 1.0 / arma::dot(yv[i], sv[i]);
        double beta = rho * arma::dot(yv[i], r);
        r += sv[i] * (alpha(i) - beta);
    }
    return -r;
}

// ---------------------------------------------------------------------------------------
// Strong Wolfe line search (Nocedal & Wright, Algorithm 3.5/3.6).
// Guarantees s^T y > 0 when slope0 < 0 and a descent direction is given.

struct WolfeStep { double scale; double f; arma::vec g; };

template<typename FG>
static WolfeStep wolfe_ls(
    const arma::vec & x0, const arma::vec & d,
    double f0, double slope0, FG fg,
    const double c1 = 1e-4, const double c2 = 0.9
) {
    auto zoom = [&](double alo, double ahi, double flo) -> WolfeStep {
        for (int j = 0; j < 20; j++) {
            double a = 0.5 * (alo + ahi);
            auto [fa, ga] = fg(x0 + a * d);
            if (fa > f0 + c1*a*slope0 || fa >= flo) { ahi = a; }
            else {
                double da = arma::dot(ga, d);
                if (std::abs(da) <= -c2 * slope0) return {a, fa, ga};
                if (da * (ahi - alo) >= 0) ahi = alo;
                alo = a;  flo = fa;
            }
        }
        double a = 0.5 * (alo + ahi);
        auto [fa, ga] = fg(x0 + a * d);
        return {a, fa, ga};
    };
    double a = 1.0, ap = 0, fp = f0;
    for (int i = 0; i < 20; i++) {
        auto [fa, ga] = fg(x0 + a * d);
        if (fa > f0 + c1*a*slope0 || (i > 0 && fa >= fp)) return zoom(ap, a, fp);
        double da = arma::dot(ga, d);
        if (std::abs(da) <= -c2 * slope0) return {a, fa, ga};
        if (da >= 0) return zoom(a, ap, fa);
        ap = a;  fp = fa;
        a  = std::min(2.0 * a, 1e6);
    }
    auto [fa, ga] = fg(x0 + a * d);
    return {a, fa, ga};
}

// ---------------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List builtin_optimize_rank(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat B  = Rcpp::as<arma::mat>(params["B"]);
    arma::mat C  = Rcpp::as<arma::mat>(params["C"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);

    const int    maxiter = config.containsElementNamed("maxeval") ? Rcpp::as<int>   (config["maxeval"])  : 10000;
    const double ftol    = config.containsElementNamed("ftol_in") ? Rcpp::as<double>(config["ftol_in"]) : 1e-9;
    const int    m_hist  = 10;

    const arma::uword n = Y.n_rows, p = Y.n_cols, q = M.n_cols, d = B.n_rows;

    // Packed-parameter offsets: x = [vec(B); vec(C); vec(M); vec(ψ)]
    const arma::uword oB = 0, oC = d*p, oM = d*p + p*q, oPsi = d*p + p*q + n*q;
    const arma::uword N  = d*p + p*q + 2*n*q;

    // Warm-start ψ with one fixed-point step
    arma::mat C2  = C % C;
    arma::mat psi = arma::log(S2);
    arma::mat Z   = O + X * B + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    psi = arma::clamp(-arma::log(1. + A * C2), -40., 40.);
    S2  = arma::exp(psi);
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    const arma::mat Xw = X.each_col() % w;

    // Joint fg evaluator for all parameters
    auto fg = [&](const arma::vec & x) -> std::pair<double, arma::vec> {
        arma::mat B_   = arma::reshape(x.subvec(oB,   oC-1   ), d, p);
        arma::mat C_   = arma::reshape(x.subvec(oC,   oM-1   ), p, q);
        arma::mat M_   = arma::reshape(x.subvec(oM,   oPsi-1 ), n, q);
        arma::mat psi_ = arma::reshape(x.subvec(oPsi, N-1    ), n, q);
        arma::mat S2_  = arma::exp(psi_);
        arma::mat C2_  = C_ % C_;
        arma::mat Z_   = O + X * B_ + M_ * C_.t();
        arma::mat A_   = arma::exp(Z_ + 0.5 * S2_ * C2_.t());
        double f = arma::accu(w.t() * (A_ - Y % Z_))
                 + 0.5 * arma::accu(w.t() * (M_ % M_ + S2_ - psi_ - 1.));
        arma::mat AmY  = A_ - Y;
        arma::mat AmYw = AmY;  AmYw.each_col() %= w;
        arma::mat gB_  = Xw.t() * AmY;
        arma::mat gC_  = AmYw.t() * M_ + (A_.t() * (S2_.each_col() % w)) % C_;
        arma::mat gM_  = AmY * C_ + M_;  gM_.each_col() %= w;
        arma::mat gPs_ = arma::diagmat(w) * (0.5 * (S2_ % (1. + A_ * C2_) - 1.));
        arma::vec g    = arma::join_cols(
            arma::join_cols(arma::vectorise(gB_), arma::vectorise(gC_)),
            arma::join_cols(arma::vectorise(gM_), arma::vectorise(gPs_)));
        return {f, g};
    };

    // Initial packed state and evaluation
    arma::vec x = arma::join_cols(
        arma::join_cols(arma::vectorise(B), arma::vectorise(C)),
        arma::join_cols(arma::vectorise(M), arma::vectorise(psi)));
    auto [f_cur, g_cur] = fg(x);

    std::deque<arma::vec> sv, yv;
    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter = 0, last_status = 5;
    const int win = 100;

    for (int it = 0; it < maxiter; it++) {
        objective_vec.push_back(f_cur);
        total_iter++;

        if (it > 0 && converged(f_cur, obj_prev, ftol)) { last_status = 3; break; }
        obj_prev = f_cur;
        if (it > 0 && it % win == 0 && (int)objective_vec.size() >= 2*win) {
            double m1 = *std::min_element(objective_vec.end()-win,   objective_vec.end());
            double m2 = *std::min_element(objective_vec.end()-2*win, objective_vec.end()-win);
            if (converged(m1, m2, ftol)) { last_status = 3; break; }
        }

        // L-BFGS direction
        arma::vec d_lbfgs;
        if (sv.empty()) {
            double gn = arma::norm(g_cur);
            d_lbfgs = (gn > 1e-20) ? arma::vec(-g_cur / gn)
                                    : arma::vec(N, arma::fill::zeros);
        } else {
            d_lbfgs = lbfgs_direction(g_cur, sv, yv);
            if (arma::dot(d_lbfgs, g_cur) >= 0) {
                sv.clear(); yv.clear();
                d_lbfgs = -g_cur / (arma::norm(g_cur) + 1e-20);
            }
        }

        double slope = arma::dot(d_lbfgs, g_cur);
        if (std::abs(slope) < 1e-20) { last_status = 4; break; }

        WolfeStep ws = wolfe_ls(x, d_lbfgs, f_cur, slope, fg);

        // Update L-BFGS history (guarded by curvature condition)
        arma::vec s_new = ws.scale * d_lbfgs;
        arma::vec y_new = ws.g - g_cur;
        double sy = arma::dot(s_new, y_new), ss = arma::dot(s_new, s_new);
        if (sy > 1e-10 * ss && ss > 1e-20) {
            sv.push_back(s_new); yv.push_back(y_new);
            if ((int)sv.size() > m_hist) { sv.pop_front(); yv.pop_front(); }
        }

        x     = x + ws.scale * d_lbfgs;
        f_cur = ws.f;
        g_cur = std::move(ws.g);
    }

    // Unpack final parameters
    B   = arma::reshape(x.subvec(oB,   oC-1  ), d, p);
    C   = arma::reshape(x.subvec(oC,   oM-1  ), p, q);
    M   = arma::reshape(x.subvec(oM,   oPsi-1), n, q);
    psi = arma::reshape(x.subvec(oPsi, N-1   ), n, q);
    S2  = arma::exp(psi);
    C2  = C % C;
    Z   = O + X * B + M * C.t();
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    const double w_bar = arma::accu(w);
    arma::mat nSig = M.t() * (M.each_col() % w) + arma::diagmat(arma::sum(S2.each_col() % w, 0));
    arma::mat Sigma = C * nSig * C.t() / w_bar;
    arma::mat Omega = C * arma::inv_sympd(nSig / w_bar) * C.t();
    arma::vec loglik = arma::sum(Y % Z - A, 1)
                     - 0.5 * arma::sum(M % M + S2 - psi - 1., 1) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("B",     B    ),
        Rcpp::Named("C",     C    ),
        Rcpp::Named("M",     M    ),
        Rcpp::Named("S2",    S2   ),
        Rcpp::Named("Z",     Z    ),
        Rcpp::Named("A",     A    ),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("Ji",    Ji   ),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     last_status  ),
            Rcpp::Named("backend",    "lbfgs"      ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter   )
        ))
    );
}

// ---------------------------------------------------------------------------------------
// VE step only (project): B and C fixed, update (M, ψ).

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_rank(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const arma::mat & B,
    const arma::mat & C,
    const Rcpp::List & config
) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]);
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]);
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]);
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);

    const int    maxiter = config.containsElementNamed("maxeval") ? Rcpp::as<int>   (config["maxeval"])  : 10000;
    const double ftol    = config.containsElementNamed("ftol_in") ? Rcpp::as<double>(config["ftol_in"]) : 1e-9;
    const int    m_hist  = 10;

    const arma::uword n = Y.n_rows, q = M.n_cols;
    const arma::uword oM = 0, oPsi = n*q, N = 2*n*q;
    const arma::mat C2 = C % C;
    const arma::mat XB = X * B;

    // Warm-start ψ
    arma::mat psi = arma::log(S2);
    arma::mat Z   = O + XB + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    psi = arma::clamp(-arma::log(1. + A * C2), -40., 40.);
    S2  = arma::exp(psi);
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    auto fg = [&](const arma::vec & x) -> std::pair<double, arma::vec> {
        arma::mat M_   = arma::reshape(x.subvec(oM,   oPsi-1), n, q);
        arma::mat psi_ = arma::reshape(x.subvec(oPsi, N-1   ), n, q);
        arma::mat S2_  = arma::exp(psi_);
        arma::mat Z_   = O + XB + M_ * C.t();
        arma::mat A_   = arma::exp(Z_ + 0.5 * S2_ * C2.t());
        double f = arma::accu(w.t() * (A_ - Y % Z_))
                 + 0.5 * arma::accu(w.t() * (M_ % M_ + S2_ - psi_ - 1.));
        arma::mat gM_  = (A_ - Y) * C + M_;  gM_.each_col() %= w;
        arma::mat gPs_ = arma::diagmat(w) * (0.5 * (S2_ % (1. + A_ * C2) - 1.));
        return {f, arma::join_cols(arma::vectorise(gM_), arma::vectorise(gPs_))};
    };

    arma::vec x = arma::join_cols(arma::vectorise(M), arma::vectorise(psi));
    auto [f_cur, g_cur] = fg(x);

    std::deque<arma::vec> sv, yv;
    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter = 0;
    const int win = 100;

    for (int it = 0; it < maxiter; it++) {
        objective_vec.push_back(f_cur);
        total_iter++;
        if (it > 0 && converged(f_cur, obj_prev, ftol)) break;
        obj_prev = f_cur;
        if (it > 0 && it % win == 0 && (int)objective_vec.size() >= 2*win) {
            double m1 = *std::min_element(objective_vec.end()-win,   objective_vec.end());
            double m2 = *std::min_element(objective_vec.end()-2*win, objective_vec.end()-win);
            if (converged(m1, m2, ftol)) break;
        }

        arma::vec d;
        if (sv.empty()) {
            double gn = arma::norm(g_cur);
            d = (gn > 1e-20) ? arma::vec(-g_cur / gn)
                              : arma::vec(N, arma::fill::zeros);
        } else {
            d = lbfgs_direction(g_cur, sv, yv);
            if (arma::dot(d, g_cur) >= 0) {
                sv.clear(); yv.clear();
                d = -g_cur / (arma::norm(g_cur) + 1e-20);
            }
        }

        double slope = arma::dot(d, g_cur);
        if (std::abs(slope) < 1e-20) break;

        WolfeStep ws = wolfe_ls(x, d, f_cur, slope, fg);

        arma::vec s_new = ws.scale * d;
        arma::vec y_new = ws.g - g_cur;
        double sy = arma::dot(s_new, y_new), ss = arma::dot(s_new, s_new);
        if (sy > 1e-10 * ss && ss > 1e-20) {
            sv.push_back(s_new); yv.push_back(y_new);
            if ((int)sv.size() > m_hist) { sv.pop_front(); yv.pop_front(); }
        }

        x     = x + ws.scale * d;
        f_cur = ws.f;
        g_cur = std::move(ws.g);
    }

    M   = arma::reshape(x.subvec(oM,   oPsi-1), n, q);
    psi = arma::reshape(x.subvec(oPsi, N-1   ), n, q);
    S2  = arma::exp(psi);
    Z   = O + XB + M * C.t();
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    arma::vec loglik = arma::sum(Y % Z - A, 1)
                     - 0.5 * arma::sum(M % M + S2 - psi - 1., 1) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("M")  = M,
        Rcpp::Named("S2") = S2,
        Rcpp::Named("Ji") = Ji,
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     3            ),
            Rcpp::Named("backend",    "lbfgs"      ),
            Rcpp::Named("objective",  objective_vec),
            Rcpp::Named("iterations", total_iter   )
        ))
    );
}

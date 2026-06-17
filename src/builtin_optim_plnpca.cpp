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
//
// Note: for datasets with large d[1]/sqrt(n) (e.g. barents), joint L-BFGS may
// converge to a local optimum inferior to nlopt-CCSAQ. The nlopt backend is
// recommended when solution quality is the priority.

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
    double f0, double slope0, FG && fg,
    const double c1 = 1e-4, const double c2 = 0.9
) {
    auto zoom = [&](double alo, double ahi, double flo) -> WolfeStep {
        for (int j = 0; j < 20; j++) {
            double a = 0.5 * (alo + ahi);
            auto res = fg(x0 + a * d);
            double fa = res.first; arma::vec ga = res.second;
            if (fa > f0 + c1*a*slope0 || fa >= flo) { ahi = a; }
            else {
                double da = arma::dot(ga, d);
                if (std::abs(da) <= -c2 * slope0) return {a, fa, ga};
                if (da * (ahi - alo) >= 0) ahi = alo;
                alo = a;  flo = fa;
            }
        }
        double a = 0.5 * (alo + ahi);
        auto res = fg(x0 + a * d);
        return {a, res.first, res.second};
    };
    double a = 1.0, ap = 0, fp = f0;
    for (int i = 0; i < 20; i++) {
        auto res = fg(x0 + a * d);
        double fa = res.first; arma::vec ga = res.second;
        if (fa > f0 + c1*a*slope0 || (i > 0 && fa >= fp)) return zoom(ap, a, fp);
        double da = arma::dot(ga, d);
        if (std::abs(da) <= -c2 * slope0) return {a, fa, ga};
        if (da >= 0) return zoom(a, ap, fa);
        ap = a;  fp = fa;
        a  = std::min(2.0 * a, 1e6);
    }
    auto res = fg(x0 + a * d);
    return {a, res.first, res.second};
}

// ---------------------------------------------------------------------------------------
// Shared L-BFGS driver: runs to convergence (relative + windowed-min plateau check),
// used identically by the joint and VE-step rank optimizers below.
// status: 3 = converged, 4 = degenerate slope, 5 = maxiter reached without internal break.

struct LbfgsResult { arma::vec x; std::vector<double> objective_vec; int total_iter; int status; };

template<typename FG>
static LbfgsResult run_lbfgs(
    arma::vec x, FG && fg, int maxiter, double ftol, int m_hist = 10
) {
    constexpr int win = 100;
    auto res0 = fg(x); double f_cur = res0.first; arma::vec g_cur = res0.second;
    const arma::uword N = x.n_elem;

    std::deque<arma::vec> sv, yv;
    std::vector<double> objective_vec;
    double obj_prev = arma::datum::inf;
    int total_iter = 0, status = 5;

    for (int it = 0; it < maxiter; it++) {
        objective_vec.push_back(f_cur);
        total_iter++;

        if (it > 0 && converged(f_cur, obj_prev, ftol)) { status = 3; break; }
        obj_prev = f_cur;
        if (it > 0 && it % win == 0 && (int)objective_vec.size() >= 2*win) {
            double m1 = *std::min_element(objective_vec.end()-win,   objective_vec.end());
            double m2 = *std::min_element(objective_vec.end()-2*win, objective_vec.end()-win);
            if (converged(m1, m2, ftol)) { status = 3; break; }
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
        if (std::abs(slope) < 1e-20) { status = 4; break; }

        WolfeStep ws = wolfe_ls(x, d_lbfgs, f_cur, slope, fg);

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

    return {std::move(x), std::move(objective_vec), total_iter, status};
}

// ---------------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List builtin_optimize_rank(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData D(data);
    arma::mat B  = Rcpp::as<arma::mat>(params["B"]);
    arma::mat C  = Rcpp::as<arma::mat>(params["C"]);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);

    const int    maxiter = config.containsElementNamed("maxeval") ? Rcpp::as<int>   (config["maxeval"])  : 10000;
    const double ftol    = config.containsElementNamed("ftol_in") ? Rcpp::as<double>(config["ftol_in"]) : 1e-9;

    const arma::uword n = D.Y.n_rows, p = D.Y.n_cols, q = M.n_cols, d = B.n_rows;

    // Packed-parameter offsets: x = [vec(B); vec(C); vec(M); vec(ψ)]
    const arma::uword oB = 0, oC = d*p, oM = d*p + p*q, oPsi = d*p + p*q + n*q;
    const arma::uword N  = d*p + p*q + 2*n*q;

    // Warm-start ψ with one fixed-point step
    arma::mat C2  = C % C;
    arma::mat psi = arma::log(S2);
    arma::mat Z   = D.O + D.X * B + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    psi = arma::clamp(-arma::log(1. + A * C2), -40., 40.);
    S2  = arma::exp(psi);
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    const arma::mat Xw = D.X.each_col() % D.w;

    // Joint fg evaluator for all parameters
    auto fg = [&](const arma::vec & x) -> std::pair<double, arma::vec> {
        arma::mat B_   = arma::reshape(x.subvec(oB,   oC-1   ), d, p);
        arma::mat C_   = arma::reshape(x.subvec(oC,   oM-1   ), p, q);
        arma::mat M_   = arma::reshape(x.subvec(oM,   oPsi-1 ), n, q);
        arma::mat psi_ = arma::reshape(x.subvec(oPsi, N-1    ), n, q);
        arma::mat S2_  = arma::exp(psi_);
        arma::mat C2_  = C_ % C_;
        arma::mat Z_   = D.O + D.X * B_ + M_ * C_.t();
        arma::mat A_   = arma::exp(Z_ + 0.5 * S2_ * C2_.t());
        double f = arma::accu(D.w.t() * (A_ - D.Y % Z_))
                 + 0.5 * arma::accu(D.w.t() * (M_ % M_ + S2_ - psi_ - 1.));
        arma::mat AmY  = A_ - D.Y;
        arma::mat AmYw = AmY; AmYw.each_col() %= D.w;
        arma::mat gB_  = Xw.t() * AmY;
        arma::mat gC_  = AmYw.t() * M_ + (A_.t() * (S2_.each_col() % D.w)) % C_;
        arma::mat gM_  = AmY * C_ + M_;  gM_.each_col() %= D.w;
        arma::mat gPs_ = arma::diagmat(D.w) * (0.5 * (S2_ % (1. + A_ * C2_) - 1.));
        arma::vec g    = arma::join_cols(
            arma::join_cols(arma::vectorise(gB_), arma::vectorise(gC_)),
            arma::join_cols(arma::vectorise(gM_), arma::vectorise(gPs_)));
        return {f, g};
    };

    // Initial packed state
    arma::vec x0 = arma::join_cols(
        arma::join_cols(arma::vectorise(B), arma::vectorise(C)),
        arma::join_cols(arma::vectorise(M), arma::vectorise(psi)));

    LbfgsResult res = run_lbfgs(x0, fg, maxiter, ftol);

    // Unpack final parameters
    B   = arma::reshape(res.x.subvec(oB,   oC-1  ), d, p);
    C   = arma::reshape(res.x.subvec(oC,   oM-1  ), p, q);
    M   = arma::reshape(res.x.subvec(oM,   oPsi-1), n, q);
    psi = arma::reshape(res.x.subvec(oPsi, N-1   ), n, q);
    S2  = arma::exp(psi);
    C2  = C % C;
    Z   = D.O + D.X * B + M * C.t();
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    const double w_bar = arma::accu(D.w);
    arma::mat nSig = M.t() * (M.each_col() % D.w) + arma::diagmat(arma::sum(S2.each_col() % D.w, 0));
    arma::mat Sigma = C * nSig * C.t() / w_bar;
    arma::mat Omega = C * arma::inv_sympd(nSig / w_bar) * C.t();
    arma::vec loglik = arma::sum(D.Y % Z - A, 1)
                     - 0.5 * arma::sum(M % M + S2 - psi - 1., 1) + ki(D.Y);

    return Rcpp::List::create(
        Rcpp::Named("B",     B    ),
        Rcpp::Named("C",     C    ),
        Rcpp::Named("M",     M    ),
        Rcpp::Named("S2",    S2   ),
        Rcpp::Named("Z",     Z    ),
        Rcpp::Named("A",     A    ),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("Ji",    loglik),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status",     res.status        ),
            Rcpp::Named("backend",    "lbfgs"           ),
            Rcpp::Named("objective",  res.objective_vec ),
            Rcpp::Named("iterations", res.total_iter    )
        ))
    );
}

// ---------------------------------------------------------------------------------------
// VE step only (project): B and C fixed, update (M, ψ).

// [[Rcpp::export]]
Rcpp::List builtin_optimize_vestep_rank(
    const Rcpp::List & data  ,
    const Rcpp::List & params,
    const Rcpp::List & config
) {
    const PlnData D(data);
    arma::mat M  = Rcpp::as<arma::mat>(params["M"]);
    arma::mat S2 = Rcpp::as<arma::mat>(params["S2"]);
    arma::mat B  = Rcpp::as<arma::mat>(params["B"]);
    arma::mat C  = Rcpp::as<arma::mat>(params["C"]);

    const int    maxiter = config.containsElementNamed("maxeval") ? Rcpp::as<int>   (config["maxeval"])  : 10000;
    const double ftol    = config.containsElementNamed("ftol_in") ? Rcpp::as<double>(config["ftol_in"]) : 1e-9;

    const arma::uword n = D.Y.n_rows, q = M.n_cols;
    const arma::uword oM = 0, oPsi = n*q, N = 2*n*q;
    const arma::mat C2 = C % C;
    const arma::mat XB = D.X * B;

    // Warm-start ψ
    arma::mat psi = arma::log(S2);
    arma::mat Z   = D.O + XB + M * C.t();
    arma::mat A   = arma::exp(Z + 0.5 * S2 * C2.t());
    psi = arma::clamp(-arma::log(1. + A * C2), -40., 40.);
    S2  = arma::exp(psi);
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    auto fg = [&](const arma::vec & x) -> std::pair<double, arma::vec> {
        arma::mat M_   = arma::reshape(x.subvec(oM,   oPsi-1), n, q);
        arma::mat psi_ = arma::reshape(x.subvec(oPsi, N-1   ), n, q);
        arma::mat S2_  = arma::exp(psi_);
        arma::mat Z_   = D.O + XB + M_ * C.t();
        arma::mat A_   = arma::exp(Z_ + 0.5 * S2_ * C2.t());
        double f = arma::accu(D.w.t() * (A_ - D.Y % Z_))
                 + 0.5 * arma::accu(D.w.t() * (M_ % M_ + S2_ - psi_ - 1.));
        arma::mat gM_  = (A_ - D.Y) * C + M_;  gM_.each_col() %= D.w;
        arma::mat gPs_ = arma::diagmat(D.w) * (0.5 * (S2_ % (1. + A_ * C2) - 1.));
        return {f, arma::join_cols(arma::vectorise(gM_), arma::vectorise(gPs_))};
    };

    arma::vec x0 = arma::join_cols(arma::vectorise(M), arma::vectorise(psi));
    LbfgsResult res = run_lbfgs(x0, fg, maxiter, ftol);

    M   = arma::reshape(res.x.subvec(oM,   oPsi-1), n, q);
    psi = arma::reshape(res.x.subvec(oPsi, N-1   ), n, q);
    S2  = arma::exp(psi);
    Z   = D.O + XB + M * C.t();
    A   = arma::exp(Z + 0.5 * S2 * C2.t());

    arma::vec loglik = arma::sum(D.Y % Z - A, 1)
                     - 0.5 * arma::sum(M % M + S2 - psi - 1., 1) + ki(D.Y);

    // status hardcoded to 3 (converged), matching prior behavior: this VE-step
    // never surfaced run_lbfgs's internal status (4 = degenerate slope, 5 = maxiter).
    return make_vestep_result(M, S2, loglik, 3, "lbfgs", res.objective_vec, res.total_iter);
}

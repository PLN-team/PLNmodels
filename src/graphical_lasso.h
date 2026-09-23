#ifndef PLNMODELS_GRAPHICAL_LASSO_H
#define PLNMODELS_GRAPHICAL_LASSO_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

// In-package graphical lasso: minimize over a positive definite Theta
//
//   -log det(Theta) + tr(S Theta) + || L o Theta ||_1
//
// by the block coordinate descent of Friedman, Hastie & Tibshirani (2008),
// in the bookkeeping of Sustik & Calderhead (2012), "GLASSOFAST: An efficient
// GLASSO implementation" (TR-12-29, UT Austin).
//
// This is a direct port of the `glassofast` Fortran subroutine the glassoFast
// package ships, which PLNnetwork and ZIPLNnetwork used to call at every
// M-step. It is shared with the normalblockr package (src/graphical_lasso.h
// there). The reason for bringing it in-house: the Fortran can hang forever
// on a nearly collapsed covariance (entries ~1e-8, as a rank-deficient
// residual covariance produces), in compiled code that never returns to R, so
// that no R-level timeout can stop it either.
//
// The port is faithful except on three points, all deliberate:
//
//  1. When S carries no off-diagonal mass the problem separates exactly and
//     Theta is diagonal. glassoFast returns 1 / max(L_ii, eps) there, dropping
//     S_ii entirely: with an unpenalized diagonal it hands back ~9.09e15
//     instead of 1 / S_ii. That is a bug, not a convention. We return the
//     correct 1 / (S_ii + L_ii).
//  2. The inner coordinate descent is guaranteed to terminate. In the Fortran
//     it is an unbounded `do` loop whose only exit is `dlx < thrLasso`, which
//     is never true once a NaN reaches `dlx`, or when it oscillates at the
//     rounding level, hence the hang. We reject non-finite input and
//     non-positive S_ii + L_ii up front, break on a non-finite `dlx`, and keep
//     a far-off backstop cap, so a caller gets a result it can test.
//     `converged` reports which happened.
//  3. It regularly checks for a user interrupt or an elapsed R time limit
//     (see check_interrupt() below), so it can always be stopped from R.
//
// Fortran guarantees non-aliasing arrays and vectorizes the rank-1 updates
// below on that basis; C++ has to be told, or the compiler keeps them scalar.
// The pointers it is applied to are always distinct allocations (a standalone
// arma::vec and the columns of a matrix).
#if defined(__GNUC__) || defined(__clang__)
  #define PLN_RESTRICT __restrict__
#elif defined(_MSC_VER)
  #define PLN_RESTRICT __restrict
#else
  #define PLN_RESTRICT
#endif

namespace pln_glasso {

// The Fortran's EPS parameter, kept to the digit for comparability.
constexpr double kEps = 1.1e-16;

// Last-resort backstop on the inner coordinate descent, not a working limit:
// termination is guaranteed structurally instead (non-finite input and a
// non-positive S_ii + L_ii are both rejected up front, and a non-finite dlx
// breaks the loop). It is set far above what a well-posed problem needs:
// weak penalties genuinely take thousands of passes.
constexpr int kMaxInner = 500000;

// Stagnation detection on the outer loop.
//
// On an ill-conditioned problem -- typically a rank-deficient residual
// covariance, as PLNnetwork produces when the number of species approaches the
// number of samples -- the sweeps settle into a small limit cycle: `dw` stops
// decreasing and oscillates just above `shr` forever, so `dw <= shr` is never
// met. Sweeping on does not help: the iterates wander inside the cycle instead
// of settling, so stopping at ~1100 sweeps and stopping at 50000 give answers
// the same distance apart (~1e-3 relative on oaks, a couple of borderline
// edges out of ~2900) as either is from the other. glassoFast burns its whole
// budget on these and reports success regardless; we detect the cycle, stop,
// and say so -- the useful output being the diagnosis, not a better solution.
//
// A sweep counts as progress when it improves the best dw seen so far by
// kStallRel; kStallPatience sweeps without progress mean the cycle. The
// patience is set from a census of 246 problems (Gaussian and count data,
// p from 10 to 120, along full penalty paths): problems that do converge never
// went more than 373 sweeps without improving their best dw, while stalled
// ones went 4000+. Anything in between separates them; the value below keeps a
// wide margin on the side that matters -- cutting a converging problem short.
constexpr double kStallRel = 1e-3;
constexpr int    kStallPatience = 1000;

// Interrupt checks are throttled on the work done (in flops, roughly) rather
// than on loop counts, so that how fast a solve can be stopped does not depend
// on the dimension: a single sweep takes seconds at p = 400, microseconds at
// p = 20. 1e7 flops is a few milliseconds.
constexpr double kInterruptWork = 1e7;

inline SEXP check_interrupt_callback(void*) {
  R_CheckUserInterrupt();
  return R_NilValue;
}

// Lets R process a pending user interrupt or an exceeded setTimeLimit(). Unlike
// Rcpp::checkUserInterrupt(), which reports whatever it catches as an
// interrupt, this preserves R's own condition -- a time limit stays the error
// "reached elapsed time limit", which is what R.utils::withTimeout() looks for.
// The longjmp is caught by R_UnwindProtect and rethrown as a C++ exception, so
// the stack unwinds normally; Rcpp resumes the jump at the .Call boundary.
inline void check_interrupt() {
  Rcpp::unwindProtect(check_interrupt_callback, nullptr);
}

// How a solve ended. The two ways of not converging are worth telling apart,
// because they mean different things (see the stagnation note above kStall*):
// `inner_failure` and `degenerate` are numerical trouble, `stalled` is a
// well-behaved solution the stopping rule simply cannot certify.
enum class Status { converged, stalled, max_iter, inner_failure, degenerate };

inline const char * status_name(Status s) {
  switch (s) {
    case Status::converged:     return "converged";
    case Status::stalled:       return "stalled";
    case Status::max_iter:      return "max_iter";
    case Status::inner_failure: return "inner_failure";
    default:                    return "degenerate";
  }
}

struct Result {
  arma::mat W;            // covariance estimate (glassoFast's `w`)
  arma::mat X;            // precision estimate  (glassoFast's `wi`)
  int niter = 0;          // outer sweeps performed
  bool converged = true;  // strictly: the dw <= shr criterion was met
  Status status = Status::converged;
  double delta = 0.0;     // best dw reached, relative to the threshold shr
  std::vector<double> dw_trace; // per-sweep dw, recorded only when asked for
};

// Previous (W, X) used to warm-start a solve. It is only used when it has the
// right size and is finite; anything else silently falls back to a cold start.
//
// Beware that a warm start does not reproduce a cold solve exactly: the outer
// loop stops on `dw <= shr`, how much a whole sweep moved W rather than how far
// W still is from the optimum, so starting closer exits sooner.
struct State {
  arma::mat W;
  arma::mat X;
  bool filled = false;

  bool usable_for(arma::uword n) const {
    return filled && W.n_rows == n && W.n_cols == n && X.n_rows == n && X.n_cols == n
           && W.is_finite() && X.is_finite();
  }
};

// `warm` is used only when it is `usable_for(S.n_rows)`; pass nullptr for a
// cold start. Mirrors glassoFast's defaults (thr = 1e-4, max_iter = 10000).
inline Result solve(const arma::mat& S, const arma::mat& L,
                    double thr = 1e-4, int max_iter = 10000,
                    const State* warm = nullptr, bool trace = false,
                    int stall_patience = kStallPatience) {
  const arma::uword n = S.n_rows;

  Result out;
  out.W.zeros(n, n);
  out.X.zeros(n, n);
  if (n == 0) return out;

  if (!S.is_finite() || !L.is_finite()) {
    out.W.fill(arma::datum::nan);
    out.X.fill(arma::datum::nan);
    out.converged = false;
    out.status = Status::degenerate;
    return out;
  }

  // S_ii + L_ii is what the soft-threshold below divides by, and the only way
  // a finite input could manufacture an infinity there; a non-positive entry
  // is a zero-variance coordinate, i.e. degenerate input. Checked before the
  // separable branch so that both paths answer the same way.
  const arma::vec diag_sum = S.diag() + L.diag();
  if (!arma::all(diag_sum > 0.0)) {
    out.W.fill(arma::datum::nan);
    out.X.fill(arma::datum::nan);
    out.converged = false;
    out.status = Status::degenerate;
    return out;
  }

  arma::mat& W = out.W;
  arma::mat& X = out.X;

  // Total off-diagonal absolute mass of S; sets both convergence thresholds.
  const double off_mass = arma::accu(arma::abs(S)) - arma::accu(arma::abs(S.diag()));

  if (off_mass <= 0.0) {
    // Separable: no coupling to estimate, so the exact solution is diagonal.
    // (See note 1 above -- this is where we part with glassoFast.)
    W.diag() = diag_sum;
    X.diag() = 1.0 / diag_sum;
    return out;
  }

  const double shr = thr * off_mass / static_cast<double>(n - 1);
  const double thr_lasso = std::max(shr / static_cast<double>(n), 2.0 * kEps);

  if (warm != nullptr && warm->usable_for(n)) {
    // The recursion carries X as the negated normalized regression
    // coefficients of each column, not as the precision matrix; a warm start
    // has to be pushed back into that representation first.
    W = warm->W;
    X = warm->X;
    for (arma::uword i = 0; i < n; ++i) {
      const double xii = X(i, i);
      X.col(i) /= -xii;
      X(i, i) = 0.0;
    }
    if (!X.is_finite()) { // a singular warm start (some X_ii == 0)
      W = S;
      X.zeros();
    }
  } else {
    W = S;
    X.zeros();
  }

  arma::vec Wd(n);
  for (arma::uword i = 0; i < n; ++i) {
    Wd(i) = S(i, i) + L(i, i);
    W(i, i) = Wd(i);
  }

  arma::vec WXj(n);
  int iter = 0;
  bool outer_converged = false;
  bool stalled = false;
  double dw_best = arma::datum::inf; // best sweep-to-sweep change so far
  int since_improve = 0;             // sweeps since dw_best last improved
  const double pass_work = static_cast<double>(n) * static_cast<double>(n);
  double work = 0.0; // since the last interrupt check

  // The rest of this function goes through raw column pointers rather than
  // Armadillo element access: this is the hot loop, and every `X(i, j)` would
  // otherwise carry a bounds check.
  const double* const PLN_RESTRICT Wd_p = Wd.memptr();
  double* const PLN_RESTRICT WXj_p = WXj.memptr();

  double dw = 0.0; // kept past the loop so the final value can be reported
  for (iter = 1; iter <= max_iter; ++iter) {
    dw = 0.0;

    for (arma::uword j = 0; j < n; ++j) {
      double* const PLN_RESTRICT Xj = X.colptr(j);
      const double* const PLN_RESTRICT Sj = S.colptr(j);
      const double* const PLN_RESTRICT Lj = L.colptr(j);

      // WXj = W * X.col(j), skipping the zeros X is expected to be full of
      std::fill(WXj_p, WXj_p + n, 0.0);
      for (arma::uword i = 0; i < n; ++i) {
        const double xij = Xj[i];
        if (xij != 0.0) {
          const double* const PLN_RESTRICT Wi = W.colptr(i);
          for (arma::uword k = 0; k < n; ++k) WXj_p[k] += Wi[k] * xij;
        }
      }

      int inner = 0;
      for (;;) {
        double dlx = 0.0;
        for (arma::uword i = 0; i < n; ++i) {
          if (i == j) continue;
          const double a = Sj[i] - WXj_p[i] + Wd_p[i] * Xj[i];
          const double b = std::fabs(a) - Lj[i];
          // soft-threshold; matches Fortran's sign(b, a), which returns +|b|
          // when a is zero (including -0.0)
          const double c = (b > 0.0) ? ((a >= 0.0 ? b : -b) / Wd_p[i]) : 0.0;
          const double delta = c - Xj[i];
          if (delta != 0.0) {
            Xj[i] = c;
            const double* const PLN_RESTRICT Wi = W.colptr(i);
            for (arma::uword k = 0; k < n; ++k) WXj_p[k] += Wi[k] * delta;
            const double ad = std::fabs(delta);
            if (ad > dlx) dlx = ad;
          }
        }
        work += pass_work; // an upper bound: n coordinates, O(n) per update
        if (work >= kInterruptWork) { check_interrupt(); work = 0.0; }
        if (dlx < thr_lasso) break;
        if (!std::isfinite(dlx) || ++inner >= kMaxInner) {
          out.converged = false;
          out.status = Status::inner_failure;
          break;
        }
      }

      WXj_p[j] = Wd_p[j];
      double* const PLN_RESTRICT Wj = W.colptr(j);
      double acc = 0.0;
      for (arma::uword k = 0; k < n; ++k) acc += std::fabs(WXj_p[k] - Wj[k]);
      if (acc > dw) dw = acc;
      for (arma::uword k = 0; k < n; ++k) Wj[k] = WXj_p[k];
      for (arma::uword k = 0; k < n; ++k) W.colptr(k)[j] = WXj_p[k]; // W(j, :)
    }

    if (trace) out.dw_trace.push_back(dw);
    if (dw <= shr) { outer_converged = true; break; }

    if (dw < dw_best * (1.0 - kStallRel)) { dw_best = dw; since_improve = 0; }
    else if (++since_improve >= stall_patience) { stalled = true; break; }
  }

  out.niter = std::min(iter, max_iter);
  if (!outer_converged) out.converged = false;
  out.delta = std::min(dw_best, dw) / shr;
  // an inner-loop failure is the more serious diagnosis and keeps precedence
  if (out.status != Status::inner_failure)
    out.status = outer_converged ? Status::converged
               : (stalled ? Status::stalled : Status::max_iter);

  // Back out the precision matrix from the regression coefficients. X(i,i) is
  // still 0 here, so it drops out of the dot product on its own.
  for (arma::uword i = 0; i < n; ++i) {
    const double tmp = 1.0 / (Wd(i) - arma::dot(X.col(i), W.col(i)));
    X.col(i) *= -tmp;
    X(i, i) = tmp;
  }

  const arma::mat Xt = X.t(); // averaging the two triangles leaves the diagonal as is
  X = 0.5 * (X + Xt);

  return out;
}

} // namespace pln_glasso

#undef PLN_RESTRICT

#endif // PLNMODELS_GRAPHICAL_LASSO_H

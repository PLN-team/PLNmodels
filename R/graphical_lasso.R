#' Graphical Lasso
#'
#' Sparse estimation of a precision matrix by the graphical Lasso, that is the
#' minimization over positive definite matrices \eqn{\Theta} of
#' \deqn{-\log\det(\Theta) + \mathrm{tr}(S\Theta) + \|\rho \circ \Theta\|_1.}
#' This is the solver used internally by [PLNnetwork()] and [ZIPLNnetwork()].
#'
#' The algorithm is the block coordinate descent of Friedman, Hastie and
#' Tibshirani (2008), in the implementation of Sustik and Calderhead (2012):
#' the code is a C++ port of the Fortran routine of the \pkg{glassoFast}
#' package. It first puts the problem on a unit diagonal: for any positive
#' diagonal matrix \eqn{D}, \eqn{\Theta} solves the problem for
#' \eqn{(S, \rho)} if and only if \eqn{D^{-1}\Theta D^{-1}} solves it for
#' \eqn{(DSD, D\rho D)}, and \eqn{D = \mathrm{diag}(S + \rho)^{-1/2}} is
#' used. This change of variables is exact, and on ordinary input the result is
#' the one of \pkg{glassoFast} up to the stopping rule (same support, entries
#' within the tolerance). On a covariance matrix whose variances span orders of
#' magnitude, as the residual covariance of [PLNnetwork()] does, it is what
#' makes the descent converge, in a few sweeps where the unscaled one would
#' cycle, fail in its inner loop, or return an indefinite precision matrix.
#' It also departs from \pkg{glassoFast} on degenerate input:
#' * it always terminates: non-finite input, or a coordinate with
#'   \eqn{S_{ii} + \rho_{ii} \leq 0}, is rejected (the result is filled with
#'   `NA` and `converged` is `FALSE`), and the inner coordinate descent is
#'   bounded, where \pkg{glassoFast} can loop forever on a nearly collapsed
#'   covariance matrix;
#' * failure to converge is reported through `converged` rather than silently;
#' * it can be interrupted from R;
#' * when `S` has no off-diagonal mass, the (diagonal) solution
#'   \eqn{1 / (S_{ii} + \rho_{ii})} is returned, where \pkg{glassoFast}
#'   returns \eqn{1 / \max(\rho_{ii}, \epsilon)};
#' * it detects when it is cycling rather than converging, and stops.
#'
#' That last point is a safeguard. Without the scaling to a unit diagonal, the
#' sweeps could settle into a small limit cycle on an ill-conditioned `S`: the
#' convergence criterion stops decreasing and oscillates just above its
#' threshold forever, while the solution itself no longer moves. \pkg{glassoFast}
#' spends its whole sweep budget on these and reports success regardless; here
#' the cycle would be detected after `stall_patience` sweeps without progress,
#' the solve stopped, and `status` would report `"stalled"`. With the scaling,
#' the covariance matrices on which this was observed (a `oaks`-derived one,
#' residual covariances along [PLNnetwork()] paths) all converge.
#'
#' The remaining statuses are `"max_iter"` (`maxit` reached while still
#' progressing), `"inner_failure"` and `"degenerate"` (numerical trouble, the
#' result may contain `NA`).
#'
#' The precision matrix `wi` is backed out of the regression coefficients of the
#' algorithm, so that it is the inverse of `w` only at the exact solution. Short
#' of it, on an ill-conditioned `S`, it can come out indefinite (whatever the
#' `status`), and a caller using it as a precision matrix then diverges. Its
#' diagonal is then shifted, which keeps the estimated network, so that its
#' smallest eigenvalue is the one of the inverse of `w`. The shift is reported
#' in `shift`, which is `0` otherwise.
#'
#' @param S a symmetric p x p (empirical) covariance matrix.
#' @param rho the penalty: either a non-negative scalar, applied to all entries
#'   (the diagonal included), or a symmetric p x p matrix of non-negative
#'   per-entry penalties (e.g. with a zero diagonal to leave it unpenalized).
#' @param thr convergence threshold, relative to the average absolute
#'   off-diagonal entry of `S` scaled to a unit diagonal (see Details). Default
#'   is `1e-4`, as in \pkg{glassoFast}.
#' @param maxit maximal number of outer sweeps. Default is `10000`, as in
#'   \pkg{glassoFast}.
#' @param w_init,wi_init optional warm start: the `w` and `wi` of a previous
#'   solve, typically at a nearby penalty along a regularization path. Both must
#'   be given, with the dimensions of `S`, to be used. Note that a warm start
#'   stops closer to the starting point than a cold one at the same `thr`, so it
#'   does not reproduce a cold solve exactly.
#' @param trace if `TRUE`, also return the per-sweep convergence criterion, to
#'   diagnose a solve that does not converge. Default is `FALSE`.
#' @param stall_patience number of consecutive sweeps without progress after
#'   which the algorithm concludes that it is cycling and stops (see Details).
#'   `Inf` disables the detection. Default is `1000`.
#'
#' @return a list with components
#' * `w`: the estimated covariance matrix,
#' * `wi`: the estimated precision matrix (symmetric),
#' * `niter`: the number of outer sweeps performed,
#' * `converged`: `TRUE` if the convergence criterion was met,
#' * `status`: how the solve ended, one of `"converged"`, `"stalled"`,
#'   `"max_iter"`, `"inner_failure"` or `"degenerate"` (see Details),
#' * `delta`: the best value reached by the convergence criterion, relative to
#'   the threshold it had to cross. `delta <= 1` means convergence; a stalled
#'   solve typically sits between 1 and 3, that is, just short of it,
#' * `shift`: the value added to the diagonal of `wi` to make it positive
#'   definite, `0` when it already was (see Details),
#' * `dw_trace`: the per-sweep criterion when `trace = TRUE`.
#'
#' @references
#' J. Friedman, T. Hastie and R. Tibshirani (2008). Sparse inverse covariance
#' estimation with the graphical lasso. *Biostatistics*, 9(3), 432--441.
#'
#' M. A. Sustik and B. Calderhead (2012). GLASSOFAST: An efficient GLASSO
#' implementation. UTCS Technical Report TR-12-29, The University of Texas at
#' Austin.
#'
#' @examples
#' data(trichoptera)
#' S <- cov(log1p(as.matrix(trichoptera$Abundance)))
#' fit <- graphical_lasso(S, rho = 0.1)
#' fit$converged
#' sum(fit$wi[upper.tri(fit$wi)] != 0) # number of edges
#'
#' ## penalty weights, leaving the diagonal unpenalized
#' W <- matrix(1, ncol(S), ncol(S)); diag(W) <- 0
#' fit <- graphical_lasso(S, rho = 0.1 * W)
#' @export
graphical_lasso <- function(S, rho, thr = 1e-4, maxit = 10000L,
                            w_init = NULL, wi_init = NULL, trace = FALSE,
                            stall_patience = 1000L) {
  S <- as.matrix(S)
  p <- nrow(S)
  if (!is.numeric(S) || ncol(S) != p)
    stop("`S` must be a square numeric matrix.")
  rho <- as.matrix(rho)
  if (!is.numeric(rho) || anyNA(rho) || any(rho < 0))
    stop("`rho` must be non-negative.")
  if (length(rho) == 1L) {
    rho <- matrix(rho, p, p)
  } else if (!identical(dim(rho), dim(S))) {
    stop("`rho`, when a matrix, must have the same dimensions as `S`.")
  }
  stopifnot(length(thr) == 1L, thr > 0, length(maxit) == 1L, maxit >= 1,
            length(stall_patience) == 1L, !is.na(stall_patience), stall_patience >= 1)
  cpp_graphical_lasso(S, rho, thr, as.integer(maxit),
                      if (is.null(w_init))  NULL else as.matrix(w_init),
                      if (is.null(wi_init)) NULL else as.matrix(wi_init), trace,
                      if (is.finite(stall_patience)) as.integer(stall_patience) else .Machine$integer.max)
}

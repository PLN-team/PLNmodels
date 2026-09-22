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
#' package, and returns the same result on ordinary input. It departs from it
#' on degenerate input only:
#' * it always terminates: non-finite input, or a coordinate with
#'   \eqn{S_{ii} + \rho_{ii} \leq 0}, is rejected (the result is filled with
#'   `NA` and `converged` is `FALSE`), and the inner coordinate descent is
#'   bounded, where \pkg{glassoFast} can loop forever on a nearly collapsed
#'   covariance matrix;
#' * failure to converge is reported through `converged` rather than silently;
#' * it can be interrupted from R;
#' * when `S` has no off-diagonal mass, the (diagonal) solution
#'   \eqn{1 / (S_{ii} + \rho_{ii})} is returned, where \pkg{glassoFast}
#'   returns \eqn{1 / \max(\rho_{ii}, \epsilon)}.
#'
#' @param S a symmetric p x p (empirical) covariance matrix.
#' @param rho the penalty: either a non-negative scalar, applied to all entries
#'   (the diagonal included), or a symmetric p x p matrix of non-negative
#'   per-entry penalties (e.g. with a zero diagonal to leave it unpenalized).
#' @param thr convergence threshold, relative to the average absolute
#'   off-diagonal entry of `S`. Default is `1e-4`, as in \pkg{glassoFast}.
#' @param maxit maximal number of outer sweeps. Default is `10000`, as in
#'   \pkg{glassoFast}.
#' @param w_init,wi_init optional warm start: the `w` and `wi` of a previous
#'   solve, typically at a nearby penalty along a regularization path. Both must
#'   be given, with the dimensions of `S`, to be used. Note that a warm start
#'   stops closer to the starting point than a cold one at the same `thr`, so it
#'   does not reproduce a cold solve exactly.
#'
#' @return a list with components
#' * `w`: the estimated covariance matrix,
#' * `wi`: the estimated precision matrix (symmetric),
#' * `niter`: the number of outer sweeps performed,
#' * `converged`: `TRUE` if the algorithm converged.
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
                            w_init = NULL, wi_init = NULL) {
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
  stopifnot(length(thr) == 1L, thr > 0, length(maxit) == 1L, maxit >= 1)
  cpp_graphical_lasso(S, rho, thr, as.integer(maxit),
                      if (is.null(w_init))  NULL else as.matrix(w_init),
                      if (is.null(wi_init)) NULL else as.matrix(wi_init))
}

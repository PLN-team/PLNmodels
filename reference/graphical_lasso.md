# Graphical Lasso

Sparse estimation of a precision matrix by the graphical Lasso, that is
the minimization over positive definite matrices \\\Theta\\ of
\$\$-\log\det(\Theta) + \mathrm{tr}(S\Theta) + \\\rho \circ
\Theta\\\_1.\$\$ This is the solver used internally by
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
and
[`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md).

## Usage

``` r
graphical_lasso(
  S,
  rho,
  thr = 1e-04,
  maxit = 10000L,
  w_init = NULL,
  wi_init = NULL,
  trace = FALSE,
  stall_patience = 1000L
)
```

## Arguments

- S:

  a symmetric p x p (empirical) covariance matrix.

- rho:

  the penalty: either a non-negative scalar, applied to all entries (the
  diagonal included), or a symmetric p x p matrix of non-negative
  per-entry penalties (e.g. with a zero diagonal to leave it
  unpenalized).

- thr:

  convergence threshold, relative to the average absolute off-diagonal
  entry of `S`. Default is `1e-4`, as in glassoFast.

- maxit:

  maximal number of outer sweeps. Default is `10000`, as in glassoFast.

- w_init, wi_init:

  optional warm start: the `w` and `wi` of a previous solve, typically
  at a nearby penalty along a regularization path. Both must be given,
  with the dimensions of `S`, to be used. Note that a warm start stops
  closer to the starting point than a cold one at the same `thr`, so it
  does not reproduce a cold solve exactly.

- trace:

  if `TRUE`, also return the per-sweep convergence criterion, to
  diagnose a solve that does not converge. Default is `FALSE`.

- stall_patience:

  number of consecutive sweeps without progress after which the
  algorithm concludes that it is cycling and stops (see Details). `Inf`
  disables the detection. Default is `1000`.

## Value

a list with components

- `w`: the estimated covariance matrix,

- `wi`: the estimated precision matrix (symmetric),

- `niter`: the number of outer sweeps performed,

- `converged`: `TRUE` if the convergence criterion was met,

- `status`: how the solve ended, one of `"converged"`, `"stalled"`,
  `"max_iter"`, `"inner_failure"` or `"degenerate"` (see Details),

- `delta`: the best value reached by the convergence criterion, relative
  to the threshold it had to cross. `delta <= 1` means convergence; a
  stalled solve typically sits between 1 and 3, that is, just short of
  it,

- `dw_trace`: the per-sweep criterion when `trace = TRUE`.

## Details

The algorithm is the block coordinate descent of Friedman, Hastie and
Tibshirani (2008), in the implementation of Sustik and Calderhead
(2012): the code is a C++ port of the Fortran routine of the glassoFast
package, and returns the same result on ordinary input. It departs from
it on degenerate input only:

- it always terminates: non-finite input, or a coordinate with
  \\S\_{ii} + \rho\_{ii} \leq 0\\, is rejected (the result is filled
  with `NA` and `converged` is `FALSE`), and the inner coordinate
  descent is bounded, where glassoFast can loop forever on a nearly
  collapsed covariance matrix;

- failure to converge is reported through `converged` rather than
  silently;

- it can be interrupted from R;

- when `S` has no off-diagonal mass, the (diagonal) solution \\1 /
  (S\_{ii} + \rho\_{ii})\\ is returned, where glassoFast returns \\1 /
  \max(\rho\_{ii}, \epsilon)\\;

- it detects when it is cycling rather than converging, and stops.

That last point matters on an ill-conditioned `S` – typically a
rank-deficient covariance, as arises when the number of variables
approaches the number of samples. The sweeps then settle into a small
limit cycle: the convergence criterion stops decreasing and oscillates
just above its threshold forever, while the solution itself no longer
moves. glassoFast spends its whole sweep budget on these and reports
success regardless; here the cycle is detected after `stall_patience`
sweeps without progress, the solve stops, and `status` reports
`"stalled"`.

Stopping early costs nothing, because sweeping on does not buy accuracy:
the iterates wander inside the cycle rather than settle. On a
`oaks`-derived covariance the solve stops after ~1100 sweeps instead of
10000, and both that answer and the 10000-sweep one sit within 1e-3 of a
50000-sweep grind, with the same number of edges and a couple of
borderline ones differing – the amplitude of the cycle, which no budget
removes. The useful consequence of `"stalled"` is thus not a warning
about the solution, but the information that the problem is in that
regime: usually too weak a penalty for a covariance that is (nearly)
rank-deficient.

The remaining statuses are `"max_iter"` (`maxit` reached while still
progressing), `"inner_failure"` and `"degenerate"` (numerical trouble,
the result may contain `NA`).

## References

J. Friedman, T. Hastie and R. Tibshirani (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432–441.

M. A. Sustik and B. Calderhead (2012). GLASSOFAST: An efficient GLASSO
implementation. UTCS Technical Report TR-12-29, The University of Texas
at Austin.

## Examples

``` r
data(trichoptera)
S <- cov(log1p(as.matrix(trichoptera$Abundance)))
fit <- graphical_lasso(S, rho = 0.1)
fit$converged
#> [1] TRUE
sum(fit$wi[upper.tri(fit$wi)] != 0) # number of edges
#> [1] 38

## penalty weights, leaving the diagonal unpenalized
W <- matrix(1, ncol(S), ncol(S)); diag(W) <- 0
fit <- graphical_lasso(S, rho = 0.1 * W)
```

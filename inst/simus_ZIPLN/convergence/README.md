# Convergence of the outer loop of ZIPLN

Scripts behind the SQUAREM acceleration of `ZIPLNfit$optimize()` and the default
`maxit_out = 500` (branch `zipln-convergence`). Cases: `cases.R` (13 real data sets).

- `diagnose.R`, `slow_mode.R`: the loop converges linearly (rates up to 0.996); the slow mode
  is the latent means of the zero counts (EM with much missing information).
- `squarem.R`: prototype of SQUAREM on the state (M, log S2, logit R, B0), with the reference
  optimum of each case (plain iterations to 1e-10, `squarem_adaptive.rds`).
- `profiled.R`, `profiled_R.R`: VE step with Omega and B (and R) profiled, by L-BFGS. Dropped:
  the profiled bound is too ill-conditioned for a first-order solver (no convergence in
  20 000 iterations on barents), the slow mode needs the Newton step of each cell.
- `integrated.R`: SQUAREM as integrated, tolerances 1e-6 to 1e-8.
- `trajectories.R`, `stopping_rules.R`: an extrapolated stopping rule does no better than
  the relative decrease of the objective.
- `baseline.R`, `compare_defaults.R`: candidate defaults against the former ones (maxit_out =
  200, no acceleration); `ftol_out = 1e-6`, `maxit_out = 500` chosen.
- `multistart_check.R`, `abandon_trajectories.R`, `abandon_rules.R`: the two starting points
  with the new defaults; one miss of the abandonment (0.64 on mollusk site|site, the bound
  stalls then improves again), not avoidable by extrapolation.

## Optimization settings of PLNnetwork (builtin backend, correlation scale): VE steps at
## fixed Omega and convergence of the alternating loop, judged on the recovery of a known
## network (same simulations and measures as ../penalty_scale/sim_penalty_scale.R).
##
##   Rscript sim_optim_settings.R <block> <n_rep> <cores> <settings> [lib]
## <settings>: comma-separated names among those below; with [lib], the package installed
## in that library is used (e.g. a former version for "legacy"), otherwise load_all().
## Results in sim_<block>_<settings>.rds.
args  <- commandArgs(TRUE)
block <- args[1]; n_rep <- as.integer(args[2]); cores <- as.integer(args[3])
settings_names <- strsplit(args[4], ",")[[1]]
if (length(args) > 4) {
  library(PLNmodels, lib.loc = args[5])
} else suppressMessages(devtools::load_all(Sys.getenv("PLNMODELS_PATH", "../../.."), quiet = TRUE))
RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)

## maxit_ve = 10000: VE step at fixed Omega to convergence
settings <- list(
  legacy      = list(),   # defaults of the package given in [lib]
  ve1_out5    = list(ve = 1L,     co = list(ftol_out = 1e-5, maxit_out = 20L)),
  ve1_out7    = list(ve = 1L,     co = list(ftol_out = 1e-7, maxit_out = 1000L)),
  full_out5   = list(ve = 10000L, co = list(ftol_out = 1e-5, maxit_out = 20L)),
  full_out7   = list(ve = 10000L, co = list(ftol_out = 1e-7, maxit_out = 1000L)),
  cg5_out5    = list(ve = 10000L, co = list(cg_max = 5L, ftol_out = 1e-5, maxit_out = 20L)),
  cg5_out7    = list(ve = 10000L, co = list(cg_max = 5L, ftol_out = 1e-7, maxit_out = 1000L))
)[settings_names]

## ---------------------------------------------------------------- simulation
graph <- function(type, p) {
  A <- matrix(0, p, p)
  if (type == "erdos") {
    A[upper.tri(A)] <- rbinom(p * (p - 1) / 2, 1, 2 / (p - 1))
  } else if (type == "hub") {
    A <- as.matrix(igraph::as_adjacency_matrix(igraph::sample_pa(p, m = 1, directed = FALSE)))
    A[lower.tri(A)] <- 0
  } else if (type == "cluster") {
    member <- rep(1:4, length.out = p)
    same <- outer(member, member, "==")
    A[upper.tri(A)] <- rbinom(p * (p - 1) / 2, 1, 3 / (p / 4 - 1)) * same[upper.tri(same)]
  } else if (type == "band") {
    A[cbind(1:(p - 1), 2:p)] <- 1
  }
  A + t(A)
}

simulate <- function(cfg, seed) {
  set.seed(seed)
  n <- cfg$n; p <- cfg$p
  A <- graph(cfg$graph, p)
  Omega <- 0.5 * A * matrix(sample(c(-1, 1), p * p, replace = TRUE), p, p)
  Omega[lower.tri(Omega)] <- t(Omega)[lower.tri(Omega)]
  diag(Omega) <- rowSums(abs(Omega)) + 0.1
  sd <- if (cfg$variance == "heterogeneous") exp(runif(p, log(0.5), log(2))) else rep(1, p)
  Sigma <- cov2cor(solve(Omega)) * tcrossprod(sd)            # same support as Omega
  mu <- if (cfg$abundance == "mixed") runif(p, 0, 4) else runif(p, 2, 4)
  x <- rnorm(n); group <- sample(rep(1:3, length.out = n))
  beta <- if (isTRUE(cfg$covariate)) rnorm(p, 0, 0.5) else rep(0, p)
  contaminated <- switch(cfg$contamination, none = integer(0), random = sample(p, 3), groups = sample(p, 8))
  if (cfg$contamination == "groups") mu[contaminated] <- mu[contaminated] + 1.5
  Z <- matrix(rnorm(n * p), n, p) %*% chol(Sigma) + matrix(mu, n, p, byrow = TRUE) + outer(x, beta)
  Y <- matrix(rpois(n * p, exp(Z)), n, p, dimnames = list(paste0("s", 1:n), paste0("sp", 1:p)))
  for (j in contaminated) {
    absent <- if (cfg$contamination == "groups") group != sample(3, 1) else seq_len(n) %in% sample(n, round(0.65 * n))
    Y[absent, j] <- 0
  }
  if (isTRUE(cfg$zero_inflation > 0)) Y[matrix(runif(n * p) < cfg$zero_inflation, n, p)] <- 0
  list(Y = Y, x = x, truth = A == 1, contaminated = contaminated)
}

## ---------------------------------------------------------------- evaluation
evaluate <- function(nets, sim) {
  p <- ncol(sim$Y)
  clean <- setdiff(seq_len(p), sim$contaminated)
  up <- upper.tri(matrix(0, length(clean), length(clean)))
  truth <- sim$truth[clean, clean][up]
  crit <- nets$criteria
  per_model <- t(sapply(nets$models, function(m) {
    A <- as.matrix(m$model_par$Omega) != 0; diag(A) <- FALSE
    est <- A[clean, clean][up]
    c(tp = sum(est & truth), n_clean = sum(est), n_edges = sum(A) / 2,
      cont_edges = sum(A) / 2 - sum(A[clean, clean]) / 2)
  }))
  prec <- per_model[, "tp"] / pmax(per_model[, "n_clean"], 1); rec <- per_model[, "tp"] / max(sum(truth), 1)
  f1 <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
  n_true <- sum(sim$truth) / 2
  s <- which.min(abs(per_model[, "n_edges"] - n_true))       # the model of the size of the true network
  b <- which.max(crit$BIC); e <- which.max(crit$EBIC)
  ## area under the precision-recall curve along the path, up to the recall reached
  o <- order(rec, -prec); aupr <- sum(diff(c(0, rec[o])) * prec[o])
  data.frame(maxF1 = max(f1), F1_size = f1[s], F1_BIC = f1[b], F1_EBIC = f1[e], AUPR = aupr,
             edges_BIC = per_model[b, "n_edges"], edges_EBIC = per_model[e, "n_edges"], edges_true = n_true,
             reached = max(per_model[, "n_edges"]) >= n_true,
             cont_share_size = if (length(sim$contaminated) && per_model[s, "n_edges"] > 0) per_model[s, "cont_edges"] / per_model[s, "n_edges"] else NA,
             floored = length(unique(unlist(lapply(nets$models, function(m) m$floored_species)))),
             degenerate = length(unique(unlist(lapply(nets$models, function(m) m$degenerate_species)))))
}

fit <- function(cfg, sim, s) {
  d <- suppressMessages(prepare_data(sim$Y, data.frame(x = sim$x, row.names = rownames(sim$Y)), offset = "none"))
  f <- if (isTRUE(cfg$covariate)) Abundance ~ 1 + x else Abundance ~ 1
  ctrl <- if (is.null(s$ve)) PLNnetwork_param(trace = 0, min_ratio = 0.05)
          else PLNnetwork_param(trace = 0, min_ratio = 0.05, maxit_ve = s$ve, config_optim = s$co)
  suppressWarnings(PLNnetwork(f, d, control = ctrl))
}

## ---------------------------------------------------------------- designs
base <- list(model = "PLNnetwork", contamination = "none", variance = "unit", abundance = "high")
grid <- function(...) {
  g <- expand.grid(..., stringsAsFactors = FALSE)
  lapply(seq_len(nrow(g)), function(i) modifyList(base, as.list(g[i, , drop = FALSE])))
}
sizes <- list(c(50, 40), c(100, 40), c(200, 40), c(100, 80))
with_sizes <- function(cfgs, sz = sizes)
  unlist(lapply(sz, function(s) lapply(cfgs, function(cfg) modifyList(cfg, list(n = s[1], p = s[2])))), recursive = FALSE)
designs <- list(
  clean = with_sizes(grid(graph = c("erdos", "hub", "cluster", "band"),
                          variance = c("unit", "heterogeneous"), abundance = c("high", "mixed"))),
  contaminated = with_sizes(grid(graph = c("erdos", "hub"), contamination = c("random", "groups")), sizes[1:3]),
  covariate = with_sizes(grid(graph = c("erdos", "hub"), variance = c("unit", "heterogeneous"), covariate = TRUE), sizes[1:3]),
  zi = with_sizes(grid(graph = c("erdos", "hub"), model = c("ZIPLNnetwork", "PLNnetwork"), zero_inflation = 0.2), sizes[2:3]),
  stars = with_sizes(grid(graph = c("erdos", "hub"), variance = c("unit", "heterogeneous")), sizes[2])
)
cfgs <- designs[[block]]
tasks <- expand.grid(cfg = seq_along(cfgs), seed = seq_len(n_rep))

run <- function(i) {
  cfg <- cfgs[[tasks$cfg[i]]]; seed <- tasks$seed[i]
  sim <- simulate(cfg, seed)
  do.call(rbind, lapply(names(settings), function(sn) {
    tt <- system.time(nets <- tryCatch(fit(cfg, sim, settings[[sn]]), error = function(e) conditionMessage(e)))[3]
    if (is.character(nets)) return(data.frame(as.data.frame(cfg[c("model", "graph", "n", "p", "variance", "abundance", "contamination")]),
                                               seed = seed, setting = sn, error = nets))
    res <- evaluate(nets, sim)
    res$pen_loglik <- sum(nets$criteria$pen_loglik)
    res$outer <- sum(sapply(nets$models, function(m) m$optim_par$iterations))
    data.frame(as.data.frame(cfg[c("model", "graph", "n", "p", "variance", "abundance", "contamination")]),
               seed = seed, setting = sn, res, time = unname(tt), error = NA_character_)
  }))
}
out <- parallel::mclapply(seq_len(nrow(tasks)), run, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(failed)) cat(sum(failed), "tasks failed\n")
res <- dplyr::bind_rows(out[!failed])
saveRDS(res, paste0("sim_", block, "_", gsub(",", "_", args[4]), ".rds"))
cat("block", block, ":", length(cfgs), "configurations x", n_rep, "replicates;", nrow(res), "rows;", sum(!is.na(res$error)), "errors\n")

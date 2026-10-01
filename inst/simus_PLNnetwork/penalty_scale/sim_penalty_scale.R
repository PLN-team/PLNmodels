## Which scale for the l1 penalty of network fits by default? Simulations with a
## known network, comparing penalty_scale = "covariance" and "correlation" (the
## latent floor being on, as by default) over
##   - four graphs: erdos (random, mean degree 2), hub (preferential attachment),
##     cluster (4 communities), band (chain);
##   - sample sizes and dimensions (n, p);
##   - latent variances: unit, or heterogeneous (standard deviations from 0.5 to 2);
##   - abundances: high (mean log-count in [2, 4]) or mixed ([0, 4]);
##   - contamination: none, random (3 species absent from 65 % of the samples each),
##     groups (3 groups of samples, 8 species living in one group only);
##   - model: PLNnetwork, PLNnetwork with a covariate, ZIPLNnetwork on zero-inflated data.
##
##   Rscript sim_penalty_scale.R <block> [n_rep] [cores]
## Blocks: clean, contaminated, covariate, zi, stars. Results in sim_<block>.rds.
## Run from this directory, or set PLNMODELS_PATH to the root of the package.
suppressMessages(devtools::load_all(Sys.getenv("PLNMODELS_PATH", "../../.."), quiet = TRUE))
RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
args  <- commandArgs(TRUE)
block <- if (length(args)) args[1] else "clean"
n_rep <- if (length(args) > 1) as.integer(args[2]) else 20
cores <- if (length(args) > 2) as.integer(args[3]) else 16
scales <- c("covariance", "correlation")

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

fit <- function(cfg, sim, scale) {
  d <- suppressMessages(prepare_data(sim$Y, data.frame(x = sim$x, row.names = rownames(sim$Y)), offset = "none"))
  f <- if (isTRUE(cfg$covariate)) Abundance ~ 1 + x else Abundance ~ 1
  if (cfg$model == "ZIPLNnetwork")
    suppressWarnings(ZIPLNnetwork(f, d, control = ZIPLNnetwork_param(trace = 0, min_ratio = 0.05, penalty_scale = scale)))
  else
    suppressWarnings(PLNnetwork(f, d, control = PLNnetwork_param(trace = 0, min_ratio = 0.05, penalty_scale = scale)))
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
  do.call(rbind, lapply(scales, function(scale) {
    tt <- system.time(nets <- tryCatch(fit(cfg, sim, scale), error = function(e) conditionMessage(e)))[3]
    if (is.character(nets)) return(data.frame(as.data.frame(cfg[c("model", "graph", "n", "p", "variance", "abundance", "contamination")]),
                                               seed = seed, scale = scale, error = nets))
    res <- evaluate(nets, sim)
    if (block == "stars") {
      ## stability selection, on 20 subsamples as by default
      invisible(capture.output(suppressWarnings(nets$stability_selection())))
      best <- tryCatch(suppressWarnings(nets$getBestModel("StARS")), error = function(e) NULL)
      if (!is.null(best)) {
        A <- as.matrix(best$model_par$Omega) != 0; diag(A) <- FALSE
        est <- A[upper.tri(A)]; truth <- sim$truth[upper.tri(sim$truth)]
        pr <- sum(est & truth) / max(sum(est), 1); rc <- sum(est & truth) / max(sum(truth), 1)
        res$F1_StARS <- if (pr + rc > 0) 2 * pr * rc / (pr + rc) else 0; res$edges_StARS <- sum(est)
      } else { res$F1_StARS <- NA; res$edges_StARS <- NA }
    }
    data.frame(as.data.frame(cfg[c("model", "graph", "n", "p", "variance", "abundance", "contamination")]),
               covariate = isTRUE(cfg$covariate), zero_inflation = if (is.null(cfg$zero_inflation)) 0 else cfg$zero_inflation,
               seed = seed, scale = scale, res, time = unname(tt), error = NA_character_)
  }))
}
out <- parallel::mclapply(seq_len(nrow(tasks)), run, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(failed)) cat(sum(failed), "tasks failed:", if (any(vapply(out, inherits, logical(1), "try-error"))) as.character(out[[which(failed)[1]]]) else "NULL result", "\n")
res <- dplyr::bind_rows(out[!failed])
saveRDS(res, paste0("sim_", block, ".rds"))
cat("block", block, ":", length(cfgs), "configurations x", n_rep, "replicates;", nrow(res), "rows;", sum(!is.na(res$error)), "errors\n")

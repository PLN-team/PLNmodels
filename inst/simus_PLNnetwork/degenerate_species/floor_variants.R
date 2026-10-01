## Variants of the floor on the variational means, as prototypes: the package
## is loaded as is, and its internal project_latent_floor() is replaced in the
## session by a version whose bound depends on the option proto.floor_mode:
##   "fixed"    : M_ij >= log(eps) - O_ij                    (the package's floor)
##   "relative" : M_ij >= log(c) + min over the positive counts of species j of
##                log(Y_ij) - O_ij: a fraction c of its smallest observed abundance
##   "targeted" : the fixed bound, on the species flagged as degenerate only.
##                A species is flagged when its residual variance exceeds the
##                threshold of $degenerate_species (or the option proto.latch),
##                and stays so down the path.
## Run from this directory, or set PLNMODELS_PATH to the root of the package.
suppressMessages(devtools::load_all(Sys.getenv("PLNMODELS_PATH", "../../.."), quiet = TRUE))
elbo_pln <- PLNmodels:::elbo_fixed_precision

floor_state <- new.env()
floor_bound <- function(optim_out, data, floor) {
  n <- nrow(data$Y); p <- ncol(data$Y)
  switch(getOption("proto.floor_mode", "fixed"),
    fixed    = log(floor) - data$O,
    relative = {
      rate <- log(data$Y) - data$O; rate[data$Y == 0] <- Inf
      smallest <- apply(rate, 2, min)                        # Inf for a species never observed
      matrix(ifelse(is.finite(smallest), log(floor) + smallest, -Inf), n, p, byrow = TRUE)
    },
    targeted = {
      flagged <- diag(as.matrix(optim_out$Sigma)) > getOption("proto.latch", PLNmodels:::latent_variance_threshold())
      if (!is.null(floor_state$flagged)) flagged <- flagged | floor_state$flagged
      floor_state$flagged <- flagged
      bound <- log(floor) - data$O
      bound[, !flagged] <- -Inf
      bound
    })
}

## project_latent_floor() of the package, with the bound above
proto_project <- function(optim_out, data, Omega, floor) {
  bound   <- floor_bound(optim_out, data, floor)
  clipped <- optim_out$M < bound
  optim_out$n_floor <- sum(clipped)
  optim_out$M[clipped] <- bound[clipped]
  if (optim_out$n_floor > 0) {
    z  <- (data$O + optim_out$M)[clipped]
    om <- matrix(diag(as.matrix(Omega)), nrow(clipped), ncol(clipped), byrow = TRUE)[clipped]
    lo <- rep(-30, length(z)); hi <- rep(30, length(z))
    for (k in seq_len(60)) {
      mid <- (lo + hi) / 2
      positive <- 1 / exp(mid) - om - exp(z + exp(mid) / 2) > 0
      lo <- ifelse(positive, mid, lo); hi <- ifelse(positive, hi, mid)
    }
    optim_out$S2[clipped] <- exp((lo + hi) / 2)
    w <- data$w
    optim_out$B <- solve(crossprod(data$X, w * data$X), crossprod(data$X, w * optim_out$M))
    R <- optim_out$M - data$X %*% optim_out$B
    optim_out$Sigma <- (crossprod(R, w * R) + diag(colSums(w * optim_out$S2), ncol(R))) / sum(w)
    optim_out$Z <- data$O + optim_out$M
    optim_out$A <- exp(optim_out$Z + .5 * optim_out$S2)
  }
  optim_out$Ji <- elbo_pln(data$Y, data$X, data$O, optim_out$B, optim_out$M, optim_out$S2, Omega)
  optim_out
}
local({
  ns <- asNamespace("PLNmodels")
  unlockBinding("project_latent_floor", ns)
  assign("project_latent_floor", proto_project, envir = ns)
})

## A path of PLNnetwork with one of the variants (mode = "none" for no floor)
fit_variant <- function(formula, data, mode = "none", floor = NULL, scale = "correlation", min_ratio = 0.05, latch = NULL) {
  options(proto.floor_mode = if (mode == "none") "fixed" else mode, proto.latch = latch)
  floor_state$flagged <- NULL
  on.exit(options(proto.floor_mode = NULL, proto.latch = NULL))
  nets <- suppressWarnings(PLNnetwork(formula, data, control = PLNnetwork_param(
    trace = 0, min_ratio = min_ratio, penalty_scale = scale, latent_floor = if (mode == "none") NULL else floor)))
  attr(nets, "floored_species") <- if (mode == "targeted") sum(floor_state$flagged) else NA
  nets
}

variants <- list(
  "none"            = function(f, d) fit_variant(f, d),
  "fixed 1e-2"      = function(f, d) fit_variant(f, d, "fixed", 1e-2),
  "fixed 1e-3"      = function(f, d) fit_variant(f, d, "fixed", 1e-3),
  "relative 1e-1"   = function(f, d) fit_variant(f, d, "relative", 1e-1),
  "relative 1e-2"   = function(f, d) fit_variant(f, d, "relative", 1e-2),
  "targeted 1e-2"   = function(f, d) fit_variant(f, d, "targeted", 1e-2),
  "targeted 1e-3"   = function(f, d) fit_variant(f, d, "targeted", 1e-3),
  "targeted 1e-3 latch 50" = function(f, d) fit_variant(f, d, "targeted", 1e-3, latch = 50))
## a subset can be selected, e.g. VARIANTS="none,targeted 1e-3"
if (nzchar(Sys.getenv("VARIANTS"))) variants <- variants[strsplit(Sys.getenv("VARIANTS"), ",")[[1]]]

degrees <- function(fit) { A <- as.matrix(fit$model_par$Omega) != 0; diag(A) <- FALSE; rowSums(A) }

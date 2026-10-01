## Non-invasive prototypes for the degenerate-species question: the package is
## loaded as is, and PLNnetworkfit$optimize is replaced in the session by a
## version that reads two options:
##   proto.floor : NULL, or eps > 0: the variational means are kept above
##                 log(eps) - O (expected count of a zero cell >= eps)
##   proto.corr  : FALSE/TRUE: l1 penalty on the correlation scale,
##                 rho_ij = lambda * w_ij * sqrt(S_ii S_jj)
## run from this directory, or set PLNMODELS_PATH to the root of the package
suppressMessages(devtools::load_all(Sys.getenv("PLNMODELS_PATH", "../../.."), quiet = TRUE))

## ELBO of each sample, from the parameters (unit weights)
elbo_pln <- function(Y, X, O, B, M, S2, Omega) {
  Omega <- as.matrix(Omega)
  Z <- O + M; A <- exp(Z + .5 * S2); R <- M - X %*% B
  rowSums(Y * Z - A + .5 * log(S2) - lgamma(Y + 1)) +
    .5 * as.numeric(determinant(Omega)$modulus) -
    .5 * rowSums((R %*% Omega) * R) - .5 * as.vector(S2 %*% diag(Omega)) + .5 * ncol(Y)
}

proto_optimize <- function(data, config) {
  eps      <- getOption("proto.floor")
  use_corr <- isTRUE(getOption("proto.corr"))
  nrm  <- normalize_covariates(data$X)
  B_sc <- sweep(private$B, 1, nrm$scales, "*")
  cond <- FALSE; iter <- 0
  objective <- numeric(config$maxit_em); convergence <- numeric(config$maxit_em)
  objective.old <- -self$loglik
  inner_config <- config
  if (!is.null(config$maxit_ve)) {
    if (config$backend == "builtin") inner_config$maxit_em <- as.integer(config$maxit_ve)
    else                             inner_config$maxeval  <- as.integer(config$maxit_ve) * 10L
  }
  args <- list(data = list(Y = data$Y, X = nrm$X_sc, O = data$O, w = data$w),
               params = list(B = B_sc, M = private$M, S2 = private$S2), config = inner_config)
  M_res_init <- private$M - nrm$X_sc %*% B_sc
  private$Sigma <- crossprod(M_res_init)/self$n + diag(colMeans(private$S2), self$p, self$p)
  last_glasso <- NULL; failure <- NULL; n_clipped <- 0
  floor_on_count <- isTRUE(getOption("proto.floor_count"))  # bound the expected count rather than M
  while (!cond) {
    iter <- iter + 1
    rho <- self$penalty * self$penalty_weights
    if (use_corr) rho <- rho * tcrossprod(sqrt(diag(as.matrix(private$Sigma))))
    glasso_out <- graphical_lasso(private$Sigma, rho = rho)
    if (!all(is.finite(glasso_out$wi))) { failure <- "glasso"; break }
    previous_Omega <- private$Omega
    private$Omega <- args$params$Omega <- Matrix::symmpart(glasso_out$wi)
    optim_out <- do.call(private$optimizer$main, args)
    if (!is.null(eps)) {
      ## either M >= log(eps) - O, or the expected count exp(O + M + S2/2) >= eps
      floor_M <- log(eps) - data$O - if (floor_on_count) .5 * optim_out$S2 else 0
      clipped <- optim_out$M < floor_M
      n_clipped <- sum(clipped)
      if (n_clipped > 0) {
        optim_out$M[clipped] <- floor_M[clipped]
        ## the variational variance of a bounded cell is set to its optimum given
        ## M, the root of 1/s - Omega_jj - exp(O + M + s/2) = 0 (decreasing in s)
        if (isTRUE(getOption("proto.refit_S2", TRUE))) {
          z  <- (data$O + optim_out$M)[clipped]
          om <- matrix(diag(as.matrix(args$params$Omega)), nrow(optim_out$M), ncol(optim_out$M), byrow = TRUE)[clipped]
          lo <- rep(-30, length(z)); hi <- rep(30, length(z))  # bisection on log(s)
          for (k in 1:60) {
            mid <- (lo + hi) / 2; s2 <- exp(mid)
            f <- 1 / s2 - om - exp(z + s2 / 2)
            lo <- ifelse(f > 0, mid, lo); hi <- ifelse(f > 0, hi, mid)
          }
          optim_out$S2[clipped] <- exp((lo + hi) / 2)
        }
        X <- args$data$X
        optim_out$B <- solve(crossprod(X), crossprod(X, optim_out$M))
        R <- optim_out$M - X %*% optim_out$B
        optim_out$Sigma <- crossprod(R) / nrow(R) + diag(colMeans(optim_out$S2), ncol(R))
        optim_out$Z <- data$O + optim_out$M
        optim_out$A <- exp(optim_out$Z + .5 * optim_out$S2)
        optim_out$Ji <- elbo_pln(data$Y, X, data$O, optim_out$B, optim_out$M, optim_out$S2, args$params$Omega)
      }
    }
    new_objective <- -sum(optim_out$Ji)
    if (!is.finite(new_objective)) { private$Omega <- previous_Omega; failure <- "objective"; break }
    do.call(self$update, optim_out)
    last_glasso <- glasso_out
    objective[iter] <- new_objective
    convergence[iter] <- abs(objective[iter] - objective.old)/abs(objective[iter])
    if ((convergence[iter] < config$ftol_em) | (iter >= config$maxit_em)) cond <- TRUE
    args$params <- list(B = private$B, M = private$M, S2 = private$S2)
    objective.old <- objective[iter]
  }
  if (!is.null(failure)) iter <- iter - 1
  private$B <- sweep(private$B, 1, nrm$scales, "/")
  if (!is.null(last_glasso)) private$Sigma <- Matrix::symmpart(last_glasso$w)
  private$monitoring$objective   <- objective[seq_len(iter)]
  private$monitoring$iterations  <- iter
  private$monitoring$failure     <- failure
  private$monitoring$n_clipped   <- n_clipped
  if (is.null(last_glasso)) private$Ji <- rep(NA_real_, self$n)
}
PLNnetworkfit$set("public", "optimize", proto_optimize, overwrite = TRUE)

resid_offdiag_max <- function(fit, X, corr = FALSE) {
  R <- fit$var_par$M - X %*% fit$model_par$B
  S <- crossprod(R) / nrow(R) + diag(colMeans(fit$var_par$S2), ncol(R))
  if (corr) S <- cov2cor(S)
  max(abs(S[upper.tri(S)]))
}

## One penalty path, started from the diagonal PLN (the empty network), with a
## grid from the top (largest off-diagonal residual covariance, or correlation)
## down to min_ratio times it.
fit_path <- function(formula, data, X, floor = NULL, corr = FALSE, weights = NULL,
                     n_penalties = 30, min_ratio = 0.1, inception = NULL, floor_count = FALSE) {
  options(proto.floor = floor, proto.corr = corr, proto.floor_count = floor_count)
  on.exit(options(proto.floor = NULL, proto.corr = FALSE, proto.floor_count = FALSE))
  if (is.null(inception))
    inception <- PLN(formula, data, control = PLN_param(covariance = "diagonal", trace = 0))
  R <- inception$var_par$M - X %*% inception$model_par$B
  S <- crossprod(R) / nrow(R) + diag(colMeans(inception$var_par$S2), ncol(R))
  if (corr) S <- cov2cor(S)
  W <- if (is.null(weights)) matrix(1, ncol(S), ncol(S)) else weights
  ok <- upper.tri(S) & W > 0 & is.finite(W)
  top <- max(abs(S[ok]) / W[ok])
  pens <- 10^seq(log10(top), log10(top * min_ratio), length.out = n_penalties)
  suppressWarnings(PLNnetwork(formula, data, penalties = pens,
    control = PLNnetwork_param(trace = 0, inception = inception, penalty_weights = weights)))
}

## Exclusion from the network: the edges of the species flagged as degenerate
## get a prohibitive penalty weight (the species stays in the model, as an
## isolated node), and the path is refitted until no new species is flagged.
fit_exclusion <- function(formula, data, X, max_rounds = 10, select = c("path", "BIC"), ...) {
  select <- match.arg(select)
  p <- ncol(model.response(model.frame(formula, data)))
  species <- colnames(model.response(model.frame(formula, data)))
  excluded <- character(0); rounds <- 0
  inception <- PLN(formula, data, control = PLN_param(covariance = "diagonal", trace = 0))
  repeat {
    rounds <- rounds + 1
    W <- matrix(1, p, p, dimnames = list(species, species))
    W[excluded, ] <- 1e8; W[, excluded] <- 1e8; diag(W) <- 1
    nets <- fit_path(formula, data, X, weights = W, inception = inception, ...)
    flagged <- if (select == "path") unique(unlist(lapply(nets$models, function(m) m$degenerate_species)))
               else suppressWarnings(nets$getBestModel("BIC"))$degenerate_species
    new <- setdiff(flagged, excluded)
    if (length(new) == 0 || rounds >= max_rounds) break
    excluded <- c(excluded, new)
  }
  attr(nets, "excluded") <- excluded; attr(nets, "rounds") <- rounds
  nets
}

degrees <- function(fit) { A <- as.matrix(fit$model_par$Omega) != 0; diag(A) <- FALSE; rowSums(A) }

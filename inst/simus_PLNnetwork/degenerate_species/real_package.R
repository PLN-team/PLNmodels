## Real data with the options of the package (penalty_scale, latent_floor):
## degenerate species along the path, and hubs at a fixed network size.
## Run from this directory, or set PLNMODELS_PATH to the root of the package.
suppressMessages(devtools::load_all(Sys.getenv("PLNMODELS_PATH", "../../.."), quiet = TRUE))
data(barents); data(oaks); data(mollusk); data(trichoptera)
mol <- suppressMessages(suppressWarnings(prepare_data(mollusk$Abundance, mollusk$Covariate)))
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
cases <- list(
  list("oaks ~1",        Abundance ~ 1 + offset(log(Offset)), oaks),
  list("oaks ~tree",     Abundance ~ 1 + tree + offset(log(Offset)), oaks),
  list("mollusk ~1",     Abundance ~ 1 + offset(log(Offset)), mol),
  list("barents ~1",     Abundance ~ 1 + offset(log(Offset)), barents),
  list("trichoptera ~1", Abundance ~ 1, tri))
settings <- list(
  covariance          = list(),
  "covariance+floor"  = list(latent_floor = 1e-3),
  correlation         = list(penalty_scale = "correlation"),
  "correlation+floor" = list(penalty_scale = "correlation", latent_floor = 1e-3))
degrees <- function(fit) { A <- as.matrix(fit$model_par$Omega) != 0; diag(A) <- FALSE; rowSums(A) }
out <- NULL
for (cs in cases) for (s in names(settings)) {
  p <- ncol(cs[[3]]$Abundance)
  tt <- system.time(nets <- suppressWarnings(PLNnetwork(cs[[2]], cs[[3]],
    control = do.call(PLNnetwork_param, c(list(trace = 0, min_ratio = 0.05), settings[[s]])))))[3]
  n_edges <- vapply(nets$models, function(m) m$n_edges, numeric(1))
  m <- nets$models[[which.min(abs(n_edges - p))]]          # about p edges: an average degree of 2
  out <- rbind(out, data.frame(data = cs[[1]], setting = s,
    degenerate_path = length(unique(unlist(lapply(nets$models, function(m) m$degenerate_species)))),
    max_Sigma_path = signif(max(vapply(nets$models, function(m) max(diag(as.matrix(m$model_par$Sigma))), numeric(1))), 3),
    edges = m$n_edges, max_degree = max(degrees(m)), median_degree = median(degrees(m)),
    cells_at_floor = round(m$optim_par$n_floor / length(cs[[3]]$Abundance), 3),
    floored_species = length(unique(unlist(lapply(nets$models, function(m) m$floored_species)))),
    time = round(unname(tt), 1)))
  saveRDS(out, "real_package.rds")
}
options(width = 200); print(out, row.names = FALSE)

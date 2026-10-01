## Real data, at a fixed network size (the BIC is not a reliable judge here):
## who are the hubs, and are the structurally absent species among them?
source("proto.R")
data(barents); data(oaks); data(mollusk)
mol <- suppressMessages(suppressWarnings(prepare_data(mollusk$Abundance, mollusk$Covariate)))
cases <- list(
  list("oaks ~1",    Abundance ~ 1 + offset(log(Offset)), oaks, ~ 1, oaks["tree"]),
  list("oaks ~tree", Abundance ~ 1 + tree + offset(log(Offset)), oaks, ~ 1 + tree, oaks["tree"]),
  list("mollusk ~1", Abundance ~ 1 + offset(log(Offset)), mol, ~ 1, mol["site"]),
  list("barents ~1", Abundance ~ 1 + offset(log(Offset)), barents, ~ 1, NULL))
methods <- list(
  baseline          = function(f, d, X) fit_path(f, d, X, min_ratio = 0.05),
  "floor_1e-3"      = function(f, d, X) fit_path(f, d, X, floor = 1e-3, min_ratio = 0.05),
  corr_penalty      = function(f, d, X) fit_path(f, d, X, corr = TRUE, min_ratio = 0.05),
  "corr+floor_1e-3" = function(f, d, X) fit_path(f, d, X, floor = 1e-3, corr = TRUE, min_ratio = 0.05))
out <- NULL
for (cs in cases) {
  d <- cs[[3]]; X <- model.matrix(cs[[4]], d); Y <- as.matrix(d$Abundance); p <- ncol(Y)
  structural <- if (is.null(cs[[5]])) character(0) else unique(structural_zeros(Y, cs[[5]])$species)
  for (meth in names(methods)) {
    nets <- methods[[meth]](cs[[2]], d, X)
    n_edges <- sapply(nets$models, function(m) m$n_edges)
    for (target in c(p, 2 * p)) {                 # an average degree of 2, then 4
      m <- nets$models[[which.min(abs(n_edges - target))]]
      deg <- sort(degrees(m), decreasing = TRUE)
      out <- rbind(out, data.frame(data = cs[[1]], method = meth, target = target, edges = m$n_edges,
        max_degree = deg[1], top3_share = round(sum(deg[1:3]) / max(2 * m$n_edges, 1), 2),
        structural_share = round(if (length(structural)) sum(degrees(m)[structural]) / max(2 * m$n_edges, 1) else NA, 2),
        expected_share = round(length(structural) / p, 2),
        degenerate = length(m$degenerate_species), max_Sigma = signif(max(diag(as.matrix(m$model_par$Sigma))), 3),
        top3 = paste(names(deg)[1:3], collapse = ",")))
    }
    saveRDS(out, "real_fixed_size.rds")
  }
}
options(width = 250); print(out, row.names = FALSE)

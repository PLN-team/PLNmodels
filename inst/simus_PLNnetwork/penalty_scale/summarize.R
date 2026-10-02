## Paired comparison of the two penalty scales, from the results of
## sim_penalty_scale.R: Rscript summarize.R [block ...]
suppressMessages(library(dplyr))
blocks <- commandArgs(TRUE)
if (!length(blocks)) blocks <- c("clean", "contaminated", "covariate", "zi", "stars")
options(width = 230, dplyr.summarise.inform = FALSE)

## one row per (configuration, seed): correlation minus covariance
paired <- function(res, metrics) {
  keys <- c("model", "graph", "n", "p", "variance", "abundance", "contamination", "covariate", "zero_inflation", "seed")
  cov <- res[res$scale == "covariance", c(keys, metrics)]
  cor <- res[res$scale == "correlation", c(keys, metrics)]
  inner_join(cov, cor, by = keys, suffix = c(".cov", ".cor"))
}
stat <- function(d, m) {
  a <- d[[paste0(m, ".cov")]]; b <- d[[paste0(m, ".cor")]]; ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  data.frame(metric = m, covariance = mean(a), correlation = mean(b), diff = mean(b - a),
             wins = sum(b > a + 1e-12), losses = sum(a > b + 1e-12), ties = sum(abs(a - b) <= 1e-12),
             p_value = if (all(b == a)) 1 else suppressWarnings(wilcox.test(b, a, paired = TRUE)$p.value))
}

for (block in blocks) {
  file <- paste0("sim_", block, ".rds")
  if (!file.exists(file)) next
  res <- readRDS(file)
  cat("\n\n================ block:", block, "-", nrow(res) / 2, "fits per scale;", sum(!is.na(res$error)), "errors\n")
  res <- res[is.na(res$error), ]
  metrics <- c("maxF1", "AUPR", "F1_size", "F1_BIC", "F1_EBIC", "time")
  if (block == "contaminated") metrics <- c(metrics, "cont_share_size")
  if (block == "stars") metrics <- c(metrics, "F1_StARS")
  d <- paired(res, metrics)

  cat("\n-- overall (mean over all configurations and replicates; wins = replicates where correlation is better)\n")
  print(do.call(rbind, lapply(metrics, function(m) stat(d, m))) %>% mutate(across(c(covariance, correlation, diff), ~ round(.x, 3)), p_value = signif(p_value, 2)), row.names = FALSE)

  by_factor <- function(factor, metric = "F1_size") {
    out <- d %>% group_by(across(all_of(factor))) %>%
      group_modify(~ stat(.x, metric)) %>% ungroup() %>%
      mutate(across(c(covariance, correlation, diff), ~ round(.x, 3)), p_value = signif(p_value, 2))
    cat("\n--", metric, "by", paste(factor, collapse = " x "), "\n"); print(as.data.frame(out), row.names = FALSE)
  }
  factors <- switch(block,
    clean = list("graph", c("n", "p"), "variance", "abundance", c("variance", "abundance")),
    contaminated = list(c("contamination", "n"), "graph"),
    covariate = list(c("variance", "n"), "graph"),
    zi = list(c("model", "n"), "graph"),
    stars = list(c("graph", "variance")))
  for (f in factors) by_factor(f)
  if (block %in% c("clean", "covariate")) for (f in factors[1:2]) by_factor(f, "maxF1")
  if (block == "zi") by_factor(c("model", "n"), "maxF1")
  if (block == "stars") { by_factor(c("graph", "variance"), "F1_StARS"); by_factor(c("graph", "variance"), "F1_BIC") }
  if (block == "contaminated") by_factor(c("contamination", "n"), "cont_share_size")

  cat("\n-- configurations where correlation is worse on F1_size (mean difference < -0.02)\n")
  worse <- d %>% group_by(model, graph, n, p, variance, abundance, contamination) %>%
    summarise(covariance = mean(F1_size.cov), correlation = mean(F1_size.cor), diff = mean(F1_size.cor - F1_size.cov)) %>%
    ungroup() %>% filter(diff < -0.02) %>% arrange(diff)
  if (nrow(worse)) print(as.data.frame(worse %>% mutate(across(c(covariance, correlation, diff), ~ round(.x, 3)))), row.names = FALSE) else cat("none\n")
  cat("\n-- path reached the size of the true network (share of fits): covariance",
      round(mean(res$reached[res$scale == "covariance"]), 3), " correlation", round(mean(res$reached[res$scale == "correlation"]), 3), "\n")
  cat("-- species bounded by the floor (mean per fit): covariance",
      round(mean(res$floored[res$scale == "covariance"]), 2), " correlation", round(mean(res$floored[res$scale == "correlation"]), 2), "\n")
  cat("-- mean number of edges at BIC / EBIC / true: covariance",
      round(mean(res$edges_BIC[res$scale == "covariance"])), "/", round(mean(res$edges_EBIC[res$scale == "covariance"])), " correlation",
      round(mean(res$edges_BIC[res$scale == "correlation"])), "/", round(mean(res$edges_EBIC[res$scale == "correlation"])), " true", round(mean(res$edges_true)), "\n")
}

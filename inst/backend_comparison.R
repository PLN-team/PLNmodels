## ============================================================
## Backend comparison: homemade Newton vs nlopt/CCSAQ
## Metrics: computation time, iterations, final loglik
## Datasets: trichoptera, barents, mollusk, oaks, microcosm, scRNA
## Covariances: full, diagonal, spherical (including scRNA full)
## Output: inst/benchmark/
## ============================================================

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

ctrl_newton <- function(cov) PLN_param(backend = "homemade", covariance = cov, trace = 0)
ctrl_nlopt  <- function(cov) PLN_param(backend = "nlopt",    covariance = cov, trace = 0)

## ---- Helper: fit one model with timing, return summary row ----
fit_timed <- function(formula, data, cov, backend_ctrl, backend_name, label) {
  t0 <- proc.time()
  m  <- tryCatch(
    PLN(formula, data = data, control = backend_ctrl(cov)),
    error = function(e) NULL
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]
  if (is.null(m)) {
    return(data.frame(
      label = label, backend = backend_name, covariance = cov,
      time_s = elapsed, n_iter = NA_integer_, loglik = NA_real_,
      converged = FALSE, stringsAsFactors = FALSE
    ))
  }
  data.frame(
    label     = label,
    backend   = backend_name,
    covariance = cov,
    time_s    = elapsed,
    n_iter    = m$optim_par$iterations,
    loglik    = m$loglik,
    converged = (m$optim_par$status == 3),
    stringsAsFactors = FALSE
  )
}

## ---- Helper: run both backends for a given (formula, data, cov) ----
compare_both <- function(formula, data, cov, label) {
  cat(sprintf("  %s [%s]...\n", label, cov))
  rbind(
    fit_timed(formula, data, cov, ctrl_newton, "newton", label),
    fit_timed(formula, data, cov, ctrl_nlopt,  "nlopt",  label)
  )
}

## ---- Data preparation ----
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- prepare_data(mollusk$Abundance,     mollusk$Covariate)

## ---- Run all comparisons ----
results <- list()

cat("=== trichoptera (n=49, p=17) ===\n")
for (cov in c("full", "diagonal", "spherical")) {
  results[[length(results)+1]] <- compare_both(Abundance ~ 1,                  tri, cov, "tri_nocov")
  results[[length(results)+1]] <- compare_both(Abundance ~ Wind + Temperature, tri, cov, "tri_cov")
}

cat("=== barents (n=89, p=30) ===\n")
for (cov in c("full", "diagonal", "spherical")) {
  results[[length(results)+1]] <- compare_both(Abundance ~ 1,                   barents, cov, "bar_nocov")
  results[[length(results)+1]] <- compare_both(Abundance ~ Depth + Temperature, barents, cov, "bar_cov")
}

cat("=== mollusk (n=163, p=32) ===\n")
for (cov in c("full", "diagonal", "spherical")) {
  results[[length(results)+1]] <- compare_both(Abundance ~ 1,             mol, cov, "mol_nocov")
  results[[length(results)+1]] <- compare_both(Abundance ~ site + season, mol, cov, "mol_cov")
}

cat("=== oaks (n=116, p=114) ===\n")
for (cov in c("full", "diagonal", "spherical")) {
  results[[length(results)+1]] <- compare_both(Abundance ~ 1    + offset(log(Offset)), oaks, cov, "oak_nocov")
  results[[length(results)+1]] <- compare_both(Abundance ~ tree + offset(log(Offset)), oaks, cov, "oak_cov")
}

cat("=== microcosm (n=880, p=259) ===\n")
for (cov in c("diagonal", "spherical")) {
  results[[length(results)+1]] <- compare_both(Abundance ~ 1    + offset(log(Offset)), microcosm, cov, "mic_nocov")
  results[[length(results)+1]] <- compare_both(Abundance ~ site + offset(log(Offset)), microcosm, cov, "mic_cov")
}
cat("  microcosm full (slow)...\n")
for (lbl in c("mic_nocov", "mic_cov")) {
  form <- if (lbl == "mic_nocov") Abundance ~ 1 + offset(log(Offset)) else Abundance ~ site + offset(log(Offset))
  results[[length(results)+1]] <- compare_both(form, microcosm, "full", lbl)
}

cat("=== scRNA (n=3918, p=500) ===\n")
for (cov in c("diagonal", "spherical")) {
  results[[length(results)+1]] <- compare_both(counts ~ 1         + offset(log(total_counts)), scRNA, cov, "scr_nocov")
  results[[length(results)+1]] <- compare_both(counts ~ cell_line + offset(log(total_counts)), scRNA, cov, "scr_cov")
}
cat("  scRNA full covariance (very slow)...\n")
for (lbl in c("scr_nocov", "scr_cov")) {
  form <- if (lbl == "scr_nocov") counts ~ 1 + offset(log(total_counts)) else counts ~ cell_line + offset(log(total_counts))
  results[[length(results)+1]] <- compare_both(form, scRNA, "full", lbl)
}

cat("All fits done.\n\n")

## ---- Combine results ----
df <- do.call(rbind, results)
df$dataset <- sub("_.*", "", df$label)
df$covariates <- sub(".*_", "", df$label)

## ---- Summary table (wide format) ----
cat("========== COMPARISON SUMMARY ==========\n")
wide <- df %>%
  select(label, covariance, backend, time_s, n_iter, loglik, converged) %>%
  tidyr::pivot_wider(
    names_from  = backend,
    values_from = c(time_s, n_iter, loglik, converged)
  ) %>%
  mutate(
    loglik_diff  = loglik_newton - loglik_nlopt,
    speedup      = time_s_nlopt / time_s_newton
  ) %>%
  arrange(label, covariance)

print(wide %>% select(label, covariance,
                       time_newton = time_s_newton, time_nlopt = time_s_nlopt, speedup,
                       iter_newton = n_iter_newton,  iter_nlopt = n_iter_nlopt,
                       ll_newton   = loglik_newton,   ll_nlopt   = loglik_nlopt,
                       ll_diff     = loglik_diff,
                       conv_newton = converged_newton, conv_nlopt = converged_nlopt
                      ) %>%
        mutate(across(where(is.numeric), ~ signif(., 4))),
      row.names = FALSE, width = 200)

write.csv(wide, "inst/benchmark/backend_comparison.csv", row.names = FALSE)
cat("Saved: backend_comparison.csv\n")

## ---- Plot 1: time comparison ----
p1 <- ggplot(df, aes(x = paste(label, covariance, sep="\n"), y = time_s, fill = backend)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ dataset, scales = "free", nrow = 2) +
  scale_fill_manual(values = c(newton = "#E69F00", nlopt = "#56B4E9")) +
  labs(title = "Computation time: Newton vs nlopt",
       x = NULL, y = "Elapsed time (s)", fill = "Backend") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave("inst/benchmark/backend_time.pdf", p1, width = 18, height = 10)
cat("\nSaved: backend_time.pdf\n")

## ---- Plot 2: loglik difference (newton - nlopt) ----
df_wide <- df %>%
  pivot_wider(names_from = backend, values_from = c(time_s, n_iter, loglik, converged)) %>%
  mutate(ll_diff = loglik_newton - loglik_nlopt,
         fit = paste(label, covariance, sep=" / "))

p2 <- ggplot(df_wide, aes(x = fit, y = ll_diff,
                           fill = ifelse(ll_diff > 0, "Newton better", "nlopt better"))) +
  geom_col(width = 0.7) +
  facet_wrap(~ dataset, scales = "free", nrow = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Newton better" = "#009E73", "nlopt better" = "#D55E00"),
                    name = NULL) +
  labs(title = "loglik difference: Newton minus nlopt",
       subtitle = "Positive = Newton finds better solution",
       x = NULL, y = "loglik(Newton) - loglik(nlopt)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave("inst/benchmark/backend_loglik.pdf", p2, width = 18, height = 10)
cat("Saved: backend_loglik.pdf\n")

## ---- Plot 3: speedup (nlopt_time / newton_time) ----
p3 <- ggplot(df_wide, aes(x = fit, y = time_s_nlopt / time_s_newton,
                           fill = dataset)) +
  geom_col(width = 0.7) +
  facet_wrap(~ dataset, scales = "free", nrow = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(title = "Speedup: nlopt_time / newton_time",
       subtitle = "> 1 means Newton is faster",
       x = NULL, y = "Speedup ratio") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none")

ggsave("inst/benchmark/backend_speedup.pdf", p3, width = 18, height = 10)
cat("Saved: backend_speedup.pdf\n")

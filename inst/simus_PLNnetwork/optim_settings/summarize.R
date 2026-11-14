## Summary of sim_optim_settings.R: mean measures per setting, and paired differences
## with the legacy setting (former version), per block.
suppressMessages(library(dplyr))
files <- list.files(".", "^sim_.*\\.rds$")
r <- bind_rows(lapply(files, function(f) mutate(readRDS(f), block = sub("^sim_([a-z]+)_.*", "\\1", f))))
cat("errors:", sum(!is.na(r$error)), "\n")
r <- filter(r, is.na(error))
r %>% group_by(block, setting) %>%
  summarise(F1_size = mean(F1_size), maxF1 = mean(maxF1), AUPR = mean(AUPR), F1_BIC = mean(F1_BIC),
            F1_EBIC = mean(F1_EBIC), outer = mean(outer), time = mean(time), floored = mean(floored > 0), .groups = "drop") %>%
  print(width = 200, n = 100)
key <- c("block", "graph", "n", "p", "variance", "abundance", "contamination", "seed")
leg <- filter(r, setting == "legacy") %>% select(all_of(key), F1_leg = F1_size, AUPR_leg = AUPR, t_leg = time)
inner_join(filter(r, setting != "legacy"), leg, by = key) %>% group_by(block, setting) %>%
  summarise(dF1 = mean(F1_size - F1_leg), F1_better = mean(F1_size > F1_leg + 1e-9), F1_worse = mean(F1_size < F1_leg - 1e-9),
            dAUPR = mean(AUPR - AUPR_leg), time_ratio = mean(time) / mean(t_leg), .groups = "drop") %>%
  print(width = 200, n = 100)

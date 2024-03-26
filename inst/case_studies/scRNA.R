library(PLNmodels)
library(factoextra)

data(scRNA)
# data subsample: only 500 random cell and the 200 most varying transcript
scRNA        <- scRNA[sample.int(nrow(scRNA), 500), ]
scRNA$counts <- scRNA$counts[, 1:200]
myZIPLN <- ZIPLN(counts ~ 1 + offset(log(total_counts)), zi = "col", data = scRNA)
myPLN   <- PLN(counts ~ 1 + offset(log(total_counts)), data = scRNA)

true_zeros <- scRNA$counts == 0
n_zeros <- sum(true_zeros)
n_nzeros <- sum(!true_zeros)

data.frame(
    fitted   = c(as.vector(fitted(myPLN)[true_zeros]), as.vector(fitted(myZIPLN)[true_zeros])),
    obs      = rep(1:n_zeros, 2),
    method   = factor(rep(c("PLN", "ZI (1 par per species)"), each = n_zeros))
  ) %>%
  ggplot() +
  geom_rug(aes(x = fitted, color = method)) +
  geom_point(aes(x = obs, y = fitted, color = method), size = 0.3, alpha = 0.1) +
  ggtitle("Distribution of observed zeros as fitted by the models") +
  xlab(label = "observations #") +
  scale_color_viridis_d() + theme_bw()

# distribution of non zero
data.frame(
  fitted   = c(as.vector(fitted(myPLN)[!true_zeros]), as.vector(fitted(myZIPLN)[!true_zeros])),
  observed = rep(c(scRNA$counts[!true_zeros]), 2),
  method   = factor(rep(c("PLN", "ZI (1 par per species)"), each = n_nzeros))
) %>%
  ggplot(aes(x = observed, y = fitted)) +
  geom_point(size = .5, alpha =.25 ) +
  facet_wrap( ~ method) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() + annotation_logticks()


data.frame(
  fitted   = c(as.vector(fitted(myZIPLN)), as.vector(fitted(myPLN))),
  observed = rep(as.vector(scRNA$counts), 2),
  method   = factor(rep(c("ZIPLN", "PLN"), each = length(scRNA$counts)))
) %>%
  ggplot(aes(x = observed, y = fitted)) +
    geom_point(size = .5, alpha =.25 ) +
    facet_wrap( ~ method) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() + annotation_logticks()

prcomp(myZIPLN$latent) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = scRNA$cell_line)
prcomp(myZIPLN$latent_pos) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = scRNA$cell_line)

prcomp(myPLN$latent) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = scRNA$cell_line)
prcomp(myPLN$latent_pos) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = scRNA$cell_line)

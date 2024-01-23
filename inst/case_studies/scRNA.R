library(PLNmodels)
library(factoextra)

data(scRNA)
# data subsample: only 500 random cell and the 200 most varying transcript
scRNA        <- scRNA[sample.int(nrow(scRNA), 500), ]
scRNA$counts <- scRNA$counts[, 1:200]
myZIPLN <- ZIPLN(counts ~ 1 + offset(log(total_counts)), data = scRNA)
myPLN   <- PLN(counts ~ 1 + offset(log(total_counts)), data = scRNA)

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

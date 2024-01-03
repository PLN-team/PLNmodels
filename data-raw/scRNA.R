## code to prepare `scRNA` dataset goes here

library(scran)    # bioManager::install("scran")
library(scater)   # bioManager::install("scater")
library(ggplot2)  # install.packages("ggplot2")

## get data from "https://github.com/LuyiTian/sc_mixology/raw/master/data/sincell_with_class_5cl.RData")
load("~/Downloads/sincell_with_class_5cl.RData")


sce_qc <- computeSumFactors(sce_sc_10x_5cl_qc)
sce_log_qc <- logNormCounts(sce_sc_10x_5cl_qc)

sce_log_qc$cell_line <- as.factor(sce_log_qc$cell_line)
myPCA <- scater::calculatePCA(sce_log_qc)
top_features <- rownames(attr(myPCA, "rotation"))
counts <- data.matrix(t(counts(sce_log_qc)[top_features , ]))
scRNA <- data.frame(
  counts       = NA, ## placeholder for Abundance, to avoid using I() and inheriting "AsIs" class
  cell_line    = as.factor(sce_log_qc$cell_line),
  total_counts = sce_log_qc$total_count_per_cell
)
scRNA$counts <- counts

usethis::use_data(scRNA, overwrite = TRUE)

## code to prepare `oaks` dataset goes here

counts <- as(read.table(file = "data-raw/oaks/counts.tsv"), "matrix")
metadata <- read.table(file = "data-raw/oaks/metadata.tsv")
offsets <- as(read.table(file = "data-raw/oaks/offsets.tsv"), "matrix")

metadata$tree <- factor(metadata$tree, levels = c("susceptible", "intermediate", "resistant"))

oaks <- PLNmodels::prepare_data(counts     = counts,
                                     covariates = metadata,
                                     offset     = offsets)

usethis::use_data(oaks, overwrite = TRUE)



context("test-import-utils")
data("trichoptera")
counts     <- data.matrix(trichoptera$Abundance)
covariates <- trichoptera$Covariate

## Test common_samples ------------------------------------------------------------------------------------

test_that("common_samples throws warnings on matrices with no dimension names",  {
  ## No names for abundance matrix
  expect_warning(common_samples(`dimnames<-`(counts, NULL),
                                covariates))
  ## No rownames for abundance matrix (and therefore no matching names)
  expect_warning(common_samples(`rownames<-`(counts, NULL),
                                covariates))
  ## No names for covariates matrix
  ## (must transform covariates to matrix as data.frames always have rownames)
  expect_warning(common_samples(counts,
                                data.matrix(covariates, rownames.force = FALSE)))
})

test_that("common_samples fails on matrices with no dimension names and incompatibles dimensions",  {
  ## No names for abundance matrix
  expect_error(suppressWarnings(common_samples(unname(counts) %>% t(), covariates)))
  ## No rownames for covariates matrix
  expect_error(suppressWarnings(common_samples(counts %>% t(), data.matrix(covariates, rownames.force = FALSE))))
})

test_that("common_samples succeeds on matrices with no or partial dimension names but compatible dimensions",  {
  ## No dimensions names
  expect_warning(result <- common_samples(unname(counts), data.matrix(covariates)))
  expect_equal(result,
               list(transpose_counts = FALSE,
                    common_samples   = paste0("Sample_", 1:49),
                    default_names    = TRUE))
  ## Partial dimensions names: sample names only in count table
  expect_warning(result <- common_samples(`colnames<-`(counts, NULL),
                                          data.matrix(covariates)))
  expect_equal(result,
               list(transpose_counts = FALSE,
                    common_samples   = paste0(1:49),
                    default_names    = TRUE))
  ## Partial dimensions names: samples names only in covariates data.frame
  expect_warning(result <- common_samples(`rownames<-`(counts, NULL),
                                          covariates))
  expect_equal(result,
               list(transpose_counts = FALSE,
                    common_samples   = paste0(1:49),
                    default_names    = TRUE))
})

test_that("common_samples fails on matrices with dimension names but no common samples",  {
  expect_error(
    suppressWarnings(
      common_samples(
        `rownames<-`(counts, paste0("Sample_", 1:49)),
        covariates)),
    "Conflicting samples names in count matrix and covariates data frames"
  )
})

test_that("common_samples succeeds on matrices with dimension names and different orientations",  {
  expect_equal(common_samples(t(counts), covariates),
               list(transpose_counts = TRUE,
                    common_samples   = as.character(1:49),
                    default_names    = FALSE))
})

test_that("common_samples find biggest subset of common samples and produces message.",  {
  expect_message(result <- common_samples(counts, covariates[1:35, ]),
                 "14 samples were dropped from the abundance matrix for lack of associated covariates.")
  expect_length(result$common_samples, 35)
  expect_equal(result$common_samples,
               as.character(1:35))
})

## Test offset_functions --------------------------------------------------------------------------------

## all samples are identical up to a proportionality coefficient (induce explicit simple computations)
test_that("compute_offset provides correct answers for proportional samples", {
  sizes <- c(1, 2, 3, 5, 6)
  counts <- sizes %o% 1:10
  median_scale_size <- median(sizes)
  geom_mean_size <- exp(mean(log(sizes)))
  gmpr <- sapply(seq_along(sizes), function(i) { exp(mean(log(sizes[i]/sizes[-i]))) } )

  expect_equal(compute_offset(counts, "TSS"),  sizes * sum(1:10))
  expect_equal(compute_offset(counts, "CSS"),  sizes / median_scale_size)
  expect_equal(compute_offset(counts, "RLE"),  sizes / geom_mean_size)
  expect_equal(compute_offset(counts, "GMPR"), gmpr)
  expect_null(compute_offset(counts, "none"))
})

test_that("compute_offset provides correct answers for single row matrices", {
  sizes <- c(2)
  counts <- sizes %o% 1:10
  median_scale_size <- median(sizes)
  geom_mean_size <- exp(mean(log(sizes)))

  expect_equal(compute_offset(counts, "TSS"),  sizes * sum(1:10))
  expect_equal(compute_offset(counts, "CSS"),  sizes / median_scale_size)
  expect_equal(compute_offset(counts, "RLE"),  sizes / geom_mean_size)
  expect_error(compute_offset(counts, "GMPR"), "GMPR is not defined when there is only one sample.")
  expect_null(compute_offset(counts, "none"))
})

test_that("compute_offset provides correct answers for single column matrices", {
  sizes <- c(1, 2, 5, 6)
  counts <- sizes %o% 1
  median_scale_size <- median(sizes)
  geom_mean_size <- exp(mean(log(sizes)))
  gmpr <- sapply(seq_along(sizes), function(i) { exp(mean(log(sizes[i]/sizes[-i]))) } )

  expect_equal(compute_offset(counts, "TSS"),  sizes)
  expect_equal(compute_offset(counts, "CSS"),  sizes / median_scale_size)
  expect_equal(compute_offset(counts, "RLE"),  sizes / geom_mean_size)
  expect_equal(compute_offset(counts, "GMPR"), gmpr)
  expect_null(compute_offset(counts, "none"))
})

test_that("compute_offset provides correct answers for identical samples", {
  sizes <- rep(1, 5)
  counts <- sizes %o% 1:10

  expect_equal(compute_offset(counts, "TSS"),  sizes * sum(1:10))
  expect_equal(compute_offset(counts, "CSS"),  sizes)
  expect_equal(compute_offset(counts, "RLE"),  sizes)
  expect_equal(compute_offset(counts, "GMPR"), sizes)
  expect_null(compute_offset(counts, "none"))
})

test_that("offset_rle provides correct answers when adding pseudocounts", {
  sizes <- c(1, 2)
  counts <- sizes %o% rep(1, 10)
  geom_mean_size <- exp(mean(log(sizes+1)))
  ## External consistency
  expect_equal(compute_offset(counts, "RLE", pseudocount = 1),  (sizes+1) / geom_mean_size)
  ## Internal consistency
  expect_equal(compute_offset(counts, "RLE", pseudocount = 1),
               compute_offset(counts + 1, "RLE"))
})

## low or no redundancy across samples
test_that("offset_css throws a warning for nearly samples with limited inner variations", {
  sizes <- seq(from = 0.8, to = 1.2, by = 0.1)
  counts <- sizes %o% rep(1, 10)
  expect_warning(compute_offset(counts, "CSS"),
                 "No instability detected in quantile distribution across samples, falling back to scaled TSS normalization.")
})

test_that("offset_css throws a warning when a sample has less than two positive counts", {
  #      [,1] [,2] [,3]
  # [1,]    1    1    0
  # [2,]    0    0    1
  counts <- matrix(c(1, 1, 0,
                     0, 0, 1),
                   nrow = 2, byrow = TRUE)
  expect_warning(compute_offset(counts, "CSS"),
                 "Some samples only have 1 positive values. Can't compute quantiles and fall back to TSS normalization")
})

test_that("offset_rle throws a warning when data is too sparse but samples share a common species", {
  ## One common species, 4 specific species
  # [,1] [,2] [,3] [,4] [,5]
  # [1,]    1    1    0    1    0
  # [2,]    1    0    1    0    1
  counts <- matrix(c(1, 1, 0, 1, 0,
                     1, 0, 1, 0, 1),
                   nrow = 2, byrow = TRUE)
  expect_warning(compute_offset(counts, "RLE"),
                 "Because of high sparsity, some samples have null or infinite offset.")
})

test_that("offset_rle fails when no species is shared across samples", {
  #      [,1] [,2] [,3]
  # [1,]    1    1    0
  # [2,]    0    0    1
  counts <- matrix(c(1, 1, 0,
                     0, 0, 1),
                   nrow = 2, byrow = TRUE)
  expect_error(compute_offset(counts, "RLE"), "Sample do not share any common species, RLE normalization failed.")
})

test_that("offset_gmpr fails when a sample shares no species with any other sample", {
  #      [,1] [,2] [,3]
  # [1,]    1    1    0
  # [2,]    0    0    1
  counts <- matrix(c(1, 1, 0,
                     0, 0, 1),
                   nrow = 2, byrow = TRUE)
  expect_error(compute_offset(counts, "GMPR"),
               "Some sample(s) do not share any species with other samples, GMPR normalization failed.",
               fixed = TRUE)
})

## Numeric offset
test_that("offset_numeric fails when the offsets are incompatible with the counts table", {
  #   a b
  # A 1 4
  # B 2 5
  # C 3 6
  counts <- structure(1:6, .Dim = 3:2, .Dimnames = list(c("A", "B", "C"),
                                                      c("a", "b")))
  offset <- counts
  ## No sample names
  expect_error(offset_numeric(counts, unname(offset[, 1])),
               "Rownames are used for sample matching.\nPlease specify them in the offset vector/matrix.",
               fixed = TRUE
  )
  ## No offset for some samples
  rownames(offset) <- c("A", "B", "c")
  expect_error(offset_numeric(counts, offset),
               "Sample(s) C from the count table lack an offset.\nConsider checking your offset (orientation, rownames).", fixed = TRUE
               )
  ##  Wrong number of features
  rownames(offset) <- c("A", "B", "C")
  expect_error(offset_numeric(counts, offset[, c(1, 1, 2)]),
               "There should be one offset per feature in the count table.\nYou have 2 features but 3 offsets.",
               fixed = TRUE
  )
})

test_that("offset_numeric works for vectors and column-matrices.", {
  #   a b
  # A 1 4
  # B 2 5
  # C 3 6
  counts <- structure(1:6, .Dim = 3:2, .Dimnames = list(c("A", "B", "C"),
                                                        c("a", "b")))
  offset <- c(A = 1, C = 3, B = 2)
  ## No difference between vectors and column matrices
  expect_equal(offset_numeric(counts, offset),
               offset_numeric(counts, as.matrix(offset)))
  ##  Correct dimensions
  expect_equal(dim(offset_numeric(counts, offset)),
               dim(counts))
  ## Values in the correct order
  expect_equal(rownames(offset_numeric(counts, offset)),
               rownames(counts))
  ## Correct values
  expect_equal(offset_numeric(counts, offset)[, 1],
               offset[rownames(counts)])
})

test_that("offset_numeric works for matrices.", {
  #   a b
  # A 1 4
  # B 2 5
  # C 3 6
  counts <- structure(1:6, .Dim = 3:2, .Dimnames = list(c("A", "B", "C"),
                                                        c("a", "b")))
  offset <- rbind(counts, counts)
  rownames(offset) <- LETTERS[6:1] ## too many offsets
  ## Correct values
  expect_equal(offset_numeric(counts, offset),
               offset[c("A", "B", "C"), ])
})

## Test prepare_data functions --------------------------------------------------------------------------

test_that("prepare_data drops samples with no positive counts", {
  counts[1:2, ] <- 0
  expect_warning(result <- prepare_data(counts, covariates, offset = "none"),
                 "Sample(s) 1 and 2 dropped for lack of positive counts.",
                 fixed = TRUE
                 )
  expect_equal(dim(result), c(nrow(counts)-2, ncol(covariates)+1))
  expect_equal(dim(result$Abundance), c(nrow(counts)-2, ncol(counts)))
})

test_that("prepare_data replace NA with 0s", {
  counts[1:2, 1] <- NA
  expect_warning(result <- prepare_data(counts, covariates, offset = "none"),
                 "NA values in count table replaced with 0.")
  expect_equal(dim(result), c(nrow(counts), ncol(covariates)+1))
  expect_equal(dim(result$Abundance), c(nrow(counts), ncol(counts)))
})

result <- data.frame(
  Abundance = NA,
  covariates
  )
result$Abundance <- counts

test_that("prepare_data succeeds on simple data", {
  expect_identical(prepare_data(counts, covariates, offset = "none"),
                   result)
})

test_that("prepare_data succeeds on simple data with missing names", {
  expect_warning(res <- prepare_data(`rownames<-`(counts, NULL),
                                     covariates, offset = "none"))
  expect_identical(res,
                   ## small hack to account for the fact that setting rownames via rownames<-
                   ## changes its mode to character
                   `rownames<-`(result, rownames(result))
  )
})

test_that("prepare_data succeeds on simple data with missing names and numeric offset", {
  expect_warning(res <- prepare_data(`rownames<-`(counts, NULL),
                                     covariates,
                                     offset = rowSums(counts)))
  expect_equal(dim(res$Offset), dim(counts))
  expect_equal(rownames(res$Offset), as.character(1:nrow(res)))
  expect_equal(colnames(res$Offset), colnames(counts))
})

test_that("prepare_data provides automatic variable names when missing from covariate", {
  res <- prepare_data(counts, unname(covariates))
  expect_identical(colnames(res),
                   c("Abundance", paste0("Variable", 1:ncol(covariates)), "Offset"))
})

test_that("prepare_data succeeds when specifying a numeric offset", {
  ## On simple data
  res <- prepare_data(counts, covariates, offset = "TSS")
  res$Offset <- matrix(rep(res$Offset, ncol(counts)),
                       ncol = ncol(counts),
                       dimnames = dimnames(counts))
  expect_identical(prepare_data(counts, covariates, offset = rowSums(counts)),
                   res)
  ## On strange data (transposed count matrix)
  expect_identical(prepare_data(t(counts), covariates, offset = rowSums(counts)),
                   res)
})

test_that("prepare_data works on 1 column abundance matrices", {
  expect_warning(toy_data <- prepare_data(
    counts     = matrix(c(1, 3, 1, 1), ncol = 1),
    covariates = data.frame(Cov = c("A", "B", "A", "A")),
    offset     = rep(0, 4)
  ))
  expect_equal(dim(toy_data), c(4, 3))
  expect_equal(dim(toy_data$Abundance), c(4, 1))
})

## Test prepare_data_* functions ------------------------------------------------------------------------

# test_that("prepare_data_from_biom fails when covariates data.frame is missing", {
#   biom <- biomformat::make_biom(data = t(counts), sample_metadata = NULL)
#   expect_error(prepare_data_from_biom(biom, offset = "none"))
# })
#
# test_that("prepare_data_from_biom succeeds on biom specified by file name.", {
#   biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
#   biom <- biomformat::read_biom(biom_file)
#   expect_equal(prepare_data_from_biom(biom_file),
#                prepare_data_from_biom(biom))
# })
#
# test_that("prepare_data_from_biom succeeds on proper biom objects (but converts all columns to character...)", {
#   biom <- biomformat::make_biom(data = t(counts),
#                                sample_metadata = lapply(covariates, as.character) %>% as.data.frame)
#   result_biom <- result
#   result_biom[-1] <- lapply(result_biom[-1], as.character)
#   rownames(result_biom) <- as.character(rownames(result_biom))
#   expect_equal(prepare_data_from_biom(biom, offset = "none"),
#                result_biom)
# })
#
# ## massage rownames to fit phyloseq sample naming conventions
# rownames(counts) <- rownames(covariates) <- paste0("Sample", 1:nrow(covariates))
#
# test_that("prepare_data_from_phyloseq fails when covariate data.frame is missing", {
#   physeq <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = FALSE))
#   expect_error(prepare_data_from_phyloseq(physeq, offset = "none"),
#                "physeq should be a phyloseq object.")
#   physeq <- phyloseq::phyloseq(
#     phyloseq::otu_table(counts, taxa_are_rows = FALSE),
#     phyloseq::tax_table(matrix(rep("", ncol(counts)),
#                                dimnames = list(colnames(counts), "Kingdom")))
#   )
#   expect_error(prepare_data_from_phyloseq(physeq, offset = "none"),
#                paste("No covariates detected in physeq Consider:",
#                      "- extracting count data from biom with as(otu_table(physeq), \"matrix\")",
#                      "- preparing a covariates data.frame",
#                      "- using prepare_data instead of prepare_data_from_phyloseq",
#                      sep = "\n"),
#                fixed = TRUE)
# })
#
# test_that("prepare_data_from_phyloseq succeeds on proper phyloseq class objects", {
#   rownames(counts) <- rownames(covariates) <- paste0("Sample", 1:nrow(covariates))
#   physeq <- phyloseq::phyloseq(
#     phyloseq::otu_table(counts, taxa_are_rows = FALSE),
#     phyloseq::sample_data(covariates)
#     )
#   result_physeq <- result; rownames(result_physeq$Abundance) <- rownames(result_physeq) <- paste0("Sample", 1:nrow(covariates))
#   expect_equal(prepare_data_from_phyloseq(physeq, offset = "none"), result_physeq)
# })

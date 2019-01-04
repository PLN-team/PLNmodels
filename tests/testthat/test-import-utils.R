context("test-import-utils")
data("trichoptera")

counts <- trichoptera$Abundance
# covariates <- trichoptera %>% inset2("Abundance", value = NULL)
covariates <- trichoptera; covariates$Abundance <- NULL

## Test common_samples ------------------------------------------------------------------------------------

test_that("common_samples throws warnings on matrices with no dimension names",  {
  ## No names for abundance matrix
  expect_warning(common_samples(`dimnames<-`(counts, NULL),
                                covariates))
  ## No names for covariates matrix
  ## (must transform covariates to matrix as data.frames always have rownames)
  expect_warning(common_samples(counts,
                                data.matrix(covariates, rownames.force = FALSE)))
})

test_that("common_samples fails on matrices with no dimension names and incompatibles dimensions",  {
  ## No names for abundance matrix
  expect_error(common_samples(`dimnames<-`(counts, NULL) %>% t(),
                                covariates))
  ## No rownames for covariates matrix
  expect_error(common_samples(counts %>% t(),
                              data.matrix(covariates, rownames.force = FALSE)))
})

test_that("common_samples succeeds on matrices with no dimension names but compatible dimensions",  {
  expect_length(common_samples(`dimnames<-`(counts, NULL),
                               data.matrix(covariates)),
                2)
  expect_named(common_samples(`dimnames<-`(counts, NULL),
                               data.matrix(covariates)),
               c("transpose_counts", "common_samples"))
  expect_type(common_samples(`dimnames<-`(counts, NULL),
                             data.matrix(covariates)),
              "list")
  expect_type(common_samples(`dimnames<-`(counts, NULL),
                             data.matrix(covariates))$transpose_counts,
              "logical")
  expect_type(common_samples(`dimnames<-`(counts, NULL),
                             data.matrix(covariates))$common_samples,
              "character")
  expect_equal(common_samples(`dimnames<-`(counts, NULL),
                              data.matrix(covariates)),
               list(transpose_counts = FALSE,
                    common_samples   = paste0("Sample_", 1:49)))
})

test_that("common_samples fails on matrices with dimension names but no common samples",  {
  expect_error(common_samples(`rownames<-`(counts, paste0("Sample_", 1:49)),
                               covariates))
})

test_that("common_samples succeeds on matrices with dimension names and different orientations",  {
  expect_length(common_samples(t(counts),covariates),
                2)
  expect_named(common_samples(t(counts), covariates),
               c("transpose_counts", "common_samples"))
  expect_type(common_samples(t(counts), covariates),
              "list")
  expect_type(common_samples(t(counts), covariates)$transpose_counts,
              "logical")
  expect_type(common_samples(t(counts), covariates)$common_samples,
              "character")
  expect_equal(common_samples(t(counts), covariates),
               list(transpose_counts = TRUE,
                    common_samples   = as.character(1:49)))
})

test_that("common_samples find biggest subset of common samples and produces message.",  {
  expect_message(common_samples(counts, covariates[1:35, ]),
                 "14 samples were dropped from the abundance matrix for lack of associated covariates.")
  expect_length(common_samples(counts, covariates[1:35, ])$common_samples, 35)
  expect_equal(common_samples(counts, covariates[1:35, ])$common_samples,
               as.character(1:35))
})

## Test offset_functions --------------------------------------------------------------------------------

## all samples are identical up to a proportionality coefficient (induce explicit simple computations)
test_that("compute_offset provides correct answers for proportional samples", {
  sizes <- c(1, 2, 5, 6)
  counts <- sizes %o% 1:10
  median_scale_size <- median(sizes)
  geom_mean_size <- exp(mean(log(sizes)))

  expect_equal(compute_offset(counts, "TSS"),  sizes * sum(1:10))
  expect_equal(compute_offset(counts, "CSS"),  sizes / median_scale_size)
  expect_equal(compute_offset(counts, "RLE"),  sizes / geom_mean_size)
  expect_equal(compute_offset(counts, "GMPR"), sizes / geom_mean_size)
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
  expect_equal(compute_offset(counts, "GMPR"), sizes / geom_mean_size)
  expect_null(compute_offset(counts, "none"))
})

test_that("compute_offset provides correct answers for single column matrices", {
  sizes <- c(1, 2, 5, 6)
  counts <- sizes %o% 1
  median_scale_size <- median(sizes)
  geom_mean_size <- exp(mean(log(sizes)))

  expect_equal(compute_offset(counts, "TSS"),  sizes)
  expect_equal(compute_offset(counts, "CSS"),  sizes / median_scale_size)
  expect_equal(compute_offset(counts, "RLE"),  sizes / geom_mean_size)
  expect_equal(compute_offset(counts, "GMPR"), sizes / geom_mean_size)
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

## low or no redundancy across samples
test_that("offset_css throws a warning for nearly samples with limited inner variations", {
  sizes <- seq(from = 0.8, to = 1.2, by = 0.1)
  counts <- sizes %o% rep(1, 10)
  expect_warning(compute_offset(counts, "CSS"),
                 "No instability detected in quantile distribution across samples, falling back to scaled TSS normalization.")
})

test_that("offset_css throws a warning has less than two positive counts", {
  #      [,1] [,2] [,3]
  # [1,]    1    1    0
  # [2,]    0    0    1
  counts <- matrix(c(1, 1, 0,
                     0, 0, 1),
                   nrow = 2, byrow = TRUE)
  expect_warning(compute_offset(counts, "CSS"),
                 "Some samples only have 1 positive values. Can't compute quantiles and fall back to TSS normalization")
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

test_that("prepare_data succeeds on simple data", {
  expect_identical(prepare_data(counts, covariates, offset = "none"),
                   trichoptera)
})

## Test prepare_data_* functions ------------------------------------------------------------------------



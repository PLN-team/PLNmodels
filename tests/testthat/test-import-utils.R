context("test-import-utils")
data("trichoptera")

counts <- trichoptera$Abundance
## covariates <- `$<-`(trichoptera, "Abundance", NULL)
covariates <- trichoptera; covariates$Abundance <- NULL

## Test common_samples

test_that("common_samples throws warnings on matrices with no dimension names",  {
  ## No names for abundance matrix
  expect_warning(common_samples(`dimnames<-`(counts, NULL),
                                covariates))
  ## No names for covariates matrix
  ## (must transform covariates to matrix as data.frames always have rownames)
  expect_warning(common_samples(counts,
                                data.matrix(covariates)))
})

test_that("common_samples fails on matrices with no dimension names and incompatibles dimensions",  {
  ## No names for abundance matrix
  expect_error(common_samples(`dimnames<-`(counts, NULL) %>% t(),
                                covariates))
  ## No rownames for covariates matrix
  expect_error(common_samples(counts %>% t(),
                              data.matrix(covariates)))
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
                    common_samples   = str_c("Sample_", 1:49)))
})

test_that("common_samples fails on matrices with dimension names but no common samples",  {
  expect_error(common_samples(`rownames<-`(counts, str_c("Sample_", 1:49)),
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

## Test prepare_data

## Test offset_functions

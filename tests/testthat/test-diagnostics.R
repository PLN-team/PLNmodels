###############################################################################
## Diagnostics of structural zeros: species absent from a whole group of
## samples, which a PLN model without the right covariate or zero-inflation
## fits by sending latent means to minus infinity. They are looked for in the
## data (structural_zeros(), run by prepare_data()), and reported after a fit
## through the latent variances (field degenerate_species).
###############################################################################

test_that("structural_zeros finds the species absent from a level, and only them", {
  set.seed(1)
  n <- 60
  group <- factor(rep(c("a", "b", "c"), each = n / 3))
  Y <- matrix(rpois(n * 5, 5), n, 5, dimnames = list(NULL, paste0("sp", 1:5)))
  Y[group == "a", "sp2"] <- 0                    # absent from a whole level
  Y[group != "c", "sp3"] <- 0                    # absent from two levels
  Y[, "sp4"] <- 0; Y[c(21, 45), "sp4"] <- 3      # rare: absent from level a by chance
  Y[, "sp5"] <- 0                                # never observed: nothing to report
  covariates <- data.frame(group = group, x = rnorm(n), id = paste0("s", 1:n))

  z <- structural_zeros(Y, covariates)
  expect_named(z, c("species", "covariate", "level", "n_samples", "prevalence_elsewhere", "p_value"))
  expect_setequal(paste(z$species, z$level), c("sp2 a", "sp3 a", "sp3 b"))
  expect_true(all(z$covariate == "group"))
  expect_true(all(z$n_samples == 20))
  expect_true(all(z$p_value <= 0.05))
  expect_false(is.unsorted(z$p_value))
  ## sp2 is in all the samples of the other levels, sp3 in half of them
  expect_equal(z$prevalence_elsewhere[z$species == "sp2"], 1)
  expect_equal(z$prevalence_elsewhere[z$species == "sp3"], c(0.5, 0.5))

  ## numeric covariates are not used, and nothing to report gives no row
  expect_equal(nrow(structural_zeros(Y, covariates["x"])), 0)
  expect_equal(nrow(structural_zeros(Y[, c("sp1", "sp4", "sp5")], covariates)), 0)
  ## missing values in the factor are left out
  covariates$group[1:3] <- NA
  expect_setequal(paste(structural_zeros(Y, covariates)$species), c("sp2", "sp3", "sp3"))
})

test_that("structural_zeros finds the oak species that live on a single tree", {
  data(oaks)
  z <- structural_zeros(oaks$Abundance, oaks[c("tree", "orientation")])
  expect_true(all(z$covariate == "tree"))
  expect_true(all(c("f_OTU_1011", "f_OTU_46", "f_OTU_1090") %in% z$species))
  ## f_OTU_1011 is on all the intermediate trees and on no other
  expect_setequal(z$level[z$species == "f_OTU_1011"], c("susceptible", "resistant"))
})

test_that("prepare_data reports structural zeros in a message, and leaves the data as is", {
  data(mollusk)
  expect_message(
    mol <- suppressWarnings(prepare_data(mollusk$Abundance, mollusk$Covariate)),
    "absent from all the samples of a level"
  )
  expect_equal(ncol(mol$Abundance), ncol(mollusk$Abundance))

  data(trichoptera)
  expect_no_message(prepare_data(trichoptera$Abundance, trichoptera$Covariate))
})

test_that("species with a degenerate latent variance are reported, once per call", {
  data(trichoptera)
  tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

  ## nothing to report at the default threshold
  expect_no_warning(fit <- PLN(Abundance ~ 1, tri, control = PLN_param(trace = 0)))
  expect_length(fit$degenerate_species, 0)

  ## a low threshold stands for a degenerate fit
  old <- options(PLNmodels.latent_variance_threshold = 5)
  on.exit(options(old))
  variances <- diag(sigma(fit))
  expect_setequal(fit$degenerate_species, names(variances)[variances > 5])
  expect_gt(length(fit$degenerate_species), 0)
  expect_warning(PLN(Abundance ~ 1, tri, control = PLN_param(trace = 0)), "latent variance")
  expect_warning(ZIPLN(Abundance ~ 1, tri, control = ZIPLN_param(trace = 0)), "latent variance")

  ## a single warning for a whole collection
  warns <- capture_warnings(
    nets <- PLNnetwork(Abundance ~ 1, tri, control = PLNnetwork_param(trace = 0, n_penalties = 5))
  )
  expect_length(grep("latent variance", warns), 1)
  expect_match(warns[grep("latent variance", warns)], "of the 5 models")
  expect_warning(ZIPLNnetwork(Abundance ~ 1, tri, control = ZIPLNnetwork_param(trace = 0, n_penalties = 3)),
                 "of the 3 models")
})

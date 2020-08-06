## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  INTERNAL FUNCTIONS ---------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Internal function to find the most comprehensive set of common samples between a count table and a covariates data.frame
common_samples <- function(counts, covariates) {
  ## Have default samples names been created when matching samples?
  default_names <- FALSE
  ## Sanity checks:
  name_warning <- paste("There are no matching names in the count matrix and the covariates data.frame.",
                        "Function will proceed assuming:",
                        "- samples are in the same order;",
                        "- samples are rows of the abundance matrix.", sep = "\n")
  ## no sample names in covariates: create sample names
  if (is.null(rownames(covariates))) {
    warning(name_warning)
    if (nrow(counts) != nrow(covariates)) {
      stop("Incompatible dimensions")
    }
    if (is.null(rownames(counts))) rownames(counts) <- paste0("Sample_", 1:nrow(counts))
    rownames(covariates) <- rownames(counts)
    default_names <- TRUE
  }
  ## Attempt name matching between covariates and count data
  count_names <- unlist(dimnames(counts))
  if (is.null(count_names) || !any(count_names %in% rownames(covariates))) {
    # If name matching is impossible, abort if
    ## - dimension are incompatible or
    ## - row (samples) names are conflicting
    warning(name_warning)
    if (nrow(counts) != nrow(covariates)) {
      stop("Incompatible dimensions")
    }
    if (!is.null(rownames(counts))) {
      stop("Conflicting samples names in count matrix and covariates data frames")
    }
    rownames(counts) <- rownames(covariates)
    default_names <- TRUE
  }
  ## Check whether samples are stored as columns in the abundance matrix
  ## and transpose if that's the case
  ## Based on a heuristic of matching names
  sample_are_cols <- any(colnames(counts) %in% rownames(covariates))
  if (sample_are_cols) counts <- t(counts)
  ## Ensure consistency by using only common samples
  common_samples <- intersect(rownames(counts), rownames(covariates))
  if (length(common_samples) < nrow(counts)) {
    message(paste0(nrow(counts) - length(common_samples), " samples were dropped from the abundance matrix for lack of associated covariates."))
  }
  return(list(transpose_counts = sample_are_cols,
              common_samples   = common_samples,
              default_names    = default_names))
}

## scaling functions --------

## Sanitize offset to ensure consistency with count matrix
sanitize_offset <- function(counts, offset, ...) {
  p <- ncol(counts) ## number of features
  ## Sanity check: transform vector offset and column-matrices to full matrices
  if (is.vector(offset) || (is.matrix(offset) && ncol(offset) == 1)) {
    offset_samples <- if (is.vector(offset)) {
      names(offset)
    } else {
      rownames(offset)
    }
    offset <- matrix(rep(offset, p),
                     ncol = p,
                     dimnames = list(offset_samples, colnames(counts)))
  }
  ## Sanity check: rownames
  if (is.null(rownames(offset))) {
    stop("Rownames are used for sample matching.\nPlease specify them in the offset vector/matrix.")
  }
  ## Sanity checks: offsets are available for all samples
  if (anyNA(ii <- match(rownames(counts), rownames(offset)))) {
    stop(paste("Sample(s) "),
         paste(rownames(counts)[is.na(ii)], collapse = " and "),
         " from the count table lack an offset.\nConsider checking your offset (orientation, rownames).")
  }
  offset <- offset[rownames(counts), , drop = FALSE]
  ## Sanity checks: offset are available for all species
  if (ncol(offset) != p) {
    stop(paste("There should be one offset per feature in the count table.\nYou have",
               p,
               "features but",
               ncol(offset),
               "offsets."))
  }
  offset
}



## Numeric offset
offset_numeric <- function(counts, offset, ...) {
  sanitize_offset(counts, offset, ...)
}

## No offset
offset_none <- function(counts) {
  return(NULL)
}

## Total Sum Scaling offset
offset_tss <- function(counts) {
  rowSums(counts)
}

## Geometric Mean Pairwise Ratio (GMPR) normalization (as presented in doi.org/10.7717/peerj.4600)
offset_gmpr <- function(counts) {
  if (nrow(counts) == 1) stop("GMPR is not defined when there is only one sample.")
  ## median of (non-null, non-infinite) pairwise ratios between counts of samples i and j
  pairwise_ratio <- function(i, j) { median(counts[i, ] / counts[j, ], na.rm = TRUE) }
  ## Matrix of pairwise ratios
  n <- nrow(counts)
  mat_pr <- matrix(NaN, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) mat_pr[i, j] <- pairwise_ratio(i, j)
    }
  }
  ## Geometric mean of pairwise ratio
  size_factor <- apply(mat_pr, 1, geom_mean)
  if (any(size_factor == 0 | !is.finite(size_factor))) stop("Some sample(s) do not share any species with other samples, GMPR normalization failed.")
  return(size_factor)
}

## Relative Log Expression (RLE) normalization (as used in DESeq2)
offset_rle <- function(counts, pseudocounts = 0L, type = c("ratio", "poscounts")) {
  type <- match.arg(type)
  ## Count manipulation: pseudo and replace 0s with NA (ignored in geometric mean computations)
  counts <- counts + pseudocounts
  if (type == "poscounts") counts[counts == 0] <- NA
  ## compute simple geometric mean for all otus
  geom_means <- counts %>% log() %>% colMeans(na.rm = TRUE) %>% exp
  if (all(geom_means == 0 | is.nan(geom_means))) stop("Samples do not share any common species, RLE normalization failed.")
  ## compute size factor as the median of all otus log-ratios in that sample
  robust_med <- function(cnts) { median((cnts/geom_means)[is.finite(geom_means) & cnts > 0], na.rm = TRUE) }
  size_factor <- apply(counts, 1, robust_med)
  if (any(is.infinite(size_factor) | size_factor == 0)) warning("Because of high sparsity, some samples have null or infinite offset.")
  return(size_factor)
}

## Cumulative Sum Scaling (CSS) normalization (as used in metagenomeSeq and presented in doi.org/10.1038/nmeth.2658)
offset_css <- function(counts, reference = median) {
  ## special treatment for edge case of one-column matrix (1 OTU, many samples)
  if (ncol(counts) == 1) return( counts[, 1] / median(counts) )
  ## remove 0s and check that all samples have at least two positive counts
  counts[counts == 0] <- NA
  if (any(rowSums(!is.na(counts)) < 2)) {
    warning("Some samples only have 1 positive values. Can't compute quantiles and fall back to TSS normalization")
    return(rowSums(counts, na.rm = TRUE))
  }
  ## compute sample-specific quantiles and cumulative sums up to quantiles
  cumsum_up_to <- function(counts, quantiles) {
    (counts * outer(counts, quantiles, `<=`)) %>% colSums(na.rm = TRUE)
  }
  mat_sample_quant <- apply(counts, 1, quantile, probs = seq(0, 1, length.out = ncol(counts)), na.rm = TRUE) %>% t()
  mat_sample_cumsum <- sapply(1:nrow(counts), function(i) { cumsum_up_to(counts[i, ], mat_sample_quant[i, ]) }) %>% t()
  ## reference quantiles, computed as median (nature article) or mean (metagenomeSeq::cumNormStat[Fast]) of sample_specific quantiles
  ## and MAD around the reference quantiles
  ref_quant <- apply(mat_sample_quant, 2, reference)
  ref_quant_mad <- sweep(mat_sample_quant, 2, ref_quant) %>% abs %>% apply(2, median)
  ## find smallest quantile for which high instability is detected
  ## instability for quantile l is defined as ref_quant_mad[l+1] - ref_quant_mad[l] >= 0.1 * ref_quant_mad[l]
  instable <- (diff(ref_quant_mad) >= 0.1 * head(ref_quant_mad, -1))
  if (any(instable)) {
    ## Hack to mimick package implementation: never choose quantile below 50%
    lhat <- max(min(which(instable)), ceiling(ncol(counts)/2))
  } else {
    warning("No instability detected in quantile distribution across samples, falling back to scaled TSS normalization.")
    lhat <- ncol(counts)
  }
  ## scaling factors are cumulative sums up to quantile lhat, divided by their median
  size_factors <- mat_sample_cumsum[ , lhat] / median(mat_sample_cumsum[ , lhat])
  return(size_factors %>% unname())
}

## Wrench normalization (from doi:10.1186/s12864-018-5160-5) in its simplest form
## @importFrom matrixStats rowWeightedMeans rowVars
offset_wrench <- function(counts, groups = rep(1, nrow(counts))) {

  ## n = number of samples, p = number of species
  n <- nrow(counts); p <- ncol(counts)
  if (n == 1) stop("Wrench is not defined when there is only one sample.")

  ## Proportions
  depths <- rowSums(counts)
  props  <- counts / depths
  ## Proportions in reference
  props_ref <- colMeans(props)
  ## Samplewise ratios of proportions
  sample_ratios <- props / matrix(props_ref, n, p, byrow = TRUE)

  ## Sample specific variance
  species_var <- species_variance(counts, groups)

  ## Correction by group effect (if needed)
  g = length(unique(groups)) ## number of groups
  if (g > 1) {
    groups <- as.character(groups)
    ## Groupwise proportions and ratios
    group_counts <- rowsum(counts, groups, reorder = TRUE)
    ## group_depths <- rowSums(group_counts)
    group_props <- group_counts / rowSums(group_counts)
    ## log groupwise ratios of proportions (with 0s replaced with NA)
    group_ratio <- group_props / matrix(props_ref, g, p, byrow = TRUE)

    ## groupwise scale and dispersion factors
    group_scales <- rowMeans(group_ratio)
    # log_group_ratio <- log(group_ratio)
    # log_group_ratio[!is.finite(x)] <- NA
    # group_log_var <- matrixStats::rowVars(log_group_ratio, na.rm = TRUE)
    log_var <- function(x) { var(log(x[is.finite(x) & x > 0])) }
    group_log_var <- apply(group_ratio, 1, log_var)

    ## regularized estimation of positive means.
    ## Correct by group effect
    sample_ratios[] <- sample_ratios / group_scales[groups]

    ## theta_gi: average ratio in sample i (from group g). Computed as the (geometric) mean ratio in the sample. The contributions of each ratio is shrinked inversely to (i) the species variance and (ii) the group variance
    ## weights[sample, species] = 1 / (sigma2[species] + sigma2[sample_group])
    weights <- 1 / outer(group_log_var[groups], species_var, "+")
    weights[] <- weights / rowSums(weights)
    ## Two way effect of sample i in group g
    theta_gi <- exp(rowSums(sample_ratios * weights, na.rm = TRUE))

    ## theta_gj: average deviation from the average ratio for species j (in group j), shrinked by the ratio between (i) group variance and (i) species variance plus group variance
    # shrinkage[sample, species] = sigma2[sample_group] / (sigma2[species] + sigma2[sample_group])
    shrinkage <- group_log_var[groups] / outer(group_log_var[groups], species_var, "+")
    theta_gj <- exp(shrinkage * (log(sample_ratios) - log(matrix(theta_gi, n, p))))

    ## Compute sample ratios from the previous shrinked quantities: group effect + sample effect within group + species effect within sample within group
    sample_ratios[] <- theta_gj * (theta_gi * group_scales[groups])
  }

  ## Compositional factors: robust means of ratio
  # comp_factors <- matrixStats::rowWeightedMeans(sample_ratios, 1 / species_var)
  weights <- 1 / matrix(species_var, nrow = n, ncol = p, byrow = TRUE)
  weights[] <- weights / rowSums(weights)
  comp_factors <- rowSums(sample_ratios * weights)
  comp_factors[] <- comp_factors / geom_mean(comp_factors) ## Centered in geometric scale

  ## Normalization factors
  norm_factors <- comp_factors * depths / geom_mean(depths)

  return(norm_factors)
}

# Helpers scaling functions ----

## Geometric mean (computed only on positive values)
geom_mean <- function(x, poscounts = TRUE, na.rm = TRUE) {
  x_log <- log(x)
  if (poscounts) x_log <- x_log[x > 0]
  exp(mean(x_log, na.rm = na.rm))
}

species_variance <- function(counts, groups = rep(1, nrow(counts))) {
  ## n = number of samples, p = number of species
  n <- nrow(counts); p <- ncol(counts)

  ## Centered log depths and counts corrected by offset
  log_depths <- log(rowSums(counts))
  log_depths[] <- log_depths
  log_counts <- log(counts) - log_depths
  log_counts[!is.finite(log_counts)] <- NA

  ## Design matrix
  groups <- as.character(groups)
  if (length(unique(groups)) > 1) {
    design <- model.matrix(counts[, 1] ~ 0 + groups)
  } else {
    design <- model.matrix(counts[, 1] ~ 1)
  }
  ## Manually regress log_counts against log_depths (for each species) to compute species-specific sigma.
  compute_sigma <- function(j) {
    valid_obs <- !is.na(log_counts[, j])
    model <- .lm.fit(design[valid_obs, , drop = FALSE], log_counts[valid_obs, j])
    sum(model$residuals^2) / (sum(valid_obs) - model$rank)
  }
  ## For numerical stability, set minimum variance to .Machine$double.eps
  pmax(vapply(seq(p), compute_sigma, numeric(1)),
       .Machine$double.eps)
}

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## EXPORTED FUNCTIONS ---------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Prepare data for use in PLN models
#' @name prepare_data
#'
#' @description Prepare data in proper format for use in PLN model and its variants. The function (i) merges a count table and
#' a covariate data frame in the most comprehensive way and (ii) computes offsets from the count table using one of several normalization schemes (TSS, CSS, RLE, GMPR, Wrench, etc). The function fails with informative messages when the heuristics used for sample matching fail.
#'
#' @param counts Required. An abundance count table, preferably with dimensions names and species as columns.
#' @param covariates Required. A covariates data frame, preferably with row names.
#' @param offset Optional. Normalization scheme used to compute scaling factors used as offset during PLN inference. Available schemes are "TSS" (Total Sum Scaling, default), "CSS" (Cumulative Sum Scaling, used in metagenomeSeq), "RLE" (Relative Log Expression, used in DESeq2), "GMPR" (Geometric Mean of Pairwise Ratio, introduced in Chen et al., 2018), Wrench (introduced in Kumar et al., 2018) or "none". Alternatively the user can supply its own vector or matrix of offsets (see note for specification of the user-supplied offsets).
#' @param ... Additional parameters passed on to [compute_offset()]
#'
#' @references Chen, L., Reeve, J., Zhang, L., Huang, S., Wang, X. and Chen, J. (2018) GMPR: A robust normalization method for zero-inflated count data with application to microbiome sequencing data. PeerJ, 6, e4600 \url{https://doi.org/10.7717/peerj.4600}
#' @references Paulson, J. N., Colin Stine, O., Bravo, H. C. and Pop, M. (2013) Differential abundance analysis for microbial marker-gene surveys. Nature Methods, 10, 1200-1202 \url{http://dx.doi.org/10.1038/nmeth.2658}
#' @references Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data. Genome Biology, 11, R106 \url{https://doi.org/10.1186/gb-2010-11-10-r106}
#' @references Kumar, M., Slud, E., Okrah, K. et al. (2018) Analysis and correction of compositional bias in sparse sequencing count data. BMC Genomics 19, 799 \url{https://doi.org/10.1186/s12864-018-5160-5}
#'
#' @return A data.frame suited for use in [PLN()] and its variants with two specials components: an abundance count matrix (in component "Abundance") and an offset vector/matrix (in component "Offset", only if offset is not set to "none")
#' @note User supplied offsets should be either vectors/column-matrices or have the same number of column as the original count matrix and either (i) dimension names or (ii) the same dimensions as the count matrix. Samples are trimmed in exactly the same way to remove empty samples.
#'
#'
#' @seealso [compute_offset()] for details on the different normalization schemes
#'
#' @export
#'
#' @examples
#' data(trichoptera)
#' proper_data <- prepare_data(
#'  counts     = trichoptera$Abundance,
#'  covariates = trichoptera$Covariate,
#'  offset     = "TSS"
#' )
#' proper_data$Abundance
#' proper_data$Offset
prepare_data <- function(counts, covariates, offset = "TSS", ...) {
  ## Convert counts and covariates to expected format
  counts     <- data.matrix(counts, rownames.force = TRUE)
  covariates <- as.data.frame(covariates)
  ## sanitize abundance matrix and covariates data.frame
  common <- common_samples(counts, covariates)
  samples <- common$common_samples
  ## sanitize offset
  if (is.numeric(offset) && is.vector(offset)) {
    offset <- matrix(offset, ncol = 1, dimnames = list(names(offset), NULL))
  }
  if (common$transpose_counts) counts <- t(counts)
  if (common$default_names) {
    rownames(counts) <- rownames(covariates) <- samples
    if (is.numeric(offset)) rownames(offset) <- samples
  }
  counts <- counts[samples, , drop = FALSE]
  ## Replace NA with 0s
  if (any(is.na(counts))) {
    counts[is.na(counts)] <- 0
    warning("NA values in count table replaced with 0.")
  }
  ## filter out empty samples
  empty_samples <- which(rowSums(counts) == 0)
  if (length(empty_samples)) {
    warning(paste0("Sample(s) ",
                   paste(samples[empty_samples], collapse = " and "),
                   " dropped for lack of positive counts."))
    samples <- samples[-empty_samples]
    counts <- counts[samples, ,drop = FALSE]
  }
  covariates <- covariates[samples, , drop = FALSE]
  if (is.null(names(covariates))) names(covariates) <- paste0("Variable", seq_along(covariates))
  ## compute offset
  offset     <- compute_offset(counts, offset, ...)
  ## prepare data for PLN
  result <- data.frame(Abundance = NA, ## placeholder for Abundance, to avoid using I() and inheriting "AsIs" class
                       covariates,
                       Offset    = NA ## placeholder for Offset, to avoid using I() and inheriting "AsIs" class
                       )
  result$Abundance <- counts
  result$Offset <- offset
  return(result)
}


#' @title Compute offsets from a count data using one of several normalization schemes
#' @name compute_offset
#'
#' @description Computes offsets from the count table using one of several normalization schemes (TSS, CSS, RLE, GMPR, etc) described in the literature.
#'
#' @inheritParams prepare_data
#' @param ... Additional parameters passed on to specific methods (for now CSS and RLE)
#' @inherit prepare_data references
#'
#' @details RLE has additional `pseudocounts` and `type` arguments to add pseudocounts to the observed counts (defaults to 0L) and to compute offsets using only positive counts (if `type == "poscounts"`). This mimicks the behavior of [DESeq2::DESeq()] when using `sfType == "poscounts"`. CSS has an additional `reference` argument to choose the location function used to compute the reference quantiles (defaults to `median` as in the Nature publication but can be set to `mean` to reproduce behavior of functions cumNormStat* from metagenomeSeq). Note that (i) CSS normalization fails when the median absolute deviation around quantiles does not become instable for high quantiles (limited count variations both within and across samples) and/or one sample has less than two positive counts, (ii) RLE fails when there are no common species across all samples (unless `type == "poscounts"` has been specified) and (iii) GMPR fails if a sample does not share any species with all other samples.
#'
#' @return If `offset = "none"`, `NULL` else a vector of length `nrow(counts)` with one offset per sample.
#'
#' @importFrom stats mad median quantile
#' @export
#'
#' @examples
#' data(trichoptera)
#' counts <- trichoptera$Abundance
#' compute_offset(counts)
#' ## Other normalization schemes
#' compute_offset(counts, offset = "GMPR")
#' compute_offset(counts, offset = "RLE", pseudocounts = 1)
#' ## User supplied offsets
#' my_offset <- setNames(rep(1, nrow(counts)), rownames(counts))
#' compute_offset(counts, offset = my_offset)
compute_offset <- function(counts, offset = c("TSS", "GMPR", "RLE", "CSS", "Wrench", "none"), ...) {
  ## special behavior for data.frame
  if (inherits(offset, "data.frame")) {
    stop(
  "You supplied a data.frame to compute_offset(). Did you mean to supply a numeric matrix?
  Try converting your data.frame to a matrix with as.matrix()."
  )
  }
  ## special behavior for numeric offset
  if (is.numeric(offset)) {
    return(offset_numeric(counts, offset, ...))
  }
  ## Choose offset function
  offset <- match.arg(offset)
  offset_function <- switch(offset,
                            "TSS"    = offset_tss,
                            "GMPR"   = offset_gmpr,
                            "RLE"    = offset_rle,
                            "CSS"    = offset_css,
                            "Wrench" = offset_wrench,
                            "none"   = offset_none
  )
  ## Ensure that counts is a matrix
  counts <- counts %>% data.matrix()
  ## Compute offset (with optional parameters)
  offset_function(counts, ...)
}

# Prepare data for use in PLN models from a biom object
#
# @description Wrapper around \code{\link[=prepare_data]{prepare_data}}, extracts the count table and the covariates data.frame from a "biom" class object
# before passing them to \code{\link[=prepare_data]{prepare_data}}. See \code{\link[=prepare_data]{prepare_data}} for details.
#
# @param biom Required. Either a biom-class object from which the count table and covariates data.frame are extracted or a file name where to read the biom.
# @inheritParams prepare_data
# @param ... Addtional arguments passed on to \code{\link[=compute_offset]{compute_offset}}
#
# @seealso \code{\link[=compute_offset]{compute_offset}} and \code{\link[=prepare_data]{prepare_data}}
# @export
#
# @details This functions depends on the biomformat package which is not a proper dependency of PLNmodels as it is not available on CRAN
#
# @importFrom biomformat read_biom sample_metadata biom_data
# @examples
# ## Requires the biomformat package
# \dontrun{
# library(biomformat)
# biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
# biom <- read_biom(biom_file)
# prepare_data_from_biom(biom)
# }
# prepare_data_from_biom <- function(biom, offset = "TSS", ...) {
#   if (is.character(biom)) biom <- biomformat::read_biom(biom)
#   sdf <- biomformat::sample_metadata(biom)
#   if (is.null(sdf) || all(is.na(sdf))) {
#     stop(paste("No covariates detected in biom. Consider:",
#                "- extracting count data from biom with biom_data()",
#                "- preparing a covariates data.frame",
#                "- using prepare_data instead of prepare_data_from_biom",
#                sep = "\n"))
#   }
#   prepare_data(counts     = biomformat::biom_data(biom) %>% as("matrix"),
#                covariates = sdf,
#                offset     = offset,
#                ...)
# }

# Prepare data for use in PLN models from a phyloseq object
#
# @description Wrapper around \code{\link[=prepare_data]{prepare_data}}, extracts the count table and the covariates data.frame from a "phyloseq" class object
# before passing them to \code{\link[=prepare_data]{prepare_data}}. See \code{\link[=prepare_data]{prepare_data}} for details.
#
# @param physeq Required. A phyloseq class object from which the count table and covariates data.frame are extracted.
# @inheritParams prepare_data
# @param ... Addtional arguments passed on to \code{\link[=compute_offset]{compute_offset}}
#
# @seealso \code{\link[=compute_offset]{compute_offset}} and \code{\link[=prepare_data]{prepare_data}}
# @export
#
# @details This functions depends on the phyloseq package which is not a proper dependency of PLNmodels as it is not available on CRAN
#
# @importFrom phyloseq sample_data otu_table
# @examples
# ## Requires the phyloseq package
# \dontrun{
# library(phyloseq)
# data(enterotype)
# prepare_data_from_phyloseq(enterotype)
# }
# prepare_data_from_phyloseq <- function(physeq, offset = "TSS", ...) {
#   if (!inherits(physeq, "phyloseq")) stop("physeq should be a phyloseq object.")
#   if (is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))) {
#     stop(paste("No covariates detected in physeq Consider:",
#                "- extracting count data from biom with as(otu_table(physeq), \"matrix\")",
#                "- preparing a covariates data.frame",
#                "- using prepare_data instead of prepare_data_from_phyloseq",
#                sep = "\n"))
#   }
#   prepare_data(counts     = phyloseq::otu_table(physeq) %>% as("matrix"),
#                covariates = phyloseq::sample_data(physeq) %>% as("data.frame"),
#                offset     = offset,
#                ...)
# }

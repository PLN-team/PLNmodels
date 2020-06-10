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
  ## Geometric mean of finite values
  geom_mean <- function(x) {
    x_log <- log(x)
    exp(mean(x_log[is.finite(x_log)], na.rm = TRUE))
  }
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
offset_rle <- function(counts, pseudocounts = 0) {
  ## Add pseudo.counts
  counts <- counts + pseudocounts
  ## compute geometric mean for all otus
  geom_means <- counts %>% log() %>% colMeans(na.rm = TRUE) %>% exp
  if (all(geom_means == 0)) stop("Sample do not share any common species, RLE normalization failed.")
  ## compute size factor as the median of all otus log-ratios in that sample
  robust_med <- function(cnts) { median((cnts/geom_means)[is.finite(geom_means) & cnts > 0]) }
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

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## EXPORTED FUNCTIONS ---------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Prepare data for use in PLN models
#' @name prepare_data
#'
#' @description Prepare data in proper format for use in PLN model and its variants. The function (i) merges a count table and
#' a covariate data frame in the most comprehensive way and (ii) computes offsets from the count table using one of several normalization schemes (TSS, CSS, RLE, GMPR, etc). The function fails with informative messages when the heuristics used for sample matching fail.
#'
#' @param counts Required. An abundance count table, preferably with dimensions names and species as columns.
#' @param covariates Required. A covariates data frame, preferably with row names.
#' @param offset Optional. Normalization scheme used to compute scaling factors used as offset during PLN inference. Available schemes are "TSS" (Total Sum Scaling, default), "CSS" (Cumulative Sum Scaling, used in metagenomeSeq), "RLE" (Relative Log Expression, used in DESeq2), "GMPR" (Geometric Mean of Pairwise Ratio, introduced in Chen et al., 2018) or "none". Alternatively the user can supply its own vector or matrix of offsets (see note for specification of the user-supplied offsets).
#' @param ... Additional parameters passed on to [compute_offset()]
#'
#' @references Chen, L., Reeve, J., Zhang, L., Huang, S., Wang, X. and Chen, J. (2018) GMPR: A robust normalization method for zero-inflated count data with application to microbiome sequencing data. PeerJ, 6, e4600 \url{https://doi.org/10.7717/peerj.4600}
#' @references Paulson, J. N., Colin Stine, O., Bravo, H. C. and Pop, M. (2013) Differential abundance analysis for microbial marker-gene surveys. Nature Methods, 10, 1200-1202 \url{http://dx.doi.org/10.1038/nmeth.2658}
#' @references Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data. Genome Biology, 11, R106 \url{https://doi.org/10.1186/gb-2010-11-10-r106}
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
#' @details RLE has an additional `pseudocounts` arguments to add pseudocounts to the observed counts (defaults to 0). CSS has an additional `reference` argument to choose the location function used to compute the reference quantiles (defaults to `median` as in the Nature publication but can be set to `mean` to reproduce behavior of functions cumNormStat* from metagenomeSeq). Note that (i) CSS normalization fails when the median absolute deviation around quantiles does not become instable for high quantiles (limited count variations both within and across samples) and/or one sample has less than two positive counts, (ii) RLE fails when there are no common species across all samples and (iii) GMPR fails if a sample does not share any species with all other samples.
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
compute_offset <- function(counts, offset = c("TSS", "GMPR", "RLE", "CSS", "none"), ...) {
  ## special behavior for numeric offset
  if (is.numeric(offset)) {
    return(offset_numeric(counts, offset, ...))
  }
  ## Choose offset function
  offset <- match.arg(offset)
  offset_function <- switch(offset,
                            "TSS"  = offset_tss,
                            "GMPR" = offset_gmpr,
                            "RLE"  = offset_rle,
                            "CSS"  = offset_css,
                            "none" = offset_none
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

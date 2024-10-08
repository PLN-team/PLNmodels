## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  INTERNAL FUNCTIONS ---------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Internal function to find the most comprehensive set of common samples between a count table and a covariates data.frame
common_samples <- function(counts, covariates, call = rlang::caller_env()) {
  ## Have default samples names been created when matching samples?
  default_names <- FALSE
  ## Sanity checks:
  name_warning <- c(
    "!" = "There are no matching names in {.var counts} and {.var covariates}.",
    "i" = "Function will proceed assuming:",
    "i" = "- samples are in the same order;",
    "i" = "- samples are rows of {.var counts}.", sep = "\n"
  )
  row_count_abort <- c(
    "x" = "{.var counts} and {.var covariates} have different number of row(s):",
    "i" = "{.var counts} has {.cls {nrow(counts)}} row(s);",
    "i" = "{.var covariates} has {.cls {nrow(covariates)}} row(s).",
    sep = "\n"
  )
  ## no sample names in covariates: create sample names
  if (is.null(rownames(covariates))) {
    cli::cli_warn(name_warning, call = call)
    if (nrow(counts) != nrow(covariates)) {
      cli::cli_abort(row_count_abort, call = call)
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
    cli::cli_warn(name_warning, call = call)
    if (nrow(counts) != nrow(covariates)) {
      cli::cli_abort(row_count_abort, call = call)
    }
    if (!is.null(rownames(counts))) {
      cli::cli_abort(c(
        "x" = "Conflicting sample names in {.var counts} matrix and {.var covariates} data frames",
        "i" = "Sample names in {.var counts} matrix is {.cls {rownames(counts)}} and in {.var covariates} is {.cls {rownames(covariates)}}."
      ), call = call)
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
  # Add condition when less counts than covariates information
  if (length(common_samples) < nrow(counts)) {
    cli::cli_warn(c(
      "!" = "There are less samples in {.var counts} than in {.var covariates}.",
      "i" = "{.cls {nrow(counts) - length(common_samples)}} samples were dropped from the {.var counts} matrix for lack of associated {.var covariates}.",
      "i" = "{.cls There {?is/are} {length(common_samples)}} sample{?s} in the final data.frame."
    ), call = call)
  }
  if (length(common_samples) < nrow(covariates)) {
    cli::cli_warn(c(
      "!" = "There are less samples in {.var covariates} than in {.var counts}.",
      "i" = "{.cls {nrow(covariates) - length(common_samples)}} samples were dropped from {.var covariates} for lack of associated {.var counts}.",
      "i" = "{.cls There {?is/are} {length(common_samples)}} sample{?s} in the final data.frame."
    ), call = call)
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

## Geometric Mean Pairwise Ratio (GMPR) normalization (as presented in \doi{10.7717/peerj.4600})
offset_gmpr <- function(counts) {
  if (nrow(counts) == 1) stop("GMPR is not defined when there is only one sample.")
  ## median of pairwise ratios between counts of samples i and j, limited to positive counts
  pairwise_ratio <- function(i, j) {
    c_i <- counts[i, ]; c_j <- counts[j, ]
    ratio <- c_i / c_j
    median(ratio[c_i > 0 & c_j > 0])
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
  size_factor
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
  size_factor
}

## Trimmed Mean of M-values (TMM) normalization (as used in edgeR and presented in \doi{10.1186/gb-2010-11-3-r25})
# We choose by default refColumn = NULL, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff
calcFactorTMM <- function(obs, ref, nO, nR, logratioTrim=.3, sumTrim=0.05, Acutoff=-1e10)
  #	TMM between two libraries simplified and adapted from M. Robinson (edgeR:::.calcFactorTMM)
  # The final output is different from original TMM as we directly multiply normalization factors
  # by library sizes to use them as such in the model, unlike the philosophy behind TMM
{
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance

  #	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  if(max(abs(logR)) < 1e-6) return(1)

  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS

    keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

 # doWeighting = TRUE
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  #	Results will be missing if the two libraries share no features with positive counts
  #	In this case, return unity
  if(is.na(f)) f <- 0
  2^f
}

offset_tmm <- function(counts, logratioTrim=.3, sumTrim=0.05, Acutoff=-1e10) {
  nsamples <- nrow(counts)
  ## Compute lib.size
  lib_size <- rowSums(counts)
  ## Reference sample calculated from the .75 quantile
  f75 <- apply(counts, 1, FUN = function(g) quantile(g, p=0.75)/sum(g))
  if(median(f75) < 1e-20) {
    refColumn <- which.max(rowSums(sqrt(counts)))
  } else {
    refColumn <- which.min(abs(f75-mean(f75)))
  }
  ## Compute TMM normalization factor
  ref <- as.numeric(counts[refColumn,])
  nR <- lib_size[refColumn]
  # rel_counts <- counts/lib_size
  f <- vapply(seq(nsamples), FUN.VALUE = numeric(1),FUN = function(i) calcFactorTMM(obs=counts[i,],ref=ref, nO=lib_size[i], nR=nR, logratioTrim=logratioTrim, sumTrim=sumTrim, Acutoff=Acutoff))
  #	Factors should multiple to one
  f <- f * lib_size
  f / geom_mean(f)
}


## Cumulative Sum Scaling (CSS) normalization (as used in metagenomeSeq and presented in \doi{10.1038/nmeth.2658})
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
  unname(size_factors)
}

## Wrench normalization (from \doi{10.1186/s12864-018-5160-5}) in its simplest form
## @importFrom matrixStats rowWeightedMeans rowVars
#' @importFrom stats binomial var
offset_wrench <- function(counts, groups = rep(1, nrow(counts)), type = c("wrench", "simple")) {
  ## Helpers and preprocessing
  log_var <- function(x) { var(log(x[is.finite(x) & x > 0])) }
  ## n = number of samples, p = number of species
  n <- nrow(counts); p <- ncol(counts)
  if (n == 1) stop("Wrench is not defined when there is only one sample.")
  if (p == 1) {
    depths <- rowSums(counts)
    return(depths / geom_mean(depths))
  }
  type <- match.arg(type)

  ## Proportions
  depths <- rowSums(counts)
  props  <- counts / depths
  ## Proportions in reference
  props_ref <- colMeans(props)
  ## Samplewise ratios of proportions
  sample_ratios <- props / matrix(props_ref, n, p, byrow = TRUE)

  ## Species variance
  species_vars <- species_variance(counts, groups, depths_as_offset = (type == 'simple'))

  ## Sample specific variances and scales, computed at the group level.

  ## Sample variance and scales (shared within groups)
  K = length(unique(groups)) ## number of groups
  if (K == 1) {
    global_ratio <- (colSums(counts) / sum(counts)) / props_ref
    sample_scales <- rep(mean(global_ratio), n)
    sample_vars   <- rep(log_var(global_ratio), n)
  } else {
    groups <- as.character(groups)
    ## Groupwise counts, proportions and ratios
    group_counts <- rowsum(counts, groups, reorder = TRUE)
    group_props  <- group_counts / rowSums(group_counts)
    group_ratios <- group_props / matrix(props_ref, K, p, byrow = TRUE)

    ## groupwise scale and (log)dispersion factors
    group_scales <- rowMeans(group_ratios)
    group_vars   <- apply(group_ratios, 1, log_var)

    ## sample scales and (log)dispersion factors
    sample_scales <- group_scales[groups]
    sample_vars   <- group_vars[groups]
  }

  ## Regularized estimation of positive means
  ## Group centered log-ratios
  log_sample_ratios <- log(sample_ratios / sample_scales)

  ## a_i is the mixed effect of sample i = weighted mean of group-corrected log-ratio in the sample.
  ## Weights inversely proportional to sample var and species variance
  ## weights[sample, species] = 1 / (sample_vars[sample] + species_vars[species]) (or NA is corresponding counts is NULL)
  weights <- 1 / outer(sample_vars, species_vars, "+")
  weights[] <- weights / rowSums(weights)
  weights[!is.finite(log_sample_ratios)] <- NA
  a_i <- rowSums(log_sample_ratios * weights, na.rm = TRUE)

  ## b_ij is the mixed effet of species j in sample i = shrunken group and sample corrected log-ratio
  # shrinkage[sample, species] = sample_vars[sample] / (sample_vars[sample] + species_vars[species])
  shrinkage <- sample_vars / outer(sample_vars, species_vars, "+")
  b_ij <- shrinkage * (log_sample_ratios - a_i)

  ## Sample ratios theta_{ij} = group mean * exp(a_i + b_ij)
  sample_ratios[] <- exp(b_ij + a_i) * sample_scales

  ## Averaging weights
  if (type == "wrench") {
    ## Use wrench defaults:
    ## - adjustment of ratio by variance terms
    sample_ratios[] <- sample_ratios / matrix(exp(species_vars / 2), n, p, byrow = TRUE)
    ## - computation of probability of absence by fitting a simple binomial model absence ~ log(depths) to each species
    pi0 <- apply(counts, 2, function(y) { suppressWarnings(glm.fit(cbind(log(depths)), 0L + (y == 0), family = binomial())$fitted.values) } )
    ## - weights proportional to hurdle and variance
    weights <- exp(outer(sample_vars, species_vars, "+")) - 1
    weights[] <- 1 / ( (1 - pi0) * (pi0 + weights) )
    weights[!is.finite(weights)] <- NA
  } else {
    ## Simple implementation: uniform weights
    weights <- matrix(1, nrow = n, ncol = p, byrow = TRUE)
  }

  ## Compositional factors: robust weighted means of ratio
  weights[] <- weights / rowSums(weights, na.rm = TRUE)
  comp_factors <- rowSums(sample_ratios * weights, na.rm = TRUE)
  comp_factors[] <- comp_factors / geom_mean(comp_factors) ## Centered in geometric scale

  ## Normalization factors
  norm_factors <- comp_factors * depths / geom_mean(depths)

  unname(norm_factors)
}

# Helpers scaling functions ----

## Geometric mean (computed only on positive values)
geom_mean <- function(x, poscounts = TRUE, na.rm = TRUE) {
  x_log <- log(x)
  if (poscounts) x_log <- x_log[x > 0]
  exp(mean(x_log, na.rm = na.rm))
}

## Transform scaling factor to normalized offsets (on the count scale)
sf2nf <- function(scaling_factors, lib_size){
  if (is.null(scaling_factors)) return(NULL)
  tmp <- scaling_factors / lib_size
  (tmp / geom_mean(tmp)) * lib_size
}

species_variance <- function(counts, groups = rep(1, nrow(counts)), depths_as_offset = TRUE) {
  ## n = number of samples, p = number of species
  n <- nrow(counts); p <- ncol(counts)

  ## Centered log depths and counts corrected by offset
  log_depths <- log(rowSums(counts))
  log_counts <- log(counts)
  log_counts[!is.finite(log_counts)] <- NA

  ## Design matrix
  groups <- as.character(groups)
  if (length(unique(groups)) > 1) {
    design <- model.matrix(counts[, 1] ~ 0 + groups)
  } else {
    design <- model.matrix(counts[, 1] ~ 1)
  }

  ## Depth status: offset or covariate
  if (depths_as_offset) {
    log_counts <- log_counts - log_depths ## log depths as an offset
  } else {
    design <- cbind(design, log_depths)   ## log depths as a covariate
  }


  ## Manually regress log_counts against log_depths (for each species) to compute species-specific sigma.
  #' @importFrom stats .lm.fit
  compute_sigma <- function(j) {
    valid_obs <- !is.na(log_counts[, j])
    model <- .lm.fit(design[valid_obs, , drop = FALSE], log_counts[valid_obs, j])
    if (model$rank == sum(valid_obs)) return(0)
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
#' @param call Optional. The execution environment in which to set the local error call.
#' @param ... Additional parameters passed on to [compute_offset()]
#'
#' @references Chen, L., Reeve, J., Zhang, L., Huang, S., Wang, X. and Chen, J. (2018) GMPR: A robust normalization method for zero-inflated count data with application to microbiome sequencing data. PeerJ, 6, e4600 \doi{10.7717/peerj.4600}
#' @references Paulson, J. N., Colin Stine, O., Bravo, H. C. and Pop, M. (2013) Differential abundance analysis for microbial marker-gene surveys. Nature Methods, 10, 1200-1202 \doi{10.1038/nmeth.2658}
#' @references Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data. Genome Biology, 11, R106 \doi{10.1186/gb-2010-11-10-r106}
#' @references Kumar, M., Slud, E., Okrah, K. et al. (2018) Analysis and correction of compositional bias in sparse sequencing count data. BMC Genomics 19, 799 \doi{10.1186/s12864-018-5160-5}
#' @references Robinson, M.D., Oshlack, A. (2010) A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11, R25 \doi{10.1186/gb-2010-11-3-r25}
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
#'  offset     = "GMPR",
#'  scale      = "count"
#' )
#' proper_data$Abundance
#' proper_data$Offset
prepare_data <- function(counts, covariates, offset = "TSS", call = rlang::caller_env(), ...) {
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
    cli::cli_warn(c(
      "!" = "There is at least one NA value in {.var counts}.",
      "i" = "{.cls {sum(is.na(counts))}} NA value{?s} in {.var counts} {?has/have} been replaced with 0."
    ), call = call)
    counts[is.na(counts)] <- 0
  }
  ## filter out empty samples
  empty_samples <- which(rowSums(counts) == 0)
  if (length(empty_samples) > 0) { # Add > 0
    cli::cli_warn(c(
      "!" = "There  is at least one empty sample in {.var counts}.",
      "i" = "{.cls {length(empty_samples)}} sample{?s} ({.cls {samples[empty_samples]}}) in {.var counts} {?has/have} been dropped for lack of positive counts."
    ), call = call)
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
  result
}


#' @title Compute offsets from a count data using one of several normalization schemes
#' @name compute_offset
#'
#' @description Computes offsets from the count table using one of several normalization schemes (TSS, CSS, RLE, GMPR, Wrench, TMM, etc) described in the literature.
#'
#' @inheritParams prepare_data
#' @param scale Either `"none"` (default) or `"count"`. Should the offset be normalized to be on the same scale as the counts ?
#' @param ... Additional parameters passed on to specific methods (for now CSS and RLE)
#' @inherit prepare_data references
#'
#' @details RLE has additional `pseudocounts` and `type` arguments to add pseudocounts to the observed counts (defaults to 0L) and to compute offsets using only positive counts (if `type == "poscounts"`). This mimics the behavior of \code{DESeq2::DESeq()} when using `sfType == "poscounts"`. CSS has an additional `reference` argument to choose the location function used to compute the reference quantiles (defaults to `median` as in the Nature publication but can be set to `mean` to reproduce behavior of functions cumNormStat* from metagenomeSeq). Wrench has two additional parameters: `groups` to specify sample groups and `type` to either reproduce exactly the default \code{Wrench::wrench()} behavior (`type = "wrench"`, default) or to use simpler heuristics (`type = "simple"`). Note that (i) CSS normalization fails when the median absolute deviation around quantiles does not become instable for high quantiles (limited count variations both within and across samples) and/or one sample has less than two positive counts, (ii) RLE fails when there are no common species across all samples (unless `type == "poscounts"` has been specified) and (iii) GMPR fails if a sample does not share any species with all other samples.
#' TMM code between two libraries is simplified and adapted from M. Robinson (edgeR:::.calcFactorTMM).
#' The final output is however different from the one produced by edgeR:::.calcFactorTMM as they are intended
#' to be used as such in the model (whereas they need to be multiplied by sequencing depths in edgeR)
#'
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
#' compute_offset(counts, offset = "RLE", pseudocounts = 1)
#' compute_offset(counts, offset = "Wrench", groups = trichoptera$Covariate$Group)
#' compute_offset(counts, offset = "GMPR")
#' compute_offset(counts, offset = "TMM")
#' ## User supplied offsets
#' my_offset <- setNames(rep(1, nrow(counts)), rownames(counts))
#' compute_offset(counts, offset = my_offset)
compute_offset <- function(counts, offset = c("TSS", "GMPR", "RLE", "CSS", "Wrench", "TMM", "none"), scale = c("none", "count"), ...) {
  ## special behavior for data.frame
  if (inherits(offset, "data.frame")) {
    cli::cli_abort(c(
      "{.var offset} must be an available scheme or a vector or matrix of offsets.",
      "x" = "You supplied a data.frame for {.var offset}",
      "i" = "Did you mean to supply a numeric matrix?",
      "i" = "Try converting your data.frame to a matrix with `as.matrix()`."
    ))
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
                            "TMM" = offset_tmm,
                            "none"   = offset_none
  )
  ## Ensure that counts is a matrix
  counts <- counts %>% data.matrix()
  ## Compute offset (with optional parameters)
  scale <- match.arg(scale)
  if (scale == "none") {
    offset_function(counts, ...)
  } else {
    lib_size <- offset_tss(counts)
    sf2nf(offset_function(counts, ...), lib_size = lib_size)
  }


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

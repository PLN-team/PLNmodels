###########################
##  INTERNAL FUNCTIONS   ##
###########################

## Internal function to find the most comprehensive set of common samples between a count table and a covariates data.frame
common_samples <- function(counts, covariates) {
  ## Sanity checks
  if (is.null(dimnames(counts)) || is.null(rownames(covariates))) {
    warning(paste("There are no dimension names for the abundance matrix and/or the covariates data.frame.",
                  "Function will proceed assuming:",
                  "- samples are in the same order;",
                  "- samples are rows of the abundance matrix.", sep = "\n"))
    if (nrow(counts) != nrow(covariates)) {
      stop("Incompatible dimensions")
    }
    rownames(counts) <- rownames(covariates) <- paste0("Sample_", 1:nrow(counts))
  }
  ## Check whether samples are stored as columns in the abundance matrix
  ## and transpose if that's the case
  ## Based on a heuristic of matching names
  sample_are_cols <- any(colnames(counts) %in% rownames(covariates))
  if (sample_are_cols) counts <- t(counts)
  sample_are_rows <- any(rownames(counts) %in% rownames(covariates))
  if (!sample_are_rows) stop("No matching samples in abundance matrix and covariates data.frame, please check the input.")
  ## Ensure consistency by using only common samples
  common_samples <- intersect(rownames(counts), rownames(covariates))
  if (length(common_samples) < nrow(counts)) {
    message(paste0(nrow(counts) - length(common_samples), " samples were dropped from the abundance matrix for lack of associated covariates."))
  }
  return(list(transpose_counts = sample_are_cols,
              common_samples   = common_samples))
}

## Internal functions to compute scaling factors from a count table

## No offset
offset_none <- function(counts) {
  return(NULL)
}

## Total Sum Scaling offset
offset_tss <- function(counts) {
  rowSums(counts)
}

## Geometric Mean Pairwise Ratio (GMPR) normalisation (as presented in doi.org/10.7717/peerj.4600)
offset_gmpr <- function(counts) {
  ## median of (non-null, non-infinite) pairwise ratios between counts of samples i and j
  pairwise_ratio <- function(i, j) { median(counts[i, ] / counts[j, ], na.rm = TRUE) }
  ## Robust geometric mean
  geom_mean <- function(x) { x %>% log %>% mean(na.rm = TRUE) %>% exp }
  ## Matrix of pairwise ratios
  n <- nrow(counts)
  mat_pr <- matrix(NaN, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      mat_pr[i, j] <- pairwise_ratio(i, j)
    }
  }
  ## Geometric mean of pairwise ratio
  size_factor <- apply(mat_pr, 1, geom_mean)
  return(size_factor)
}

## Relative Log Expression (RLE) normalization (as used in DESeq2)
offset_rle <- function(counts, pseudocounts = 0) {
  ## Add pseudo.counts and compute geometric for all otus
  geom_means <- (counts + pseudocounts) %>% log() %>% colMeans(na.rm = TRUE) %>% exp
  if (all(is.nan(geom_means))) stop("There are no common OTUs, RLE normalization failed")
  ## compute size factor as the median of all otus log-ratios in that sample
  size_factor <- apply(counts, 1, function(x) {median(x / geom_means, na.rm = TRUE)} )
  return(size_factor)
}

## Cumulative Sum Scaling (CSS) normalization (as used in metagenomeSeq and presented in doi.org/10.1038/nmeth.2658)
offset_css <- function(counts) {
  ## matrix of sample-specific quantiles and cumulative sums up to quantile
  mat_sample_quant <- apply(counts, 1, sort) %>% t()
  mat_sample_cumsum <- apply(mat_sample_quant, 1, cumsum) %>% t()
  ## reference quantiles, computed as median of sample_specific quantiles and MAD around the reference quantiles
  ref_quant_mad <- apply(mat_sample_quant, 1, mad, constant = 1)
  ## find smallest quantile for which high instability is detected
  ## instability for quantile l is defined as ref_quant_mad[l+1] - ref_quant_mad[l] >= 0.1 * ref_quant_mad[l]
  instable <- (diff(ref_quant_mad) >= 0.1 * head(ref_quant_mad, -1))
  lhat <- min(which(instable))
  if (!is.finite(lhat)) {
    warning("No instability detected in quantile distribution across samples, falling back to scaled TSS normalization.")
    lhat <- ncol(counts)
  }
  ## scaling factors are cumulative sums up to quantile lhat, divided by their median
  size_factors <- mat_sample_cumsum[ , lhat] / median(mat_sample_cumsum[ , lhat])
  return(size_factors)
}

###########################
##  EXPORTED FUNCTIONS   ##
###########################

prepare_data <- function(counts, covariates, offset = "TSS", ...) {
  ## sanitize abundance matrix and covariates data.frame
  common <- common_samples(counts, covariates)
  samples <- common$common_samples
  if (common$transpose_counts) counts <- t(counts)
  counts     <- counts[samples, ]
  ## Replace NA with 0s
  if (any(is.na(counts))) {
    counts[is.na(counts)] <- 0
    warning("NA values in count table replaced with 0.")
  }
  ## filter out empty samples
  empty_samples <- which(rowSums(counts) == 0)
  if (length(empty_samples)) {
    warning(paste0("Samples ", samples[empty_samples], " were dropped as they have no positive counts."))
    samples <- samples[-empty_samples]
    counts <- counts[samples, ,drop = FALSE]
  }
  covariates <- covariates[samples, ] %>% as.data.frame
  if (is.null(names(covariates))) names(covariates) <- paste0("Variable", seq_along(covariates))
  ## compute offset
  offset     <- compute_offset(counts, offset, ...)
  ## prepare data for PLN
  result <- data.frame(Abundance = I(counts),
                       covariates)
  if (!is.null(offset)) result$Offset <- I(offset)
  return(result)
}

compute_offset <- function(counts, offset = c("TSS", "GMPR", "RLE", "CSS", "none"), ...) {
  ## Choose offset function
  offset <- match.arg(offset)
  offset_function <- switch(offset,
                            "TSS"  = offset_tss,
                            "GMPR" = offset_gmpr,
                            "RLE"  = offset_rle,
                            "CSS"  = offset_css,
                            "none" = offset_none
  )
  ## Compute offset (with optional parameters)
  offset_function(counts, ...)
}

## Example biom class object
# library(biomformat)
# biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
# biom_file
# biom <- read_biom(biom_file)

# Example phyloseq class object
# library(phyloseq)
# data(enterotypes)

prepare_data_from_biom <- function(biom, offset = "TSS", ...) {
  if (is.character(biom)) biom <- read_biom(biom)
  if (is.null(sample_metadata(biom))) {
    stop(paste("No covariates detected in biom. Consider:",
               "- extracting count data from biom with biom_data()",
               "- preparing a covariates data.frame",
               "- using prepare_data instead of prepare_data_from_biom",
               sep = "\n"))
  }
  prepare_data(counts     = biom_data(biom) %>% as("matrix"),
               covariates = sample_metadata(biom),
               offset     = offset,
               ...)
}

prepare_data_from_phyloseq <- function(physeq, offset = "TSS", ...) {
  if (!inherits(physeq, "phyloseq")) stop("physeq should be a phyloseq object.")
  if (is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))) {
    stop(paste("No covariates detected in physeq Consider:",
               "- extracting count data from biom with as(otu_table(physeq), \"matrix\")",
               "- preparing a covariates data.frame",
               "- using prepare_data instead of prepare_data_from_phyloseq",
               sep = "\n"))
  }
  prepare_data(counts     = phyloseq::otu_table(physeq) %>% as("matrix"),
               covariates = phyloseq::sample_data(physeq),
               offset     = offset,
               ...)
}



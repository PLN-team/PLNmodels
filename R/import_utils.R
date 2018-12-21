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

prepare_data <- function(counts, covariates, offset = c("TSS", "GMPR", "none")) {
  ## Offset setup
  offset <- match.arg(offset)
  offset_function <- switch (offset,
                             "TSS"  = offset_tss,
                             "GMPR" = offset_gmpr,
                             "none" = offset_none
  )
  ## sanitize abundance matrix and covariates data.frame
  common <- common_samples(counts, covariates)
  if (common$transpose_counts) counts <- t(counts)
  counts     <- counts[common$common_samples, ]
  covariates <- covariates[common$common_samples, ] %>% as.data.frame
  if (is.null(names(covariates))) names(covariates) <- paste0("Variable", seq_along(covariates))
  ## compute offset
  offset     <- offset_function(counts)
  ## prepare data for PLN
  result <- data.frame(Abundance = I(counts),
                       covariates)
  if (!is.null(offset)) results$offset <- I(offset)
  return(result)
}

offset_none <- function(counts) { return(NULL) }

offset_tss <- offset_gmpr <- offset_none

## prepare_data(counts, covariates) %>% str()


## TODO: write functions to prepare data from a biom or a phyloseq-class object.

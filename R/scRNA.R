#' Single cell RNA-seq data
#'
#' A dataset containing the counts of the 500 most varying transcripts in the mixtures of
#' 5 cell lines in human liver (obtained with standard 10x scRNAseq Chromium protocol).
#'
#' @format A data frame named 'scRNA' with 3918 rows (the cells) and 3 variables:
#' \describe{
#'   \item{counts}{a 500 trancript by 3918 count matrix}
#'   \item{cell_line}{factor, the cell line of the current row (among 5)}
#'   \item{total_counts}{Total number of reads for that cell}
#'   ...
#' }
#' @source \url{https://github.com/LuyiTian/sc_mixology/}
"scRNA"

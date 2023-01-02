#' Barents fish data set
#'
#' @name barents
#'
#' @description This data set gives the abundance of 30 fish species observed in 89 sites in the Barents sea.
#' For each site, 4 additional covariates are known. Subsample of the original datasets studied by Fossheim et al, 2006.
#'
#' @format A data frame with 6 variables:
#' * Abundance: A 30 fish species by 89 sites count matrix
#' * Offset: A 30 fish species by 116 samples offset matrix, measuring the sampling effort in each site
#' * 4 covariates for latitude, longitude, depth (in meters), temperature (in Celsius degrees).
#'
#' @references Fossheim, Maria, Einar M. Nilssen, and Michaela Aschan. "Fish assemblages in the Barents Sea." Marine Biology Research 2.4 (2006). \doi{10.1080/17451000600815698}
#' @examples
#' data(barents)
#' @source Data from M. Fossheim and coauthors.
"barents"


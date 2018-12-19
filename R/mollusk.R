##' Mollusk data set
##'
##' This data set gives the abundance of 32 mollusk species in 163 samples.
##' For each sample, 4 additional informations (covariates) are knowne.
##'
##' @format A data frame with 163 rows and 5 variables (4 vectors and 1 matrix):
##' \describe{
##'   \item{Abundance}{a 163 x 32 matrix of abundancies/counts (163 samples and 32 mollusk species)}
##'   \item{site}{a factor with 8 levels indicating the sampling site}
##'   \item{season}{a factor with 4 levels indicating the season}
##'   \item{method}{a factor with 2 levels for the method of sampling (wood or string)}
##'   \item{duration}{a numeric with 3 levels for the time of exposure in week}
##' }
##'
##' This format is convenient for using formula in multivariate analysis (multiple outputs and inputs).
##' Original data set has been extracted ade4 and formatted for PLNmodels.
##'
##' @examples
##' data(mollusk)
##' str(mollusk)
##' ## also see the package case study about mollusk
##'
##' @references Richardot-Coulet, M., Chessel D. and Bournaud M. (1986) Typological value of the benthos of old beds of a large river. Methodological approach. Archiv fùr Hydrobiologie, 107, 363–383.
##' @source Data from Richardot-Coulet, Chessel and Bournaud.
"mollusk"


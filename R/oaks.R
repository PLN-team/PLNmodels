#' Oaks amplicon data set
#'
#' @name oaks
#'
#' @description This data set gives the abundance of 114 taxa (66 bacterial OTU, 48 fungal OTUs) in 116 samples.
#' For each sample, 11 additional covariates are known.
#'
#' @format A data frame with 13 variables:
#' * Abundance: A 114 taxa by 116 samples count matrix
#' * Offset: A 114 taxa by 116 samples offset matrix
#' * Sample: Unique sample id
#' * tree: Tree status with respect to the pathogen (susceptible, intermediate or resistant)
#' * branch: Unique branch id in each tree (4 branches were sampled in each tree, with 10 leaves per branch)
#' * leafNO: Unique leaf id in each tree (40 leaves were sampled in each tree)
#' * distTObase: Distance of the sampled leaf to the base of the branch
#' * distTOtrunk: Distance of the sampled leaf to the base of the tree trunk
#' * distTOground: Distance of the sampled leaf to the base of the ground
#' * pmInfection: Powdery mildew infection, proportion of the upper leaf area displaying mildew symptoms
#' * orientation: Orientation of the branch (South-West SW or North-East NE)
#' * readsTOTfun: Total number of ITS1 reads for that leaf
#' * readsTOTbac: Total number of 16S reads for that leaf
#'
#'
#' @seealso [prepare_data()]
#' @references Jakuschkin, B., Fievet, V., Schwaller, L. et al. Deciphering the Pathobiome: Intra- and Interkingdom Interactions Involving the Pathogen Erysiphe alphitoides . Microb Ecol 72, 870–880 (2016). \doi{10.1007/s00248-016-0777-x}
#' @examples
#' data(oaks)
#' \dontrun{
#' oaks_networks <- PLNnetwork(formula = Abundance ~ 1 + offset(log(Offset)), data = oaks)
#' }
#' @source Data from B. Jakuschkin and coauthors.
"oaks"


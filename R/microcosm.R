#' Cow microbiome data set
#'
#' @name microcosm
#'
#' @description This data gives the evolution of the microbiota of 45 lactating cows
#' before and after calving in terms of counts of Amplicon Sequence Variants (ASV) in various
#' body sites. Three body sites (vagina, mouth, nose) in addition to the milk of the 4 teats were sampled
#' at 4 times points: 1 week before calving (except for the milk), 1 month, 3 months and 7 months after calving.
#' We present a reduced version of the orginal data consisting of 880 samples with
#' sequencing depths ranging from 1,001 to 71,881 reads per sample (see [Mariadassou et al. (2023)](https://doi.org/10.1186/s42523-023-00252-w) for details). ASVs with prevalence in the dataset
#' lower than 5% were filtered out and samples for which the total depth (after ASV filtering were removed),
#' resulting in a count table of n = 880 samples with p = 259 ASV and a mean proportion of zeroes of 90.3%.
#'
#' @format A data frame with 6 variables:
#' * Abundance: A 880 samples by 259 taxa count matrix
#' * sample: Unique sample id
#' * Offset: sequencing depth
#' * site: sampling site (`O`: oral; `N`: nasal; `V`: vaginam; `M`: milk)
#' * time: sampling time (`-1W`: 1 week before calving; `1M`: 1 month after calving; `3M`: 3 months after calving; `7M`: 7 months after calving)
#' * site_time: factor of possible pairs of (site, time). The combination `M_-1W` is absent.
#'
#' @seealso [prepare_data()]
#' @references  Mariadassou, M., Nouvel, L.X., Constant, F. et al. Microbiota members from body sites of dairy cows are largely shared within individual hosts throughout lactation but sharing is limited in the herd. anim microbiome 5, 32 (2023). \doi{10.1186/s42523-023-00252-w}
#' @examples
#' data(microcosm)
#' \dontrun{
#' my_ZIPLN <- ZIPLN(formula = Abundance ~ 0 + site + offset(log(Offset)) | 1, data = microcosm)
#' }
#' @source Data from M. Mariadassou and coauthors.
"microcosm"


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ZIPLNfit-S3methods.R
\name{predict.ZIPLNfit}
\alias{predict.ZIPLNfit}
\title{Predict counts of a new sample}
\usage{
\method{predict}{ZIPLNfit}(
  object,
  newdata,
  responses = NULL,
  level = 1,
  type = c("link", "response", "deflated"),
  ...
)
}
\arguments{
\item{object}{an R6 object with class \code{\link{ZIPLNfit}}}

\item{newdata}{A data frame in which to look for variables and offsets with which to predict}

\item{responses}{Optional data frame containing the count of the observed variables (matching the names of the provided as data in the PLN function), assuming the interest in in testing the model.}

\item{level}{Optional integer value the level to be used in obtaining the predictions. Level zero corresponds to the population predictions (default if \code{responses} is not provided) while level one (default) corresponds to predictions after evaluating the variational parameters for the new data.}

\item{type}{Scale used for the prediction. Either \code{"link"} (default, predicted positions in the latent space), \code{"response"} (predicted average counts, accounting for zero-inflation) or \code{"deflated"} (predicted average counts, not accounting for zero-inflation and using only the PLN part of the model).}

\item{...}{additional parameters for S3 compatibility. Not used}
}
\description{
Predict counts of a new sample
}
\details{
Note that \code{level = 1} can only be used if responses are provided,
as the variational parameters can't be estimated otherwise. In the absence of responses, \code{level} is ignored and the fitted values are returned

Note also that when \code{type = "response"} corresponds to predicting
values with \eqn{(1 - \pi)A}, where \eqn{A} is the average count in
the PLN part of the model and \eqn{\pi} the probability of zero-inflation,
whereas \code{type = "deflated"} corresponds to \eqn{A}.
}

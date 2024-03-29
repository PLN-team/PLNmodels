% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLNfit-S3methods.R
\name{coef.PLNfit}
\alias{coef.PLNfit}
\title{Extract model coefficients}
\usage{
\method{coef}{PLNfit}(object, type = c("main", "covariance"), ...)
}
\arguments{
\item{object}{an R6 object with class \code{\link{PLNfit}}}

\item{type}{type of parameter that should be extracted. Either "main" (default) for \deqn{B} or "covariance" for \deqn{\Sigma}}

\item{...}{additional parameters for S3 compatibility. Not used}
}
\value{
A matrix of coefficients extracted from the PLNfit model.
}
\description{
Extracts model coefficients from objects returned by \code{\link[=PLN]{PLN()}} and its variants
}
\examples{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
coef(myPLN) ## B
coef(myPLN, type = "covariance") ## Sigma
}
\seealso{
\code{\link[=sigma.PLNfit]{sigma.PLNfit()}}, \code{\link[=vcov.PLNfit]{vcov.PLNfit()}}, \code{\link[=standard_error.PLNfit]{standard_error.PLNfit()}}
}

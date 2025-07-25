% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLNfit-class.R
\name{PLNfit}
\alias{PLNfit}
\title{An R6 Class to represent a PLNfit in a standard, general framework}
\description{
The function \code{\link[=PLN]{PLN()}} fit a model which is an instance of a object with class \code{\link{PLNfit}}.
Objects produced by the functions \code{\link[=PLNnetwork]{PLNnetwork()}}, \code{\link[=PLNPCA]{PLNPCA()}}, \code{\link[=PLNmixture]{PLNmixture()}} and \code{\link[=PLNLDA]{PLNLDA()}} also enjoy the methods of \code{\link[=PLNfit]{PLNfit()}} by inheritance.

This class comes with a set of R6 methods, some of them being useful for the user and exported as S3 methods.
See the documentation for \code{\link[=coef]{coef()}}, \code{\link[=sigma]{sigma()}}, \code{\link[=predict]{predict()}}, \code{\link[=vcov]{vcov()}} and \code{\link[=standard_error]{standard_error()}}.

Fields are accessed via active binding and cannot be changed by the user.
}
\examples{
\dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- PLN(Abundance ~ 1, data = trichoptera)
class(myPLN)
print(myPLN)
}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{n}}{number of samples}

\item{\code{q}}{number of dimensions of the latent space}

\item{\code{p}}{number of species}

\item{\code{d}}{number of covariates}

\item{\code{nb_param}}{number of parameters in the current PLN model}

\item{\code{model_par}}{a list with the matrices of the model parameters: B (covariates), Sigma (covariance), Omega (precision matrix), plus some others depending on the variant)}

\item{\code{var_par}}{a list with the matrices of the variational parameters: M (means) and S2 (variances)}

\item{\code{optim_par}}{a list with parameters useful for monitoring the optimization}

\item{\code{latent}}{a matrix: values of the latent vector (Z in the model)}

\item{\code{latent_pos}}{a matrix: values of the latent position vector (Z) without covariates effects or offset}

\item{\code{fitted}}{a matrix: fitted values of the observations (A in the model)}

\item{\code{vcov_coef}}{matrix of sandwich estimator of the variance-covariance of B (need fixed -ie known- covariance at the moment)}

\item{\code{vcov_model}}{character: the model used for the residual covariance}

\item{\code{weights}}{observational weights}

\item{\code{loglik}}{(weighted) variational lower bound of the loglikelihood}

\item{\code{loglik_vec}}{element-wise variational lower bound of the loglikelihood}

\item{\code{AIC}}{variational lower bound of the AIC}

\item{\code{BIC}}{variational lower bound of the BIC}

\item{\code{entropy}}{Entropy of the variational distribution}

\item{\code{ICL}}{variational lower bound of the ICL}

\item{\code{R_squared}}{approximated goodness-of-fit criterion}

\item{\code{criteria}}{a vector with loglik, BIC, ICL and number of parameters}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-PLNfit-new}{\code{PLNfit$new()}}
\item \href{#method-PLNfit-update}{\code{PLNfit$update()}}
\item \href{#method-PLNfit-optimize}{\code{PLNfit$optimize()}}
\item \href{#method-PLNfit-optimize_vestep}{\code{PLNfit$optimize_vestep()}}
\item \href{#method-PLNfit-postTreatment}{\code{PLNfit$postTreatment()}}
\item \href{#method-PLNfit-predict}{\code{PLNfit$predict()}}
\item \href{#method-PLNfit-predict_cond}{\code{PLNfit$predict_cond()}}
\item \href{#method-PLNfit-show}{\code{PLNfit$show()}}
\item \href{#method-PLNfit-print}{\code{PLNfit$print()}}
\item \href{#method-PLNfit-clone}{\code{PLNfit$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-new"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-new}{}}}
\subsection{Method \code{new()}}{
Initialize a \code{\link{PLNfit}} model
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$new(responses, covariates, offsets, weights, formula, control)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{responses}}{the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{covariates}}{design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{offsets}}{offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{weights}}{an optional vector of observation weights to be used in the fitting process.}

\item{\code{formula}}{model formula used for fitting, extracted from the formula in the upper-level call}

\item{\code{control}}{a list-like structure for controlling the fit, see \code{\link[=PLN_param]{PLN_param()}}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-update"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-update}{}}}
\subsection{Method \code{update()}}{
Update a \code{\link{PLNfit}} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$update(
  B = NA,
  Sigma = NA,
  Omega = NA,
  M = NA,
  S = NA,
  Ji = NA,
  R2 = NA,
  Z = NA,
  A = NA,
  monitoring = NA
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{B}}{matrix of regression matrix}

\item{\code{Sigma}}{variance-covariance matrix of the latent variables}

\item{\code{Omega}}{precision matrix of the latent variables. Inverse of Sigma.}

\item{\code{M}}{matrix of variational parameters for the mean}

\item{\code{S}}{matrix of variational parameters for the variance}

\item{\code{Ji}}{vector of variational lower bounds of the log-likelihoods (one value per sample)}

\item{\code{R2}}{approximate R^2 goodness-of-fit criterion}

\item{\code{Z}}{matrix of latent vectors (includes covariates and offset effects)}

\item{\code{A}}{matrix of fitted values}

\item{\code{monitoring}}{a list with optimization monitoring quantities}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Update the current \code{\link{PLNfit}} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-optimize"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-optimize}{}}}
\subsection{Method \code{optimize()}}{
Call to the NLopt or TORCH optimizer and update of the relevant fields
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$optimize(responses, covariates, offsets, weights, config)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{responses}}{the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{covariates}}{design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{offsets}}{offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{weights}}{an optional vector of observation weights to be used in the fitting process.}

\item{\code{config}}{part of the \code{control} argument which configures the optimizer}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-optimize_vestep"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-optimize_vestep}{}}}
\subsection{Method \code{optimize_vestep()}}{
Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (Sigma, B). Intended to position new data in the latent space.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$optimize_vestep(
  covariates,
  offsets,
  responses,
  weights,
  B = self$model_par$B,
  Omega = self$model_par$Omega,
  control = PLN_param(backend = "nlopt")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{covariates}}{design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{offsets}}{offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{responses}}{the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{weights}}{an optional vector of observation weights to be used in the fitting process.}

\item{\code{B}}{Optional fixed value of the regression parameters}

\item{\code{Omega}}{precision matrix of the latent variables. Inverse of Sigma.}

\item{\code{control}}{a list-like structure for controlling the fit, see \code{\link[=PLN_param]{PLN_param()}}.}

\item{\code{Sigma}}{variance-covariance matrix of the latent variables}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with three components:
\itemize{
\item the matrix \code{M} of variational means,
\item the matrix \code{S2} of variational variances
\item the vector \code{log.lik} of (variational) log-likelihood of each new observation
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-postTreatment"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-postTreatment}{}}}
\subsection{Method \code{postTreatment()}}{
Update R2, fisher and std_err fields after optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$postTreatment(
  responses,
  covariates,
  offsets,
  weights = rep(1, nrow(responses)),
  config_post,
  config_optim,
  nullModel = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{responses}}{the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{covariates}}{design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{offsets}}{offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class}

\item{\code{weights}}{an optional vector of observation weights to be used in the fitting process.}

\item{\code{config_post}}{a list for controlling the post-treatments (optional bootstrap, jackknife, R2, etc.). See details}

\item{\code{config_optim}}{a list for controlling the optimization (optional bootstrap, jackknife, R2, etc.). See details}

\item{\code{nullModel}}{null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
The list of parameters \code{config} controls the post-treatment processing, with the following entries:
\itemize{
\item jackknife boolean indicating whether jackknife should be performed to evaluate bias and variance of the model parameters. Default is FALSE.
\item bootstrap integer indicating the number of bootstrap resamples generated to evaluate the variance of the model parameters. Default is 0 (inactivated).
\item variational_var boolean indicating whether variational Fisher information matrix should be computed to estimate the variance of the model parameters (highly underestimated). Default is FALSE.
\item sandwich_var boolean indicating whether sandwich estimator should be computed to estimate the variance of the model parameters (highly underestimated). Default is FALSE.
\item trace integer for verbosity. should be > 1 to see output in post-treatments
}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-predict"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-predict}{}}}
\subsection{Method \code{predict()}}{
Predict position, scores or observations of new data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$predict(
  newdata,
  responses = NULL,
  type = c("link", "response"),
  level = 1,
  envir = parent.frame()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{newdata}}{A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}

\item{\code{responses}}{Optional data frame containing the count of the observed variables (matching the names of the provided as data in the PLN function), assuming the interest is in testing the model.}

\item{\code{type}}{Scale used for the prediction. Either \code{link} (default, predicted positions in the latent space) or \code{response} (predicted counts).}

\item{\code{level}}{Optional integer value the level to be used in obtaining the predictions. Level zero corresponds to the population predictions (default if \code{responses} is not provided) while level one (default) corresponds to predictions after evaluating the variational parameters for the new data.}

\item{\code{envir}}{Environment in which the prediction is evaluated}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Note that \code{level = 1} can only be used if responses are provided,
as the variational parameters can't be estimated otherwise. In the absence of responses, \code{level} is ignored and the fitted values are returned
}

\subsection{Returns}{
A matrix with predictions scores or counts.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-predict_cond"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-predict_cond}{}}}
\subsection{Method \code{predict_cond()}}{
Predict position, scores or observations of new data, conditionally on the observation of a (set of) variables
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$predict_cond(
  newdata,
  cond_responses,
  type = c("link", "response"),
  var_par = FALSE,
  envir = parent.frame()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{newdata}}{a data frame containing the covariates of the sites where to predict}

\item{\code{cond_responses}}{a data frame containing the count of the observed variables (matching the names of the provided as data in the PLN function)}

\item{\code{type}}{Scale used for the prediction. Either \code{link} (default, predicted positions in the latent space) or \code{response} (predicted counts).}

\item{\code{var_par}}{Boolean. Should new estimations of the variational parameters of mean and variance be sent back, as attributes of the matrix of predictions. Default to \code{FALSE}.}

\item{\code{envir}}{Environment in which the prediction is evaluated}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A matrix with predictions scores or counts.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-show"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-show}{}}}
\subsection{Method \code{show()}}{
User friendly print method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$show(
  model = paste("A multivariate Poisson Lognormal fit with", self$vcov_model,
    "covariance model.\\n")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model}}{First line of the print output}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-print"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-print}{}}}
\subsection{Method \code{print()}}{
User friendly print method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PLNfit-clone"></a>}}
\if{latex}{\out{\hypertarget{method-PLNfit-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PLNfit$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

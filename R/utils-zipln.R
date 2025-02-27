.extract_terms_zi <- function(formula) {

  ## Check if a ZI specific formula has been provided
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
    zicovar <- TRUE
    ff_zi <-  ~. ; ff_zi[[3]]  <- formula[[3]][[3]] ; ff_zi[[2]]  <- NULL
    ff_pln <- ~. ; ff_pln[[3]] <- formula[[3]][[2]] ; ff_pln[[2]] <- NULL
    tt_zi  <- terms(ff_zi)  ; attr(tt_zi , "offset") <- NULL
    tt_pln <- terms(ff_pln) ; attr(tt_pln, "offset") <- NULL
    formula[[3]][1] <- call("+")
  } else {
    ff_pln <- formula
    ff_zi <- NULL
    tt_pln <- terms(ff_pln) ; attr(tt_pln, "offset") <- NULL
    tt_zi  <- NULL
    zicovar <- FALSE
  }

  list(ZI = tt_zi, PLN = tt_pln, formula = formula, zicovar = zicovar)
}

#' @importFrom stats .getXlevels
extract_model_zi <- function(call, envir) {

  ## create the call for the model frame
  call_args  <- call[match(c("formula", "data", "subset", "weights"), names(call), 0L)]
  call_args <- c(as.list(call_args), list(xlev = attr(call$formula, "xlevels"), na.action = NULL))

  ## Extract terms for ZI and PLN components
  terms <- .extract_terms_zi(as.formula(eval(call$formula, envir = envir)))
  ## eval the call in the parent environment with adjustement due to ZI terms
  call_args$formula <- terms$formula
  frame <- do.call(stats::model.frame, call_args, envir = envir)

  ## Save level for predict function
  xlevels <- list(PLN = .getXlevels(terms$PLN, frame))
  if (!is.null(terms$ZI)) xlevels$ZI = .getXlevels(terms$ZI, frame)
  if (!is.null(xlevels$PLN)) attr(call$formula, "xlevels") <- xlevels

  ## Create the set of matrices to fit the PLN model
  X  <- model.matrix(terms$PLN, frame, xlev = xlevels$PLN)
  if (terms$zicovar) X0 <- model.matrix(terms$ZI, frame, xlev = xlevels$ZI) else X0 <- matrix(NA,0,0)

  Y <- model.response(frame)
  ## model.response oversimplifies into a numeric when a single variable is involved
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = length(Y), ncol = 1)
  if (ncol(Y) == 1 & is.null(colnames(Y))) colnames(Y) <- "Y"

  # Offsets are only considered for the PLN component
  O <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))
  if (is.vector(O)) O <- O %o% rep(1, ncol(Y))

  # Model weights
  w <- model.weights(frame)
  if (is.null(w)) {
    w <- rep(1.0, nrow(Y))
  } else {
    stopifnot(all(w > 0) && length(w) == nrow(Y))
  }

  list(Y = Y, X = X, X0 = X0, O = O, w = w, formula = call$formula, zicovar = terms$zicovar)
}

# Test convergence for a named list of parameters
# oldp, newp: named list of parameters
# xtol_rel: double ; negative or NULL = disabled
# xtol_abs: double ; negative or NULL = disabled
# Returns boolean
parameter_list_converged <- function(oldp, newp, xtol_abs = NULL, xtol_rel = NULL) {
  # Strategy is to compare each pair of list elements with matching names.
  stopifnot(is.list(oldp), is.list(newp))
  oldp <- oldp[order(names(oldp))]
  newp <- newp[order(names(newp))]
  stopifnot(all(names(oldp) == names(newp)))

  # Check convergence with xtol_rel if enabled
  if(is.double(xtol_rel) && xtol_rel > 0) {
    if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_rel * abs(o)) }, oldp, newp))) {
      return(TRUE)
    }
  }

  # Check convergence with xtol_abs (homogeneous) if enabled
  if(is.double(xtol_abs) && xtol_abs > 0) {
    if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_abs) }, oldp, newp))) {
      return(TRUE)
    }
  }

  # If no criteria has triggered, indicate no convergence
  FALSE
}

#' #' @importFrom glmnet glmnet
#' optim_zipln_B <- function(M, X, Omega, config) {
#'
#'   if(config$lambda > 0) {
#'     if (!is.null(config$ind_intercept)) {
#'       m_bar <- colMeans(M)
#'       x_bar <- colMeans(X[, -config$ind_intercept])
#'       X <- scale(X[, -config$ind_intercept], x_bar, FALSE)
#'       M <- scale(M, m_bar, FALSE)
#'     }
#'     p <- ncol(M); d <- ncol(X)
#'     if (d > 0) {
#'       Omega12 <- chol(Omega)
#'       y <- as.vector(M %*% t(Omega12))
#'       x <- kronecker(Omega12, X)
#'       glmnet_out <- glmnet(x, y, lambda = config$lambda, intercept = FALSE, standardize = FALSE)
#'       B <- matrix(as.numeric(glmnet_out$beta), nrow = d, ncol = p)
#'     } else {
#'       B <- matrix(0, nrow = d, ncol = p)
#'     }
#'
#'     if (!is.null(config$ind_intercept)) {
#'       mu0 <- m_bar - as.vector(crossprod(B, x_bar))
#'       B <- rbind(mu0, B)
#'     }
#'
#'   } else {
#'     B <- optim_zipln_B_dense(M, X)
#'   }
#'   B
#' }
#'

#' #' Helper function for ZIPLN initialization.
#' #'
#' #' @description
#' #' Barebone function to compute starting points for B, M and S when fitting a PLN. Mostly intended for internal use.
#' #'
#' #' @param Y Response count matrix
#' #' @param X Design matrix for the count (PLN) component
#' #' @param X0 Design matrix for the zero (Bernoulli) component
#' #' @param O Offset matrix (in log-scale)
#' #' @param w Weight vector (defaults to 1)
#' #' @param s Scale parameter for S (defaults to 0.1)
#' #' @return a named list of starting values for model parameter B and variational parameters M and S used in the iterative optimization algorithm of [PLN()]
#' #'
#' #' @details The default strategy to estimate B and M is to fit a linear model with covariates `X` to the response count matrix (after adding a pseudocount of 1, scaling by the offset and taking the log). The regression matrix is used to initialize `B` and the residuals to initialize `M`. `S` is initialized as a constant conformable matrix with value `s`.
#' #'
#' #' @rdname compute_ZIPLN_starting_point
#' #' @examples
#' #' \dontrun{
#' #' data(barents)
#' #' Y <- barents$Abundance
#' #' X  <- model.matrix(Abundance ~ 1 + Temperature, data = barents)
#' #' X0 <- model.matrix(Abundance ~ 0 + Longitude, data = barents)
#' #' O <- log(barents$Offset)
#' #' w <- rep(1, nrow(Y))
#' #' compute_ZIPLN_starting_point(Y, X, O, w)
#' #' }
#' #'
#' #' @importFrom stats lm.fit
#' #' @export
#' compute_ZIPLN_starting_point <- function(Y, X, X0, O, w, s = 0.1) {
#'   # Y = responses, X = covariates, O = offsets (in log scale), w = weights
#'   n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
#'   fits <- lm.fit(w * X, w * log((1 + Y)/exp(O)))
#'   list(B = matrix(fits$coefficients, d, p),
#'        M = matrix(fits$residuals, n, p),
#'        S = matrix(s, n, p))
#' }

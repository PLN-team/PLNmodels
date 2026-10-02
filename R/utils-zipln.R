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

  ## The formula may be passed as a variable, in which case call$formula is a symbol
  formula_obj <- as.formula(eval(call$formula, envir = envir))

  ## Extract terms for ZI and PLN components
  terms <- .extract_terms_zi(formula_obj)
  ## eval the call in the parent environment with adjustement due to ZI terms
  call_args$formula <- terms$formula
  frame <- do.call(stats::model.frame, call_args, envir = envir)

  ## Save level for predict function
  xlevels <- list(PLN = .getXlevels(terms$PLN, frame))
  if (!is.null(terms$ZI)) xlevels$ZI = .getXlevels(terms$ZI, frame)
  if (!is.null(xlevels$PLN)) attr(formula_obj, "xlevels") <- xlevels

  ## Create the set of matrices to fit the PLN model
  X  <- model.matrix(terms$PLN, frame, xlev = xlevels$PLN)
  if (terms$zicovar) X0 <- model.matrix(terms$ZI, frame, xlev = xlevels$ZI) else X0 <- matrix(NA,0,0)

  ## Offsets are only considered for the PLN component
  YOw <- extract_response_offset_weights(frame)

  list(Y = YOw$Y, X = X, X0 = X0, O = YOw$O, w = YOw$w, formula = formula_obj, zicovar = terms$zicovar)
}

# Test convergence between two named lists of parameters, element by element
# (a tolerance which is NULL or not positive is disabled)
parameter_list_converged <- function(oldp, newp, xtol_abs = NULL, xtol_rel = NULL) {
  stopifnot(is.list(oldp), is.list(newp))
  oldp <- oldp[order(names(oldp))]
  newp <- newp[order(names(newp))]
  stopifnot(all(names(oldp) == names(newp)))

  if(is.double(xtol_rel) && xtol_rel > 0) {
    if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_rel * abs(o)) }, oldp, newp))) {
      return(TRUE)
    }
  }

  if(is.double(xtol_abs) && xtol_abs > 0) {
    if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_abs) }, oldp, newp))) {
      return(TRUE)
    }
  }

  FALSE
}

# C++ VE step of ZIPLN for a backend and a covariance structure (`suffix`)
zipln_MS_fn <- function(backend, suffix) {
  if (backend == "builtin") {
    switch(suffix,
      full      = builtin_optimize_vestep_zipln_full,
      diagonal  = builtin_optimize_vestep_zipln_diagonal,
      spherical = builtin_optimize_vestep_zipln_spherical,
      fixed     = builtin_optimize_vestep_zipln_fixed)
  } else {
    switch(suffix,
      full      = nlopt_optimize_vestep_zipln_full,
      diagonal  = nlopt_optimize_vestep_zipln_diagonal,
      spherical = nlopt_optimize_vestep_zipln_spherical,
      fixed     = nlopt_optimize_vestep_zipln_fixed)
  }
}

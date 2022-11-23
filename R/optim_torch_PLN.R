#' @import torch
optimize_torch_PLN <- function(responses, covariates, offsets, weights, init_parameters, control) {

  ## problem dimensions
  p <- ncol(responses)

  ## conversion to torch tensor
  X   <- torch_tensor(covariates)
  Y   <- torch_tensor(responses)
  O   <- torch_tensor(offsets)
  w   <- torch_tensor(weights)
  w_bar <- sum(w)

  get_objective <- function() {
    S2 <- S * S
    Z <- O + M + torch_matmul(X, Theta)
    log_det_Sigma <- switch(control$covariance,
      "spherical" = p * log(sum(torch_matmul(w, M * M + S2)) / (w_bar * p)),
      "diagonal"  = sum(torch_log(torch_matmul(w, M * M + S2) / w_bar)),
      { # default value
      Mw <- torch_matmul(torch_diag(torch_sqrt(w)), M)
      Sigma <- (torch_matmul(torch_transpose(Mw, 2, 1), Mw) + torch_diag(torch_matmul(w, S2))) / w_bar
      torch_logdet(Sigma)
      })
    neg_ELBO <- .5 * w_bar * log_det_Sigma -
      sum(torch_matmul(w , Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
    neg_ELBO
  }

  get_Sigma <- function(M, S) {
    switch(control$covariance,
           "spherical" = torch_eye(p) * sum(torch_matmul(w, M * M + S*S)) / (w_bar * p),
           "diagonal"  = torch_diag(torch_matmul(w, M * M + S*S) / w_bar),
           "full"      = {
             Mw <- torch_matmul(torch_diag(torch_sqrt(w)), M)
             (torch_matmul(torch_transpose(Mw, 2, 1), Mw) + torch_diag(torch_matmul(w, S*S))) / w_bar
           }
    )
  }

  ## Initialization
  Theta <- torch_tensor(init_parameters$Theta, requires_grad = TRUE)
  M     <- torch_tensor(init_parameters$M    , requires_grad = TRUE)
  S     <- torch_tensor(init_parameters$S    , requires_grad = TRUE)
  optimizer <- optim_rprop(c(Theta = Theta, M = M, S = S), lr = control$learning_rate)
  Theta_old <- as.numeric(optimizer$param_groups[[1]]$params$Theta)

  ## Optimization loop
  message <- "failure"
  objective <- double(length = control$maxeval + 1)
  for (iterate in seq.int(control$maxeval)) {

    ## Optimization
    optimizer$zero_grad()   # reinitialize gradients
    loss <- get_objective() # compute current ELBO
    loss$backward()         # backward propagation
    optimizer$step()        # optimization

    ## assess convergence
    objective[iterate + 1] <- loss$item()
    Theta_new <- as.numeric(optimizer$param_groups[[1]]$params$Theta)
    delta_f   <- abs(objective[iterate] - objective[iterate + 1]) / abs(objective[iterate + 1])
    delta_x   <- sum(abs(Theta_old - Theta_new))/sum(abs(Theta_new))
    Theta_old <- Theta_new

    ## Display progress
    if (control$trace >  1 && (iterate %% 50 == 0))
      cat('\niteration: ', iterate, 'objective', objective[iterate + 1],
          'delta_f'  , round(delta_f, 6), 'delta_x', round(delta_x, 6))

    ## Check for convergence
    if (delta_f < control$ftol_rel | delta_x < control$xtol_rel) {
      objective <- objective[1:iterate + 1]
      message <- "converged"
      break
    }
  }

  Sigma <- get_Sigma(M, S)
  Omega <- torch::torch_inverse(Sigma)
  S2 <- torch::torch_multiply(S,S)
  Z  <- O + M + torch::torch_matmul(X, Theta)
  A  <- torch::torch_exp(Z + S2/2)
  KY <- rowSums(.logfactorial(responses))

  Ji <- as.numeric(
    .5 * torch_logdet(Omega) +
      torch_sum(Y * Z - A + .5 * torch_log(S2), dim = 2) -
      .5 * torch_sum(torch_matmul(M, Omega) * M + S2 * torch_diag(Omega), dim = 2) +
      .5 * p - KY
  )
  attr(Ji, "weights") <- weights

  out <- list(
    Theta      = Theta,
    Sigma      = Sigma,
    M          = M,
    S          = as.matrix(S),
    Z          = Z,
    A          = A,
    Ji         = Ji,
    monitoring = list(
      objective  = objective,
      iterations = iterate,
      message    = message)
  )
  out
}

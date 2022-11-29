#' @import torch
optimize_torch_PLN <- function(Y, X, O, w, init_parameters, configuration) {

  ## problem dimensions
  p <- ncol(Y)

  ## conversion to torch tensor
  KY  <- rowSums(.logfactorial(Y))
  X   <- torch_tensor(X)
  Y   <- torch_tensor(Y)
  O   <- torch_tensor(O)
  w   <- torch_tensor(w)
  w_bar <- sum(w)

  get_objective <- function() {
    S2 <- S * S
    Z <- O + M + torch_matmul(X, Theta)
    log_det_Sigma <- switch(configuration$covariance,
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

  ## Initialization
  Theta <- torch_tensor(t(init_parameters$Theta), requires_grad = TRUE)
  M     <- torch_tensor(init_parameters$M       , requires_grad = TRUE)
  S     <- torch_tensor(init_parameters$S       , requires_grad = TRUE)
  optimizer <- optim_rprop(c(Theta = Theta, M = M, S = S), lr = configuration$learning_rate)
  Theta_old <- as.numeric(optimizer$param_groups[[1]]$params$Theta)

  ## Optimization loop
  status <- 5
  objective <- double(length = configuration$maxeval + 1)
  for (iterate in seq.int(configuration$maxeval)) {

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
    if (configuration$trace >  1 && (iterate %% 50 == 0))
      cat('\niteration: ', iterate, 'objective', objective[iterate + 1],
          'delta_f'  , round(delta_f, 6), 'delta_x', round(delta_x, 6))

    ## Check for convergence
    if (delta_f < configuration$ftol_rel | delta_x < configuration$xtol_rel) {
      objective <- objective[1:iterate + 1]
      if (delta_f < configuration$ftol_rel) status <- 3
      if (delta_x < configuration$xtol_rel) status <- 4
      break
    }
  }

  Z  <- O + M + torch::torch_matmul(X, Theta)
  A  <- torch::torch_exp(Z + S*S/2)

  out <- list(
    Theta      = t(as.matrix(Theta)),
    M          = as.matrix(M),
    S          = as.matrix(S),
    Z          = as.matrix(Z),
    A          = as.matrix(A),
    objective  = objective,
    iterations = iterate,
    status     = status
  )
  out
}

#' @import torch
optimize_PLNPCA <- function(responses, covariates, offsets, weights, init_parameters, control) {

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
    log_det_Sigma <- switch(configuration$covariance,
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
    switch(configuration$covariance,
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
  optimizer <- optim_rprop(c(Theta = Theta, M = M, S = S), lr = configuration$learning_rate)
  Theta_old <- as.numeric(optimizer$param_groups[[1]]$params$Theta)

  ## Optimization loop
  message <- "failure"
  objective <- double(length = configuration$maxeval + 1)
  for (iterate in seq.int(configuration$maxeval)) {

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
    if (configuration$trace >  1 && (iterate %% 50 == 0))
      cat('\niteration: ', iterate, 'objective', objective[iterate + 1],
          'delta_f'  , round(delta_f, 6), 'delta_x', round(delta_x, 6))

    ## Check for convergence
    if (delta_f < configuration$ftol_rel | delta_x < configuration$xtol_rel) {
      objective <- objective[1:iterate + 1]
      message <- "converged"
      break
    }
  }

  # Sigma <- get_Sigma(M, S)
  # Omega <- torch::torch_inverse(Sigma)
  # S2 <- torch::torch_multiply(S,S)
  # Z  <- O + M + torch::torch_matmul(X, Theta)
  # A  <- torch::torch_exp(Z + S2/2)
  # KY <- rowSums(.logfactorial(responses))
  #
  # Ji <- as.numeric(
  #   .5 * torch_logdet(Omega) +
  #     torch_sum(Y * Z - A + .5 * torch_log(S2), dim = 2) -
  #     .5 * torch_sum(torch_matmul(M, Omega) * M + S2 * torch_diag(Omega), dim = 2) +
  #     .5 * p - KY
  # )
  # attr(Ji, "weights") <- w

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

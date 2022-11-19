#' @import torch
optimize_torch_PLN <- function(responses, covariates, offsets, weights, init_parameters, control) {

  ## problem dimensions
  n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates)

  ## conversion to torch tensor
  X   <- torch_tensor(covariates)
  Y   <- torch_tensor(responses)
  O   <- torch_tensor(offsets)
  w   <- torch_tensor(weights)
  W   <- torch_diag(w)
  w_bar <- sum(w)

  get_objective <-
    switch(control$covariance,
           "spherical" = function(Theta, M, S) {
             S2 <- S*S
             Z <- O + M + torch_matmul(X, Theta)
             sigma2 <- sum(torch_matmul(w, M * M + S2)) / (w_bar * p)
             .5 * p * w_bar * log(sigma2) - sum(torch_matmul(w , Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
           },
           "diagonal"  = function(Theta, M, S) {
             S2 <- S*S
             Z <- O + M + torch_matmul(X, Theta)
             D <-  torch_matmul(w, M * M + S2) / w_bar;
             .5 * w_bar * sum(torch_log(D)) - sum(torch_matmul(w , Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
           },
           "full"      = function(Theta, M, S) {
             S2 <- S*S
             Z <- O + M + torch_matmul(X, Theta)
             Mw <- torch_matmul(torch_diag(torch_sqrt(w)), M)
             Sigma <- (torch_matmul(torch_transpose(Mw, 2, 1), Mw) + torch_diag(torch_matmul(w, S2))) / w_bar
             .5 * w_bar * torch_logdet(Sigma) - sum(torch_matmul(w, Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
           },
           "fixed"      = function(Theta, M, S) {
             S2 <- S*S
             Z <- O + M + torch_matmul(X, Theta)
             Mw <- torch_matmul(torch_diag(torch_sqrt(w)), M)
             Sigma <- (torch_matmul(torch_transpose(Mw, 2, 1), Mw) + torch_diag(torch_matmul(w, S2))) / w_bar
             .5 * w_bar * torch_logdet(Sigma) - sum(torch_matmul(w, Y * Z - torch_exp(Z + .5 * S2) + .5 * torch_log(S2)))
           }

    )

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
  Theta <- torch::torch_tensor(init_parameters$Theta, requires_grad = TRUE)
  M     <- torch::torch_tensor(init_parameters$M, requires_grad = TRUE)
  S     <- torch::torch_tensor(init_parameters$S, requires_grad = TRUE)
  objective <- double(length = control$maxeval + 1)
  optimizer <- torch::optim_rprop(c(Theta, M, S), lr = control$learning_rate)
  message <- "failure"

  ## Optimization loop
  for (iterate in seq.int(control$maxeval)) {
    ## Optimization
    optimizer$zero_grad()                # reinitialize gradients
    loss <- get_objective(Theta, M, S)   # compute current ELBO
    loss$backward()                      # backward propagation and optimization
    optimizer$step()
    ## Display progress
    objective[iterate + 1] <- loss$item()
    if (control$trace >  1 && (iterate %% 50 == 0))
      cat(' \n iteration : ', iterate, 'Objective', objective[iterate + 1])
    ## Check for convergence
    if (abs(objective[iterate] - objective[iterate + 1]) / abs(objective[iterate + 1]) < control$ftol_rel) {
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
    S2         = S2,
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

library(torch)
library(R6)

log_stirling <- function(n_){
  n_ <- n_+ (n_==0)
  torch_log(torch_sqrt(2*pi*n_)) + n_*log(n_/exp(1))
}

PLN <-
  R6Class("PLN",
          public = list(
            Y = NULL,
            O = NULL,
            X = NULL,
            n = NULL,
            p = NULL,
            d = NULL,
            M = NULL,
            S = NULL,
            A = NULL,
            B = NULL,
            Sigma = NULL,
            Omega = NULL,
            loglik = NULL,
            ELBO_list = NULL,

            ## Constructor
            initialize = function(Y, O, X){
              self$Y <- torch_tensor(Y)
              self$O <- torch_tensor(O)
              self$X <- torch_tensor(X)
              self$n <- nrow(Y)
              self$p <- ncol(Y)
              self$d <- ncol(X)
              ## Variational parameters
              self$M <- torch_zeros(self$n, self$p, requires_grad = TRUE)
              self$S <- torch_ones(self$n , self$p, requires_grad = TRUE)
              ## Model parameters
              self$B <- torch_full(c(self$d, self$p), -8.0, requires_grad = TRUE)
              self$Sigma <- torch_eye(self$p)
              self$Omega <- torch_eye(self$p)
              ## Monitoring
              self$ELBO_list <- c()
            },

            get_Sigma = function(M, S){
              1/self$n * (torch_mm(torch_t(M),M) + torch_diag(torch_sum(S**2, dim = 1)))
            },

            get_ELBO = function(B, Phi, M, S, Sigma){

              S2 <- torch_square(S)

              M_shift <- torch_roll(M)
              M_diff - M - torch_mm(Phi, M_shift)
              mu <- self$O + torch_mm(self$X[index], B)
              mu_eps <- mu - torch_mm(Phi, mu)
              Omega_eps <- torch_inverse(Sigma - torch_mm(torch_mm(Phi, Sigma), torch_transpose(Phi)))

              A <- torch_exp(self$O + M + .5 * S2)

              res <-
                torch_sum(self$Y * (self$O + M) - A + .5 * torch_log(S2))


              - self$n/2 * torch_logdet(Sigma) -
                   +

                 +

                .5 * torch_omega_squared_norm(torch_index_select(M, 2, 1:self$p) - mu, params$Omega) +
                .5 * torch_frobenius(torch_index_select(S2, 2, 1:self$p), torch_diag(params$Omega), dim = 2) -
                .5 * sum(data$w[index] - 1) * torch_logdet(params$Omega_eps[index]) +
                .5 * torch_omega_squared_norm(M_diff - mu_eps, Omega_eps) +
                .5 * torch_frobenius(torch_index_select(S2, 2, 2:data$w[index]), torch_diag(Omega_eps))

              + .5 * torch_log(S2)
              res



              S2 <- torch_square(S)
              XB <- torch_mm(self$X, B)

              A  <- torch_exp(self$O + M + XB + S2/2)
              self$n/2 * torch_logdet(Omega) -
                torch_sum(A - self$Y * (self$O + M + XB) - .5 * torch_log(S2))
            },

            get_loglik = function(B, M, S, Omega) {
              S2 <- S**2
              XB <- torch_mm(self$X, B)
              A  <- torch_exp(self$O + M + XB + S2/2)
              J <- self$n/2 * torch_logdet(Omega) +
                .5 * self$n * self$p - torch_sum(log_stirling(self$Y)) -
                torch_sum(A - self$Y * (self$O + M + XB) - .5 * torch_log(S2)) -
                .5 * torch_sum(torch_mm(M, Omega) * M + S2 * torch_diag(Omega))
              J
            },

            fit = function(N_iter, lr, tol = 1e-8, verbose = FALSE){
              self$ELBO_list <- double(length = N_iter)
              optimizer <- optim_rprop(c(self$B, self$M, self$S), lr = lr)
              objective0 <- Inf
              for (i in 1:N_iter){
                ## reinitialize gradients
                optimizer$zero_grad()

                ## compute current ELBO
                loss <- - self$get_ELBO(self$B, self$M, self$S, self$Omega)

                ## backward propagation and optimization
                loss$backward()
                optimizer$step()

                ## update parameters with close form
                self$Sigma <- self$get_Sigma(self$M, self$S)
                self$Omega <- torch_inverse(self$Sigma)

                objective <- -loss$item()
                if(verbose && (i %% 50 == 0)){
                  pr('i : ', i )
                  pr('ELBO', objective)
                }
                self$ELBO_list[i] <- objective
                if (abs(objective0 - objective)/abs(objective) < tol) {
                  self$ELBO_list <- self$ELBO_list[1:i]
                  break
                } else {
                  objective0 <- objective
                }
              }
              self$loglik <- self$get_loglik(self$B, self$M, self$S, self$Omega)
            },

            plotLogNegElbo = function(from = 10){
              plot(log(-self$ELBO_list[from:length(self$ELBO_list) ]), type = "l")
            }
          )
  )

Y <- oaks$Abundance
X <- cbind(rep(1, nrow(Y)))
O <- log(oaks$Offset)
myPLN <- PLN$new(Y = Y, O = O, X = X)

system.time(myPLN$fit(1000, lr = 0.05, tol = 1e-9))

plot(-myPLN$ELBO_list, type="l")



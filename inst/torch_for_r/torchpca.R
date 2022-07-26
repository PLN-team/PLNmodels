library(torch)
library(R6)




ELBO_PCA <- function(Y, O, covariates, M, S, C, Theta){
  ## compute the ELBO with a PCA parametrization'''
  n = Y$shape[1]
  q = C$shape[2]
  # Store some variables that will need to be computed twice
  A = O + torch_mm(covariates, Theta) + torch_mm(M, C$t())
  SrondS = torch_multiply(S, S)
  
  YA = torch_sum(torch_multiply(Y, A))
  moinsexpAplusSrondSCCT = torch_sum(-torch_exp(A + 1 / 2 *
                                                  torch_mm(SrondS, torch_multiply(C, C)$t())))
  moinslogSrondS = 1 / 2 * torch_sum(torch_log(SrondS))
  MMplusSrondS = torch_sum(-1 / 2 * (torch_multiply(M, M) + torch_multiply(S, S)))
  log_stirlingY = torch_sum(log_stirling(Y))
  return(YA + moinsexpAplusSrondSCCT + moinslogSrondS + MMplusSrondS - log_stirlingY + n * q / 2)
}

VEM_PLNPCA <- R6Class("VEM_PLNPCA", 
                      public = list(
                        Y = NULL, 
                        O = NULL,
                        covariates = NULL, 
                        p = NULL, 
                        q = NULL,
                        n = NULL, 
                        d = NULL, 
                        M = NULL, 
                        S = NULL, 
                        A = NULL, 
                        C = NULL,
                        Sigma = NULL, 
                        Theta = NULL,
                        good_init = NULL,
                        fitted = NULL, 
                        ELBO_list = NULL, 
                        initialize = function(Y,O,covariates,q, good_init = TRUE){
                          self$Y <- Y
                          self$O <- O
                          self$covariates <- covariates
                          self$q = q
                          self$good_init = good_init 
                          self$p <- Y$shape[2]
                          self$n <- Y$shape[1]
                          self$d <- covariates$shape[2]
                          ## Variational parameters
                          self$M <- torch_zeros(self$n, self$q, requires_grad = TRUE)
                          self$S <- torch_ones(self$n, self$q, requires_grad = TRUE)
                          ## Model parameters 
                          print('Initialization ...')
                          if (self$good_init){
                            self$Theta <- Poisson_reg(Y,O,covariates)$detach()$clone()$requires_grad_(TRUE)
                            self$C <- init_C(Y,O,covariates, self$Theta,q)$detach()$clone()$requires_grad_(TRUE)
                          }
                          else{
                            self$Theta <- torch_zeros(self$d, self$p, requires_grad = TRUE)
                            self$C <- torch_randn(self$p, self$q, requires_grad = TRUE)
                          }
                          print('Initialization finished.')
                          self$Sigma <- torch_eye(self$p)
                          self$fitted = FALSE
                          self$ELBO_list = c()
                        },
                        getSigma = function(){
                          return(torch_mm(self$C, self$C$t()))
                        },
                        fit = function(N_iter, lr, verbose = FALSE){
                          optimizer = optim_rprop(c(self$Theta, self$C, self$M, self$S), lr = lr)
                          for (i in 1:N_iter){
                            optimizer$zero_grad()
                            loss = - ELBO_PCA(self$Y, self$O, self$covariates, self$M, self$S, self$C, self$Theta)
                            loss$backward()
                            optimizer$step()
                            if(verbose && (i%%50 == 0)){
                              print('i :')
                              print(i)
                              print('ELBO :')
                              print(-loss$item()/(self$n))
                            }
                            self$ELBO_list = c(self$ELBO_list, -loss$item()/(self$n))
                          }
                          self$fitted <- TRUE 
                        }, 
                        plotLogNegElbo = function(from = 10){
                          plot(log(-self$ELBO_list[from:length(self$ELBO_list) ]))
                        }
                      )
)




sample_PLN <- function(C,Theta,O,covariates, B_zero = None, ZI = FALSE){
  #Sample Poisson log Normal variables. If ZI is True, then the model will
  #be zero inflated.
  #The sample size n is the the first size of O, the number p of variables
  #considered is the second size of O. The number d of covariates considered
  #is the first size of beta.
  
  #Args:
  #  C: torch.tensor of size (p,q). The matrix c of the PLN model
  #  Theta: torch.tensor of size (d,p).
  #  0: torch.tensor. Offset, size (n,p)
  #  covariates : torch.tensor. Covariates, size (n,d)
  #  B_zero: torch.tensor of size (d,p), optional. If ZI is True,
  #      it will raise an error if you don't set a value. Default is None.
  #      ZI: Bool, optional. If True, the model will be Zero Inflated. Default is False.
  #  Returns :
  #      Y: torch.tensor of size (n,p), the count variables.
  #      Z: torch.tensor of size (n,p), the gaussian latent variables.
  #      ksi: torch.tensor of size (n,p), the bernoulli latent variables.
  
  n_ <- dim(covariates)[1] 
  q_ <- dim(C)[2]
  
  XB = torch_matmul(covariates$unsqueeze(2),Theta$unsqueeze(1))$squeeze()
  
  Z <- torch_matmul(torch_randn(n,q), torch_transpose(C, 2,1)) + XB + O
  parameter = torch_exp(Z)
  if (ZI == TRUE){
    ZI_covariates <- as.matrix(covariates,Theta)
    ksi <- torch_tensor(matrix(rbinom(n = n*p,prob = 1/(1+ exp(-as.matrix(ZI_covariates))), size = 1), nrow = n, ncol = p))
  }
  else{
    ksi <- 0
  }
  Y = torch_tensor((1-ksi)*matrix(rpois(n = n*p, lambda = as.matrix(parameter)), nrow = n, ncol = p))
  return(Y)
}
log_stirling <- function(n_){
  n_ <- n_+ (n_==0)
  return(torch_log(torch_sqrt(2*pi*n_)) + n_*log(n_/exp(1))) 
}




d = 1L
n = 200L
p = 2000L
q = 10L

### Sampling some data ###
# O <-  torch_tensor(matrix(0,nrow = n, ncol = p))
# covariates <- torch_tensor(matrix(rnorm(n*d),nrow = n, ncol = d))
# true_Theta <- torch_tensor(matrix(rnorm(d*p),nrow = d, ncol = p))
# true_C <- torch_tensor(matrix(rnorm(p*q), nrow = p, ncol = q) )/3
# true_Sigma <- torch_matmul(true_C,torch_transpose(true_C, 2,1))
# true_Theta <- torch_tensor(matrix(rnorm(d*p),nrow = d, ncol = p))/2
# Y <- sample_PLN(true_C,true_Theta,O,covariates)



#plnpca = VEM_PLNPCA$new(Y,O,covariates, q, good_init = FALSE)
#plnpca$fit(300,0.1,verbose = TRUE)
#plnpca$plotLogNegElbo()




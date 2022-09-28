# Attempts of PLN network, with prescribed network backbone
# Exemple JC :

library(PLNmodels)
set.seed(1)

# Data
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance[1:20, 1:5], trichoptera$Covariate[1:20, ])
# trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
p <- ncol(trichoptera$Abundance)

# Metaweb
# G <- matrix(rbinom(p^2, 1, 1/p), p, p); G <- G + t(G); G[which(G>1)] <- 1
G <- matrix(c(0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1,
              0, 0, 0, 0, 0, 1, 0, 1, 0, 0), p, p)
colnames(G) <- rownames(G) <- colnames(trichoptera$Abundance)

# Penalties
penalties <- 10^seq(6, -6, length.out=53); penaltyMin <- min(penalties)/100
list_W <- lapply(penalties, function(rho){penaltyMin*G + rho*(1-G)})
WedgeList <- c()
for(j in 1:(p-1)){for(k in (j+1):p){
  if(G[j, k]==1){WedgeList <- c(WedgeList, paste0(colnames(G)[j], '|', colnames(G)[k]))}
}}

# Fit
myPLN <- PLNnetwork(Abundance ~ 1, data=trichoptera,
                    control_init=list(penalty_weights = list_W))
plot(myPLN)

# Actual penalties
path <- myPLN$coefficient_path()
penalties
unique(path$Penalty)

# Coef path
plot(penalties, rep(0, length(myPLN$penalties)), log='x', ylim=range(path$Coeff), col=0)
for(e in unique(path$Edge)){
  # lines(path$Penalty[which(path$Edge==e)], path$Coeff[which(path$Edge==e)], col=1+1*(e%in%WedgeList))
  lines(penalties, path$Coeff[which(path$Edge==e)], col=1+1*(e%in%WedgeList))
}

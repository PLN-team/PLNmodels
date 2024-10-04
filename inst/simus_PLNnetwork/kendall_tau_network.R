library(dplyr)
library(PLNmodels)

###############################################################################
## Weighted Kendall's tau
# Cf Shi98-StatProbLetters, Vig14-ArXiv
F_WeightedTau <- function(X, Y, weight.method='unif'){
  R = rank(X); S = rank(Y)
  n = length(R); N = n*(n-1)/2
  deltaR = deltaS = deltaR2 = deltaS2 = rep(0, N)
  weight = rep(1, N)
  ij = 0
  for (i in 1:(n-1)) {
    deltaR[ij+(1:(n-i))] <- sign(R[i]-R[(i+1):n])
    deltaS[ij+(1:(n-i))] <- sign(S[i]-S[(i+1):n])
    ij <- ij + (n-i)
  }
  if(weight.method == 'inverse'){# wij = vi vj; vi = 1/(1+sqrt(R[i]S[i]))
    v = 1/(1+sqrt(R*S)); vprod = v%o%v; weight = vprod[lower.tri(vprod, diag = FALSE)]
  }
  tau = sum(deltaR*deltaS*weight) / sum(weight)
  sigma = sqrt(mean(mean(rowMeans(vprod)^2)))
  stat = 1.5*sqrt(n)*tau*mean(weight) / sigma
  list(tau=tau, stat=stat, pval=pnorm(abs(stat), lower.tail=FALSE), v=v)
}

#########################################
## Extracting the maximal penalty for each edge from PLnetwork output
# Cf Shi98-StatProbLetters, Vig14-ArXiv
F_GetEdgePenalty <- function(PLNnet){
  PLNnet_path <- coefficient_path(PLNnet)
  nodeList<- colnames(coefficients(PLNnet$models[[1]]))
  maxPenalty <- PLNnet_path %>% filter(Coeff != 0) %>% group_by(Edge) %>% summarize(max_pen = max(Penalty))
  nodeNb = length(nodeList); edgePenalty = matrix(0, nodeNb, nodeNb)
  for (r in 1:nrow(maxPenalty)) {
    # r = 13
    pairName = as.character(maxPenalty$Edge[r]); pos <- gregexpr(pattern='\\|', pairName)[[1]][[1]]
    node1 <- substr(pairName, 1, (pos-1)); node2 <- substr(pairName, (pos+1), nchar(pairName))
    j <- which(nodeList == node1); k <- which(nodeList == node2)
    edgePenalty[j, k] <- edgePenalty[k, j] <- maxPenalty$max_pen[r]
  }
  edgePenalty[lower.tri(edgePenalty, diag=FALSE)]
}

## Mollusc data
data("mollusk")
mollusc <- prepare_data(mollusk$Abundance, mollusk$Covariate)

## Adjust min_ratio penalty so that a comparable number of deges is inferred for the lower penalty
network_mollusc_raw  <- PLNnetwork(Abundance ~ 1        + offset(log(Offset)), data = mollusc, control = PLNnetwork_param(n_penalties = 50, min_ratio = 0.005))
network_mollusc_site <- PLNnetwork(Abundance ~ 0 + site + offset(log(Offset)), data = mollusc, control = PLNnetwork_param(n_penalties = 50, min_ratio = 0.001))

edgePenalty_raw  <- F_GetEdgePenalty(network_mollusc_raw)
edgePenalty_site <- F_GetEdgePenalty(network_mollusc_site)

comparison <- F_WeightedTau(edgePenalty_raw, edgePenalty_site, weight.method = 'inverse')
comparison[c("tau", "stat", "pval")]




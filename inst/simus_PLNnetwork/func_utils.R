mse <- function(theta, theta.star) {
  mean((theta - theta.star)^2, na.rm = TRUE)
}

ham <- function(theta, theta.star) {
  return( sum( ((as.vector(theta) != 0)*1) != ((as.vector(theta.star)!=0)*1) , na.rm = TRUE) )
}

perf_mse <- function(theta.hat, theta.star) {

  if (is.list(theta.hat)) {
    return(sapply(theta.hat, mse, theta.star))
  } else {
    return(mse(theta.hat, theta.star))
  }
}

perf_ham <- function(theta.hat, theta.star) {

  if (is.list(theta.hat)) {
    return(sapply(theta.hat, ham, theta.star))
  } else {
    return(ham(theta.hat, theta.star))
  }
}

trapeze <- function(x, y) {
  dx <- abs(diff(x))
  area <- sum(c(y[-1]*dx, y[-length(y)]*dx))
  area
}

roc_stabpath <- function(path, true_net) {

  edges_stab <- order(apply(path$prob, 1, function(y) trapeze(path$penalties, y)), decreasing = TRUE)

  edges <- which(true_net[upper.tri(true_net)] != 0)
  P <- length(edges)
  d <- ncol(true_net)
  N <- d*(d - 1)/2 - P

  TP <- cumsum(edges_stab %in% edges)
  FP <- cumsum(!(edges_stab %in% edges))
  recall    <- TP/P         ; recall[P == 0] <- NA
  fallout   <- FP/N         ; fallout[N == 0] <- NA
  precision <- TP/(TP + FP) ; precision[TP + FP == 0] <- NA

  res <- data.frame(
          fallout   = round(fallout  , 3),
          recall    = round(recall   , 3),
          precision = round(precision, 3)
        )
  res
}

perf_roc <- function(theta.hat, theta.star) {

  diag(theta.star) <- 0
  roc <- function(theta) {
    diag(theta) <- 0

    nzero <- which(theta != 0)
    zero  <- which(theta == 0)

    true.nzero <- which(theta.star != 0)
    true.zero  <- which(theta.star == 0)

    TP <- sum(nzero %in% true.nzero)
    TN <- sum(zero %in%  true.zero)
    FP <- sum(nzero %in% true.zero)
    FN <- sum(zero %in%  true.nzero)

    recall    <- TP/(TP + FN) ## also recall and sensitivity
    fallout   <- FP/(FP + TN) ## also 1 - specificit
    precision <- TP/(TP + FP) ## also PPR
    recall[TP + FN == 0] <- NA
    fallout[TN + FP == 0] <- NA
    precision[TP + FP == 0] <- NA

    res <-  round(c(fallout,recall,precision),3)
    res[is.nan(res)] <- 0
    names(res) <- c("fallout","recall", "precision")
    res
  }

  if (is.list(theta.hat)) {
    return(as.data.frame(do.call(rbind, lapply(theta.hat, roc))))
  } else {
    return(roc(theta.hat))
  }
}


perf_auc <- function(roc, threshold = 1) {
  cut <- (roc$fallout < threshold) & (roc$recall < threshold)
  fallout <- c(0, roc$fallout[cut], threshold)
  recall  <- c(0, roc$recall[cut] , threshold)
  dx <- diff(fallout)
  res <- sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2
  res <- ifelse(is.character(res), NA, res)
  res
}

perf_aupr <- function(roc) {
  no_na <- !is.na(roc$precision)
  precision <- roc$precision[no_na]
  recall <- roc$recall[no_na]
  precision <- c(precision,precision[length(precision)])
  recall    <- c(recall, 1)
  dx <- diff(recall)
  res <- sum(c(precision[-1]*dx, precision[-length(precision)]*dx))/2
  res <- ifelse(is.character(res), NA, res)
  res
}

## Ajout par Mahendra
library(igraph)

edge_to_node <- function(x, n = max(x)) {
  x <- x-1 ## easier for arithmetic to number edges starting from 0
  n.node <- round((1 + sqrt(1+8*n)) / 2) ## n.node * (n.node -1) / 2 = n (if integer)
  j.grid <- cumsum(0:n.node)
  j <- findInterval(x, vec = j.grid)
  i <- x - j.grid[j]
  ## Renumber i and j starting from 1 to stick with R convention
  return(data.frame(node1 = i+1, node2 = j+1))
}

node_pair_to_egde <- function(x, y, node.set = union(x, y)) {
  ## Convert node labels to integers (starting from 0)
  x <- match(x, node.set) - 1
  y <- match(y, node.set) - 1
  ## For each pair (x,y) return, corresponding edge number
  n <- length(node.set)
  j.grid <- cumsum(0:(n-1))
  x + j.grid[y] + 1
}

format_stabsel <- function(stabsel, node.names = NULL) {
  res <- stabsel$prob
  colnames(res) <- stabsel$penalties
  res <- res %>%
    as.data.frame() %>%
    mutate(Edge = 1:n()) %>%
    gather(key = "Penalty", value = "Prob", -Edge) %>%
    mutate(Penalty = as.numeric(Penalty),
           Node1   = edge_to_node(Edge)$node1,
           Node2   = edge_to_node(Edge)$node2)
  if (!is.null(node.names)) {
    res <- res %>%
      mutate(Node1 = factor(node.names[Node1], levels = node.names),
             Node2 = factor(node.names[Node2], levels = node.names))
  }
  res
}

plot_stabsel <- function(stabsel, node.names = NULL) {
  p <- format_stabsel(stabsel, node.names) %>%
    ggplot(aes(x = Penalty, y = Prob, group = Edge)) +
    geom_line()
  p
}

entropy <- function(x) {
  x <- x[x>0 & x<1]
  probs <- cbind(x, 1-x)
  -sum(probs * log(probs))
}

stability <- function(x, type = c("entropy", "stars")) {
  type <- match.arg(type)
  n <- length(x)
  metric <- switch(type,
                   entropy = 1 - entropy(x) / (n*log(2)),
                   stars   = 1 - mean(4*x*(1-x)))
  metric
}

select_penalty <- function(stabsel, beta = 0.05) {
  df <- format_stabsel(stabsel) %>%
    group_by(Penalty) %>%
    summarize(Stability = stability(Prob, type = "stars")) %>%
    arrange(-Penalty) %>% ## FRAGILE, penalty in same order as in original network model
    mutate(Penalty.Nb = 1:n()) %>%
    filter(Stability >= 1 - 2*beta)
  c("Penalty" = min(df$Penalty), "Penalty.Nb" = as.integer(max(df$Penalty.Nb)))
}

coefficient_path <- function(networks, precision = TRUE, corr = TRUE, node.labels = NULL) {
  penalties <- networks$penalties
  lapply(penalties, function(x) {
    if (precision) {
      G <- networks$getModel(x)$model_par$Omega
    } else {
      G <- networks$getModel(x)$model_par$Sigma
      dimnames(G) <- dimnames(networks$getModel(x)$model_par$Omega)
    }
    if (corr) {
      G <- ifelse(precision, -1, 1) * G / tcrossprod(sqrt(diag(G)))
    }
    G <- G %>% melt(value.name = "Coeff", varnames = c("Node1", "Node2"))
    if (!is.null(node.labels)) {
     full.names <- unique(G$Node1)
     G <- G %>%
       mutate(Node1 = factor(Node1, levels = full.names, labels = node.labels),
              Node2 = factor(Node2, levels = full.names, labels = node.labels))
    }
    G %>%
      mutate(Penalty = x,
             Node1   = as.character(Node1),
             Node2   = as.character(Node2)) %>%
      filter(Node1 < Node2) %>%
      mutate(Edge = node_pair_to_egde(Node1, Node2))
  }) %>% bind_rows()
}

density_path <- function(networks) {
  penalties <- networks$penalties
  data.frame(Penalty = penalties,
             Density = sapply(penalties,
                              function(x) { networks$getModel(x)$latent_network() %>% mean() }),
             stringsAsFactors = FALSE)
}

## alternative plot method for graph of partial correlation
plot_partial_corr <- function(x, edge.color = c("#F8766D", "#00BFC4"),
                              remove.isolated = FALSE, node.labels = NULL) {
  net <- x$model_par$Omega
  net <- - net / tcrossprod(sqrt(diag(net)))
  G <- graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)
  if (!is.null(node.labels)) {
    V(G)$label <- node.labels
  } else {
    V(G)$label <- colnames(net)
  }
  ## Nice nodes
  V.deg <- degree(G)/sum(degree(G))
  V(G)$label.cex <- V.deg / max(V.deg) + .5
  V(G)$size <- V.deg * 100
  V(G)$label.color <- rgb(0, 0, .2, .8)
  V(G)$frame.color <- NA
  ## Nice edges
  E(G)$color <-ifelse(E(G)$weight > 0, edge.color[1], edge.color[2])
  E(G)$width <- 15*abs(E(G)$weight)

  if (remove.isolated) {
    G <- delete.vertices(G, which(degree(G) == 0))
  }

  ## Plot
  plot(G, layout = layout_in_circle)
  invisible(G)
}

auc <- function(fallout, recall, threshold = 1) {

  if (tail(fallout, 1) - head(fallout, 1) < 0) {
    fallout <- rev(fallout)
    recall  <- rev(recall)
  }

  cut <- (fallout < threshold) & (recall < threshold)
  fallout <- c(0, fallout[cut], threshold)
  recall  <- c(0, recall[cut] , threshold)
  dx <- diff(fallout)
  res <- sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2
  res <- ifelse(is.character(res), NA, res)
  res
}

aupr <- function(recall, precision) {

  if (tail(recall, 1) - head(recall, 1) < 0) {
    precision <- rev(precision)
    recall  <- rev(recall)
  }

  no_na <- !is.na(precision)
  precision <- precision[no_na]
  recall <- recall[no_na]
  precision <- c(precision,precision[length(precision)])
  recall    <- c(recall, 1)
  dx <- diff(recall)
  res <- sum(c(precision[-1]*dx, precision[-length(precision)]*dx))/2
  res <- ifelse(is.character(res), NA, res)
  res
}



## This R script is dedicated to graph analysis
## It encompasses: 
## Nodal Measures: strength, betweenness centrality, clustering coefficient, efficiency
## Global Measures: strength, modularity, efficiency
## Community Structure: spectural clustering

library(igraph)
## Nodal Measures -------------------------------------------------
nodal_strength <- function(g) {
  strength_vec <- strength(g, mode = "all", weights = E(g)$weight)
  print(head(strength_vec))
  
  strength_vec ## Returns a numeric vector
}

nodal_betweenness <- function(g) {
  distance <- 1/E(g)$weight ## Distance between nodes as inverse edge weight
  betw_centrality <- betweenness(g, directed = F,
                                 weights = inv_weights, 
                                 normalized = T)
  print(head(betw_centrality)) 
  
  betw_centrality ## Returns a numeric vector
}

nodal_efficiency <- function(g) {
  sapply(V(g), function(v) {
   neighbors_v <- neighbors(g, v) 
   subgraph <- induced_subgraph(g, neighbors_v)
   
   if (vcount(subgraph) < 2) {
     return(NA) ## efficiency undefined for < 2 nodes
   }
   
   ## The distance between any node pair, distance is inverse weight
   sp_lengths <- distances(subgraph, weights = 1/E(subgraph)$weight)
   inv_sp <- 1/abs(sp_lengths) 
   diag(inv_sp) <- 0 #ignore self_distances
   
   n <- vcount(subgraph)
   eff <- sum(inv_sp) / n*(n-1)
   return(eff)
  })
}

nodal_clust_coef <- function(g) {
  transitivity(g, type = "local", weight = E(g)$weight) ## Weighted coef based on A. Barrat
}

## Global Measures ------------------------------------------------




### Spectural Clustering Code ----------------------------------------
## Compute the Similarity Matrix with Gaussian Radial Kernel
S <- function(X) {
  exp(-0.5 * as.matrix(dist(X)))
}

## Weighted adjacency matrix with 10 nearest neighbors
knn = function(values, K){
  index = order(values, decreasing = TRUE)[(1:K) + 1]
  values[-index] = NA
  return(values)
}


## Compute the weight Matrix
W <- function(S, type = "Full", extra = NA) {
  if (type == "Soft") {
    lambda <- extra$lambda 
    W <- S - lambda
    W[W < 0 ] <- 0
    
    ## Epsilon neighborhood
  } else if (type == "Epsilon") {
    S > epsilon
    
    ## Fully Connected Graph
  } else {
    W <- apply(S, 2, function(x) knn(x, K = 10))
    W = (W + t(W)) / 2
    W[is.na(W)] = 0
    W
  }
}

## Compute the Laplacian
L <- function(W) {
  D <- diag(rowSums(W))
  D - W
}

## Pick the K smallest eigenvalues of L and corrsponding eigenvectors
bottomK <- function(L, K) {
  n <- nrow(L)
  E <- eigen(L)
  U <- E$vectors[,(n-K+1):n]
  
  U
}

## K-mean Cluster rows of U
KM <- function(U, K) {
  kmeans(U, centers = K, iter.max = 100, nstart = 10)$cluster
}

## Wrapper Function
spectralClust <- function(X, K, cancerType = "All", extra = NA) {
  s <- S(X) ## Similarity
  w <- W(s, "Full", extra) ## Weight
  l <- L(w) ## Laplacian
  u <- bottomK(l, K) ## Eigen values
  clustVec <- KM(u, K)
  
  par(mfrow = c(1,2))
  ## Plot first two eigenvectors, colored by cluster
  library(ggplot2)
  U2 <- data.frame(cbind(u[,1:2], clustVec))
  colnames(U2) <- c("U1",  "U2", "C")
  
  plot(U2[,1], U2[,2], col = factor(clustVec), 
       main = "Clustering in top 2 eigenvectors",
       xlab = "U1", ylab = "U2")
  
  ## Plot Grpah Rrepresentation
  G = graph_from_adjacency_matrix(w, mode="undirected", weighted = TRUE)
  layout_circle <- layout_in_circle(G)
  plot(G, vertex.size=3, vertex.label=NA, 
       layout = layout_circle, 
       main = paste("Spectral Graph for ", cancerType))
  
  mtext(cancerType, 
        side = 3, line = -2, outer = TRUE, cex = 1.5)
  
  list("clustVec" = clustVec, 
       "G" = G)
}

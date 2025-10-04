## This script is dedicated for partial correlation calculation through GLasso

## Create the covariance matrix with input X matrix
createCov <- function(X) {cov(X)}

## glasso for finding the penalized inverse covariance matrix
## convert inverse covariance to partial correlation matrix
## @S: an n*n matrix, the covariance matrix
## @rho: a float, the penalization coefficient
## @n: an integer, the sample size
parCorr_from_precision <- function(precision, rho, n) {
  ## Compute partial correlation matrix from inverse covariance matrix
  invCovarMat <- precision
  invCovarDiag <- diag(invCovarMat)
  partCovarMat <- -invCovarMat / sqrt(outer(invCovarDiag, invCovarDiag))
  
  symmPartCovarMat <- (partCovarMat + t(partCovarMat))/2
  
  sparseCovInvGraph <-
    graph_from_adjacency_matrix(symmPartCovarMat,
                                mode = "undirected",
                                weighted = T) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    filter(from != to)
}

## Finds the optimal rho value based on extended BIC for graphs
## @S: an n*n matrix, the covariance matrix
## @rhoList: a numeric vector, a list of penalization coefficients
## @n: an, integer, the sample size
find_optimal_rho <- function(S, rhoList, n) {
  
}


## Plot the partial correlation graph in a circle for one sample
## g: igraph class, the input graph object
## rho: a float, the penalization coefficient
## status: a character, the disease status, healthy or any cancer type
plot_parcorr_graph <- function(g, rho, status) {
  ggraph(g, 
         layout = 'linear', circular = TRUE) + 
    geom_node_point()+
    geom_edge_arc(aes(colour = weight)) + 
    scale_edge_colour_gradient(low="grey",
                               high="black",
                               limits = c(0, 0.4))+
    coord_fixed()+
    labs(subtitle = paste("rho = 1,", status),
         colour = "Pcorr")
}

## Plot a list of ggplots in 3*3 grid, save the plot into the image directory
## @graphPlotList: a list of ggplot object
## @criteria: a character string, the criteria by which plots are created separately.
plotgrid_graph <- function(graphPlotList, criteria){
  plot_grid(plotlist = graphPlotList,
            nrow = 3, ncol = 3)
  ggsave(paste0("../Images/PartialCorrGraphs_by", criteria, ".png"))
}


pcaCombo <- function(X) {
  PCA <- prcomp(X, center = T, scale = F, rank. = min(dim(X)))
  
  ## Compute Scree Plot
  totalVar <- sum(PCA$sdev^2)
  varExpl <- PCA$sdev^2/totalVar
  
  # PCs that explain > 1% of variances
  minPC <- which(varExpl > 0.01)
  percPCmin <- varExpl[minPC]
  scree <- data.frame(
    "k" = 1:min(dim(X)),
    "perc" = varExpl
  )
  
  limitPC <- which.min(-diff(percPCmin, 1) < 0.03)
  
  
  A <- ggplot(scree) + 
    geom_line(aes(x = k, y = perc, 
                  color = adjustcolor("red", alpha = 0.4)), 
              stat = "identity")+
    geom_hline(linetype = "dashed", yintercept = percPCmin[limitPC]) + 
    geom_vline(linetype = "dashed", xintercept = limitPC) +
    xlab("# of PCs") + 
    ylab("Proportion of variance explained")+
    ggtitle("Scree Plot") + 
    theme(legend.position = "none") 
  
  ## Compute PC1 vs PC 2 Scatter Plot
  Z3 <- scale(X) %*% PCA$rotation[,1:3]
  Z3 <- data.frame(Z3)
  colnames(Z3) <- c("PC1", "PC2", "PC3")
  
  Z3$cancerType <- cancerTypes
  B <- ggplot(Z3) + 
    geom_point(aes(x = PC1, y = PC2, color = cancerTypes, alpha = 0.3))+
    ggtitle("PC 1,2 Scatter Plot")+
    scale_color_brewer(palette = "Paired")
  
  ## PC 1 and PC 2 Trace Plot
  V3 <- data.frame(PCA$rotation[,1:3])
  colnames(V3) <- c("PC1", "PC2")
  V3$k <- 1:dim(V3)[1]
  
  C <- ggplot(V3) + 
    geom_bar(aes(x = k, y = PC1), color = adjustcolor("red", alpha = 0.4),  stat = "identity")+
    geom_hline(linetype = "dashed", yintercept = 0.05) +
    geom_hline(linetype = "dashed", yintercept = -0.05) +
    xlab("Genes") + 
    ylab("PC 1") + 
    theme(legend.position = "none")+
    ggtitle("PC 1 Trace Plot")
  
  D <- ggplot(V3) + 
    geom_bar(aes(x = k, y = PC2), color = adjustcolor("blue", alpha = 0.4), stat = "identity")+
    geom_hline(linetype = "dashed", yintercept = 0.05) +
    geom_hline(linetype = "dashed", yintercept = -0.05) +
    xlab("Genes") + 
    ylab("PC 2") + 
    theme(legend.position = "none") + 
    ggtitle("PC 2 Trace Plot")
  
  cowplot::plot_grid(plotlist = list(A,B,C,D),
                     nrow = 2, ncol = 2,
                     labels = c("A", "B", "C", "D")) %>%
    print()
  
  return(V3)
  
}

sparsePcaCombo <- function(X, budget, Kmax, cancerTypes) {
  library(PMA)
  PCA <- SPC(X, sumabsv = budget, K = Kmax, trace = FALSE, 
             compute.pve = TRUE, vpos = FALSE, vneg = FALSE)
  
  ## Compute Scree Plot
  varExpl <- PCA$prop.var.explained
  
  # # PCs that explain > 1% of variances
  # minPC <- which.max(varExpl > 0.01)
  # percPCmin <- varExpl[minPC]
  scree <- data.frame(
    "k" = 1:Kmax,
    "perc" = varExpl
  )
  
  # limitPC <- which.min(-diff(percPCmin, 1) < 0.03)
  
  A <- ggplot(scree) + 
    geom_line(aes(x = k, y = perc, 
                  color = adjustcolor("red", alpha = 0.4)), 
              stat = "identity")+
    # geom_hline(linetype = "dashed", yintercept = percPCmin) + 
    # geom_vline(linetype = "dashed", xintercept = minPC) +
    xlab("# of PCs") + 
    ylab("Proportion of variance explained")+
    ggtitle("Scree Plot") + 
    theme(legend.position = "none") 
  
  ## Compute PC1 vs PC 2 Scatter Plot
  Vall <- PCA$v
  V3 <- Vall[,1:3]
  Z3 <- scale(X) %*% V3
  Z3 <- data.frame(Z3)
  
  colnames(V3) <- c("PC1", "PC2", "PC3")
  colnames(Z3) <- c("PC1", "PC2", "PC3")
  
  Z3$Types <- cancerTypes
  
  B <- ggplot(Z3) + 
    geom_point(aes(x = PC1, y = PC2, color = Types, alpha = 0.3))+
    ggtitle("PC 1,2 Scatter Plot")+
    scale_color_brewer(palette = "Paired")
  
  ## PC 1 and PC 2 Trace Plot
  V3 <- data.frame(V3)
  V3$k <- 1:dim(V3)[1]
  
  C <- ggplot(V3) + 
    geom_bar(aes(x = k, y = PC1), color = adjustcolor("red", alpha = 0.4),  stat = "identity")+
    geom_hline(linetype = "dashed", yintercept = 0.05) +
    geom_hline(linetype = "dashed", yintercept = -0.05) +
    xlab("Genes") + 
    ylab("PC 1") + 
    theme(legend.position = "none")+
    ggtitle("PC 1 bar plot")
  
  D <- ggplot(V3) + 
    geom_bar(aes(x = k, y = PC2), color = adjustcolor("blue", alpha = 0.4), stat = "identity")+
    geom_hline(linetype = "dashed", yintercept = 0.05) +
    geom_hline(linetype = "dashed", yintercept = -0.05) +
    xlab("Genes") + 
    ylab("PC 2") + 
    theme(legend.position = "none") + 
    ggtitle("PC 2 bar plot")
  
  cowplot::plot_grid(plotlist = list(A,B,C,D),
                     nrow = 2, ncol = 2,
                     labels = c("A", "B", "C", "D")) %>%
    print()
  
  score3D(Z3, paste('Sparse PCA, K = ', Kmax, "L1 norm =", budget), "cancerTypes")
  
  return(Vall)
}
makeScore3D <- function(df, PCM, cancerTypes) {
  ## df: n*p matrix, rows are observations, columns are features
  ## PCM: PC matrix with the top 3 PCs
  score <- as.matrix(df) %*% as.matrix(PCM)
  colnames(score) <- paste0("PC", 1:ncol(score))
  score <- cbind(score, 
                 "Healthy" = ifelse(cancerTypes == "Normal", 
                                    "Healthy", "Cancer"))
  score <- cbind(score, "Types"= cancerTypes)
  score <- data.frame(score)
  
  score
}


score3D <- function(df, title = 'Score Scatter Plot of ', 
                    colorFactor = NA) {
  library(plotly)
  if (is.na(colorFactor)) {
    fig <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, 
                   size = 10) 
  } else {
    fig <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, 
                   color = ~eval(sym(colorFactor)), 
                   colors = "Paired",
                   size = 10) 
  }
  
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC 1'),
                                     yaxis = list(title = 'PC 2'),
                                     zaxis = list(title = 'PC 3')),
                        title = title,
                        legend = list(title=list(text='<b> Categories </b>'))
                        )
  
  fig
}

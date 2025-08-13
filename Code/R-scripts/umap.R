UMAP.2d <- function(X, cancerTypes) {
  UMAP = umap(X, config = umap.defaults)
  
  colnames(UMAP$layout) = c("x1", "x2")
  df = data.frame(UMAP$layout)
  df$Types = as.factor(cancerTypes)
  
  write.csv(df, file = "../Data/umap-output-2d.csv")
  
  ggplot(data = df)+
    geom_point(aes(x = x1, y = x2, color = Types)) +
    labs(title = "UMAP Scatter Plot on Top 200 Genes") +
    xlab("UMAP 1") + xlab("UMAP 2") +
    scale_color_brewer(palette = "Paired") %>%
    show()
  
}

UMAP.3d <- function(X, cancerTypes) {
  UMAP = umap(X, config = umap.defaults)
  browser()
  colnames(UMAP$layout) = c("x1", "x2", "x3")
  df = data.frame(UMAP$layout)
  df$Types = as.factor(cancerTypes)
  
  write.csv(df, file = "../Data/umap-output-3d.csv")
  
  fig <- plot_ly(df, x = ~x1, y = ~x2, z = ~x3, 
                 color = ~Types, 
                 colors = "Paired",
                 size = 10)
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'UMAP 1'),
                                     yaxis = list(title = 'UMAP 2'),
                                     zaxis = list(title = 'UMAP 3')),
                        title = title,
                        legend = list(title=list(text='<b> Categories </b>'))
  )
  
  fig
}

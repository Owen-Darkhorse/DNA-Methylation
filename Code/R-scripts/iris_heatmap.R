heatmap = function(A){
  n = nrow(A)
  p = ncol(A)  
  df = data.frame(value = c(A),  i = 1:n, j = rep(1:p, rep(n, p)))
  ggplot(df, aes(j, i, fill = value)) + 
    geom_tile()+
    scale_fill_viridis_c()+
    # scale_fill_gradient(low="darkblue", high="yellow", mid = "darkgreen")+
    scale_y_reverse()+
    theme_void()
}

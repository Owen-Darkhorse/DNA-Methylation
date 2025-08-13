IVT <- function(X, k) {
  apply(X, MARGIN = 2, function(x) IVT.column(x, k), simplify = T)
}

IVT.column <- function(x, k) {
  n <- length(x)
  x.rank <- rank(x)
  r <- (x.rank-k)/(n-2*k + 1)
  qnorm(r)
}
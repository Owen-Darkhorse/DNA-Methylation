INT <- function(X, k) {
  apply(X, MARGIN = 2, function(x) INT.column(x, k), simplify = T)
}

## Monotonically transform the raw data into the normal distribution's scale
INT.column <- function(x, k) {
  n <- length(x)
  nonZeroCount <- sum(x>0) + 1
  ranks <- rep(1, length(x)) ## Initialize all ranks as 1
  ranks[x > 0] <- rank(x[x>0]) + 1 # For positive methylation, assign positive ranks to them
  
  r <- (ranks-k)/(nonZeroCount-2*k + 1)
  qnorm(r)
}

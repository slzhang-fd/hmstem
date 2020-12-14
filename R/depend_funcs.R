cov2corr <- function(x){
  K <- sqrt(length(c(x)))
  x <- matrix(c(x), K, K)
  scale <- diag(1 / sqrt(diag(x)), nrow = K, ncol = K)
  return(c(scale %*% x %*% scale))
}
summary_res_basic <- function(x, name, rownum, colnum, dimnames){
  res <- list(matrix(colMeans(matrix(x, ncol = rownum * colnum)), rownum, colnum, dimnames=dimnames),
              matrix(apply(matrix(x, ncol = rownum * colnum), 2, sd), rownum, colnum, dimnames=dimnames),
              matrix(apply(matrix(x, ncol = rownum * colnum), 2, quantile, 0.025), rownum, colnum, dimnames=dimnames),
              matrix(apply(matrix(x, ncol = rownum * colnum), 2, quantile, 0.975), rownum, colnum, dimnames=dimnames))
  names(res) <- paste0(name, c('', '_sd', '_low2.5', '_high97.5'))
  return(res)
}

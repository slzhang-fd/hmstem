sigma_hessian_expect <- function(sigma_inv, residu, calcu_diag){
  K <- nrow(sigma_inv)
  off_num <- K*(K-1)/2 + calcu_diag * K
  residu_num <- nrow(residu)
  h_res <- matrix(0, off_num, off_num)
  for(k1 in 1:off_num){
    I_tmp1 <- matrix(0, K, K)
    I_tmp1[lower.tri(I_tmp1, diag = calcu_diag)][k1] <- 1
    I_tmp1[upper.tri(I_tmp1)] <- t(I_tmp1)[upper.tri(I_tmp1)]

    for(k2 in k1:off_num){
      I_tmp2 <- matrix(0, K, K)
      I_tmp2[lower.tri(I_tmp2, diag = calcu_diag)][k2] <- 1
      I_tmp2[upper.tri(I_tmp2)] <- t(I_tmp2)[upper.tri(I_tmp2)]

      h_res[k1,k2] <- h_res[k2,k1] <-
        - sum(diag(residu_num * sigma_inv %*% I_tmp1 %*% sigma_inv %*% I_tmp2)) / 2
    }
  }
  h_res
}

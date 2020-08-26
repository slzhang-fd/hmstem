sample_Y_star_3fac <- function(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e, mu_v_u_w){
  K <- ncol(Y_star)
  # mu_v_u_w <- Xbeta + V_all[j_ind,] + U_all[i_ind,] + W_all[jt_ind,]
  for(k in 1:K){
    sigma11 <- Sigma_e[k,k]
    sigma12 <- Sigma_e[k,-k]
    sigma22_inv <- solve(Sigma_e[-k,-k])
    sigma_Y_cond <- sigma11 - sigma12 %*% sigma22_inv %*% sigma12
    mu_Y_cond <- mu_v_u_w[,k] + (Y_star[,-k] - mu_v_u_w[,-k]) %*% sigma22_inv %*% sigma12

    Y_star[Y_pos_ind[,k],k] <- rtruncnorm(1, a=0, mean = mu_Y_cond[Y_pos_ind[,k]], sd = sqrt(sigma_Y_cond))
    Y_star[Y_zero_ind[,k],k] <- rtruncnorm(1, b=0, mean = mu_Y_cond[Y_zero_ind[,k]], sd = sqrt(sigma_Y_cond))
  }
  Y_star
}

sample_X_all1_3fac <- function(X_all, Sigma_e_inv, Sigma_x_inv, x_len, temp){
  # temp <- rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T)
  for(pp in unique(x_len)){
    sigma_X_cond <- solve(pp * Sigma_e_inv + Sigma_x_inv)
    x_loc <- which(x_len==pp)
    mu_X_cond <- temp[x_loc,] %*% Sigma_e_inv %*% sigma_X_cond
    X_all[x_loc,] <- mu_X_cond + rmvnorm(length(x_loc), sigma = sigma_X_cond)
  }
  X_all
}
Optim_step_3fac <- function(x_covs, Y_star, U_all, V_all, W_all, coeffs,
                            B_e, Sigma_e_inv, B_u, B_v, B_w, j_ind, i_ind, jt_ind){
  K <- ncol(Y_star)
  p <- nrow(coeffs)
  ## update coeffs
  params <- optim(par = c(coeffs), fn = neg_loglik_Y_star_coeff, gr = neg_loglik_Y_star_coeff_deriv,
                  x_covs=x_covs, temp0=Y_star - V_all[j_ind,] - U_all[i_ind,] - W_all[jt_ind,],
                  Sigma_e_inv=Sigma_e_inv,
                  method = "BFGS")$par
  coeffs <- matrix(params,p,K)
  ## update Sigma_e
  params <- optim(par = c(B_e[lower.tri(B_e)]), fn = neg_loglik_Y_star_sigma,
                  temp=Y_star - x_covs %*% coeffs - V_all[j_ind,] - U_all[i_ind,] - W_all[jt_ind,],
                  method = "BFGS")$par
  B_e[lower.tri(B_e)] <- params
  # B_e <- diag(1/sqrt(rowSums(B_e^2))) %*% B_e

  ## update Sigma_u
  params <- optim(par = B_u[lower.tri(B_u,diag = T)], fn = neg_loglik_mvnorm,
                  x = U_all, method = "BFGS")$par
  B_u[lower.tri(B_u, diag = T)] <- params

  ## update Sigma_v
  params <- optim(par = B_v[lower.tri(B_v,diag = T)], fn = neg_loglik_mvnorm,
                  x = V_all, method = "BFGS")$par
  B_v[lower.tri(B_v, diag = T)] <- params

  ## update Sigma_w
  params <- optim(par = B_w[lower.tri(B_w,diag = T)], fn = neg_loglik_mvnorm,
                  x = W_all, method = "BFGS")$par
  B_w[lower.tri(B_w, diag = T)] <- params

  list('coeffs' = coeffs,
       'B_e' = B_e,
       'B_u' = B_u,
       'B_v' = B_v,
       'B_w' = B_w)
}



#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
sample_Y_star <- function(Y_pos_ind, Y_zero_ind, Y_star, U_all, V_all, Xbeta, Sigma_e, j_ind, i_ind){
  K <- ncol(Y_star)
  NN <- nrow(Y_star)
  mu_Y_cond <- rep(0, NN)
  sigma_Y_cond <- 1
  mu_v_u <- Xbeta + V_all[j_ind,] + U_all[i_ind,]

  for(k in 1:K){
    sigma11 <- Sigma_e[k,k]
    sigma12 <- Sigma_e[k,-k]
    sigma22_inv <- solve(Sigma_e[-k,-k])
    sigma_Y_cond <- sigma11 - sigma12 %*% sigma22_inv %*% sigma12
    mu_Y_cond <- mu_v_u[,k] + (Y_star[,-k] - mu_v_u[,-k]) %*% sigma22_inv %*% sigma12

    Y_star[Y_pos_ind[,k],k] <- rtruncnorm(1, a=0, mean = mu_Y_cond[Y_pos_ind[,k]], sd = sqrt(sigma_Y_cond))
    Y_star[Y_zero_ind[,k],k] <- rtruncnorm(1, b=0, mean = mu_Y_cond[Y_zero_ind[,k]], sd = sqrt(sigma_Y_cond))
  }
  Y_star
}

sample_U_all <- function(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_u_inv, t_cums, j_ind){
  sigma_U_cond <- list()
  for(pp in unique(diff(t_cums))){
    sigma_U_cond[[pp]] <- solve(pp * Sigma_e_inv + Sigma_u_inv)
  }
  temp <- Y_star - Xbeta - V_all[j_ind,]
  for(i in 1:nrow(U_all)){
    T_i <- t_cums[i+1] - t_cums[i]
    n_i <- t_cums[i] + 1
    mu_U_cond <- sigma_U_cond[[T_i]] %*% Sigma_e_inv %*%
      colSums(matrix(temp[n_i:(n_i+T_i-1),],nrow=T_i))
    U_all[i,] <- rmvnorm(1, mean = mu_U_cond, sigma = sigma_U_cond[[T_i]]) ## U_all[i_ind[n_i],] <-
  }
  U_all
}
sample_U_all1 <- function(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_u_inv, t_len, j_ind, i_ind){
  temp <- rowsum(Y_star - Xbeta - V_all[j_ind,], i_ind, reorder = T)
  for(pp in unique(t_len)){
    sigma_U_cond <- solve(pp * Sigma_e_inv + Sigma_u_inv)
    i_loc <- which(t_len==pp)
    mu_U_cond <- temp[i_loc,] %*% Sigma_e_inv %*% sigma_U_cond
    U_all[i_loc,] <- mu_U_cond + rmvnorm(length(i_loc), sigma = sigma_U_cond)
  }
  U_all
}

sample_V_all <- function(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_v_inv, it_cums, i_ind){
  sigma_V_cond <- list()
  for(pp in unique(diff(it_cums))){
    sigma_V_cond[[pp]] <- solve(pp * Sigma_e_inv + Sigma_v_inv)
  }
  temp <- Y_star - Xbeta - U_all[i_ind,]
  for(j in 1:nrow(V_all)){
    T_j <- it_cums[j+1] - it_cums[j]
    n_j <- it_cums[j] + 1
    mu_V_cond <- sigma_V_cond[[T_j]] %*% Sigma_e_inv %*%
      colSums(matrix(temp[n_j:(n_j+T_j-1),],nrow = T_j))
    V_all[j,] <- rmvnorm(1, mean = mu_V_cond, sigma_V_cond[[T_j]])
  }
  V_all
}
sample_V_all1 <- function(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_v_inv, it_len, j_ind, i_ind){
  temp <- rowsum(Y_star - Xbeta - U_all[i_ind,], j_ind, reorder = T)
  for(pp in unique(it_len)){
    sigma_V_cond <- solve(pp * Sigma_e_inv + Sigma_v_inv)
    j_loc <- which(it_len==pp)
    mu_V_cond <- temp[j_loc,] %*% Sigma_e_inv %*% sigma_V_cond
    V_all[j_loc,] <- mu_V_cond + rmvnorm(length(j_loc), sigma = sigma_V_cond)
  }
  V_all
}
## Sigma = B %*% t(B) or B = t(chol(Sigma))
## return -2 * loglik / NN
neg_loglik_Y_star_coeff <- function(params, x_covs, temp0, Sigma_e_inv){
  NN <- nrow(temp0)
  K <- ncol(temp0)
  p <- ncol(x_covs)

  coeffs <- matrix(params,p,K)

  temp <- temp0 - x_covs %*% coeffs
  sum(diag(Sigma_e_inv %*% t(temp) %*% temp)) / NN
}
neg_loglik_Y_star_coeff_deriv <- function(params, x_covs, temp0, Sigma_e_inv){
  NN <- nrow(temp0)
  K <- ncol(temp0)
  p <- ncol(x_covs)

  coeffs <- matrix(params,p,K)
  c(-2 * t(x_covs) %*% (temp0 - x_covs %*% coeffs) %*% Sigma_e_inv / NN)
}
neg_loglik_Y_star_sigma <- function(params, temp){
  NN <- nrow(temp)
  K <- ncol(temp)
  B <- diag(rep(1,K))
  B[lower.tri(B)] <- params
  Sigma_e <- B %*% t(B)
  Sigma_e_inv <- chol2inv(t(B))

  # temp <- Y_star - (x_covs %*% coeffs + V_all[j_ind,] + U_all[i_ind,])
  log(det(Sigma_e)) + sum(diag(Sigma_e_inv %*% t(temp) %*% temp)) / NN
}
## multivariate normal density with zero mean and B
## return -2 * loglik / N
neg_loglik_mvnorm <- function(params, x){
  K <- ncol(x)
  N <- nrow(x)
  B <- diag(rep(1,K))
  B[lower.tri(B,diag = T)] <- params
  Sigma <- B %*% t(B)
  # Sigma_inv <- solve(Sigma)
  Sigma_inv <- chol2inv(t(B))
  log(det(Sigma)) + sum(diag(Sigma_inv %*% t(x) %*% x)) / N
}

Optim_step <- function(x_covs, Y_star, U_all, V_all, coeffs, B_e, Sigma_e_inv, B_u, B_v, j_ind, i_ind){
  K <- ncol(Y_star)
  p <- nrow(coeffs)
  ## update coeffs and Sigma_e
  params <- optim(par = c(coeffs), fn = neg_loglik_Y_star_coeff, gr = neg_loglik_Y_star_coeff_deriv,
                  x_covs=x_covs, temp0=Y_star - V_all[j_ind,] - U_all[i_ind,],
                  Sigma_e_inv=Sigma_e_inv,
                  method = "BFGS")$par
  coeffs <- matrix(params,p,K)
  ## update Sigma_e
  params <- optim(par = c(B_e[lower.tri(B_e)]), fn = neg_loglik_Y_star_sigma,
                  temp=Y_star - (x_covs %*% coeffs + V_all[j_ind,] + U_all[i_ind,]),
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

  list('coeffs' = coeffs,
       'B_e' = B_e,
       'B_u' = B_u,
       'B_v' = B_v)
}

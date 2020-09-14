#' @importFrom Matrix bdiag
Optim_step_1fac <- function(x_covs, Y_star, U_all, coeffs,
                            B_e, Sigma_e_inv, B_u, i_ind){
  K <- ncol(Y_star)
  p <- nrow(coeffs)
  ## update coeffs
  params <- optim(par = c(coeffs), fn = neg_loglik_Y_star_coeff, gr = neg_loglik_Y_star_coeff_deriv,
                  x_covs=x_covs, temp0=Y_star - U_all[i_ind,],
                  Sigma_e_inv=Sigma_e_inv,
                  method = "BFGS")$par
  coeffs <- matrix(params,p,K)
  ## update Sigma_e
  if(K>1){
    params <- optim(par = c(B_e[lower.tri(B_e)]), fn = neg_loglik_Y_star_sigma,
                    temp=Y_star - x_covs %*% coeffs - U_all[i_ind,],
                    method = "BFGS")$par
    B_e[lower.tri(B_e)] <- params
  }
  # B_e <- diag(1/sqrt(rowSums(B_e^2))) %*% B_e

  ## update Sigma_u
  params <- optim(par = B_u[lower.tri(B_u,diag = T)], fn = neg_loglik_mvnorm,
                  x = U_all, method = "BFGS")$par
  B_u[lower.tri(B_u, diag = T)] <- params

  list('coeffs' = coeffs,
       'B_e' = B_e,
       'B_u' = B_u)
}
#' @export hiera_1factor
hiera_1factor <- function(Y, x_covs, i_ind, max_steps){
  K <- ncol(Y)
  rcd_num <- nrow(Y)
  ind_num <- max(i_ind)
  cov_num <- ncol(x_covs)

  t_len <- aggregate(i_ind, by=list(i_ind), length)$x
  ## initialize random variables
  Y_star <- matrix(0, rcd_num, K)
  U_all <- matrix(rnorm(ind_num*K),ind_num,K)
  ## initialize parameters
  coeffs <- matrix(0, cov_num, K)
  B_e <- B_u <- diag(rep(1,K))

  coeffs_all <- matrix(0, max_steps, cov_num*K)
  B_e_all <- B_u_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  for(iter in 1:max_steps){
    cat('\r iter: ', iter, 'mu=',coeffs[1,])
    Sigma_e <- B_e %*% t(B_e)
    # Sigma_u <- B_u %*% t(B_u)
    Sigma_e_inv <- chol2inv(t(B_e))
    Sigma_u_inv <- chol2inv(t(B_u))
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e,
                                 Xbeta + U_all[i_ind,])
    cat(' sample Y done |')
    # sample U
    U_all <- sample_X_all1_3fac(U_all, Sigma_e_inv, Sigma_u_inv, t_len,
                                rowsum(Y_star - Xbeta, i_ind, reorder = T))
    cat('sample U done |')
    # update parameters
    params <- Optim_step_1fac(x_covs, Y_star, U_all, coeffs,
                              B_e, Sigma_e_inv, B_u, i_ind)
    coeffs <- params$coeffs
    B_e <- params$B_e
    B_u <- params$B_u
    # store results
    scale_mat <- diag(1/sqrt(rowSums(B_e^2)))
    coeffs_all[iter,] <- coeffs %*% scale_mat
    B_e_all[iter,] <- scale_mat %*% B_e
    B_u_all[iter,] <- scale_mat %*% B_u
  }
  return(list('coeffs_all'=coeffs_all,
              'B_e_all'=B_e_all,
              'B_u_all'=B_u_all,
              'Y_star' = Y_star,
              'U_all' = U_all))
}

#' @export log_lik_1fac
log_lik_1fac <- function(params, x_covs, y_u, u_all) {
  cov_num <- ncol(x_covs)
  k_dims <- ncol(y_u)
  y_nums <- nrow(y_u)
  u_nums <- nrow(u_all)

  # decompose params
  coeffs <- matrix(params[1:(cov_num * k_dims)], cov_num, k_dims)
  sigma_e <- sigma_u <- diag(rep(1, k_dims))
  lower_num <- k_dims * (k_dims - 1) / 2
  sigma_e[lower.tri(sigma_e, diag = F)] <-
    params[cov_num * k_dims + 1:lower_num]
  sigma_u[lower.tri(sigma_u, diag = T)] <-
    params[cov_num * k_dims + lower_num + 1:(lower_num + k_dims)]
  sigma_e[upper.tri(sigma_e)] <- t(sigma_e)[upper.tri(sigma_e)]
  sigma_u[upper.tri(sigma_u)] <- t(sigma_u)[upper.tri(sigma_u)]
  temp <- y_u - x_covs %*% coeffs
  sigma_e_inv <- solve(sigma_e)
  sigma_u_inv <- solve(sigma_u)
  res <- 0.5 * log(det(sigma_e_inv)) * y_nums -
    0.5 * sum(diag(sigma_e_inv %*% t(temp) %*% temp)) +
    0.5 * log(det(sigma_u_inv)) * u_nums -
    0.5 * sum(diag(sigma_u_inv %*% t(u_all) %*% u_all))
  res
}
#' @export log_lik_1fac_deriv
log_lik_1fac_deriv <- function(params, x_covs, y_u, u_all) {
  cov_num <- ncol(x_covs)
  k_dims <- ncol(y_u)
  y_nums <- nrow(y_u)
  u_nums <- nrow(u_all)
  # decompose params
  coeffs <- matrix(params[1:(cov_num * k_dims)], cov_num, k_dims)
  sigma_e <- sigma_u <- diag(rep(1, k_dims))
  lower_num <- k_dims * (k_dims - 1) / 2
  sigma_e[lower.tri(sigma_e, diag = F)] <-
    params[cov_num * k_dims + 1:lower_num]
  sigma_u[lower.tri(sigma_u, diag = T)] <-
    params[cov_num * k_dims + lower_num + 1:(lower_num+k_dims)]

  sigma_e[upper.tri(sigma_e)] <- t(sigma_e)[upper.tri(sigma_e)]
  sigma_u[upper.tri(sigma_u)] <- t(sigma_u)[upper.tri(sigma_u)]
  sigma_e_inv <- solve(sigma_e)
  sigma_u_inv <- solve(sigma_u)

  temp <- y_u - x_covs %*% coeffs
  sigma_e_deriv <- - y_nums * sigma_e_inv +
    sigma_e_inv %*%  t(temp) %*% temp %*% sigma_e_inv
  sigma_u_deriv <- - u_nums * sigma_u_inv +
    sigma_u_inv %*% t(u_all) %*% u_all %*% sigma_u_inv
  diag(sigma_u_deriv) <- diag(sigma_u_deriv) / 2
  c(t(x_covs) %*% temp %*% sigma_e_inv,
    sigma_e_deriv[lower.tri(sigma_e_deriv, diag = F)], sigma_u_deriv[lower.tri(sigma_u_deriv, diag = T)])
}

sigma_hessian <- function(sigma_inv, residu, calcu_diag){
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

      Vij <- matrix(0, K,K) #I_tmp1 * I_tmp2
      V_upperij <- sigma_inv %*% (I_tmp1 %*% sigma_inv %*% I_tmp2 + I_tmp2 %*% sigma_inv %*% I_tmp1 - Vij) %*% sigma_inv
      h_res[k1,k2] <- h_res[k2,k1] <-
        sum(diag(- residu_num * sigma_inv %*% Vij +
                   residu_num * sigma_inv %*% I_tmp1 %*% sigma_inv %*% I_tmp2 -
                   V_upperij %*% t(residu) %*% residu))/2
    }
  }
  h_res
}
beta_sigma_deriv <- function(sigma_inv, residu, x_covs){
  K <- nrow(sigma_inv)
  off_num <- K*(K-1)/2
  cov_num <- ncol(x_covs)

  deriv_res <- matrix(0, cov_num * K, off_num)

  temp_mat <- t(x_covs) %*% residu %*% sigma_inv
  for(k1 in 1:off_num){
    I_tmp1 <- matrix(0, K, K)
    I_tmp1[lower.tri(I_tmp1, diag = F)][k1] <- 1
    I_tmp1[upper.tri(I_tmp1)] <- t(I_tmp1)[upper.tri(I_tmp1)]
    deriv_res[,k1] <- c(-temp_mat %*% I_tmp1 %*% sigma_inv)
  }
  deriv_res
}

#' @export log_lik_1fac_hessian
log_lik_1fac_hessian <- function(coeffs, sigma_e_inv, sigma_u_inv, x_covs, y_u, u_all) {
  K <- nrow(sigma_e_inv)
  h1 <- - kronecker(sigma_e_inv, t(x_covs) %*% x_covs)
  h3 <- sigma_hessian(sigma_u_inv, u_all, calcu_diag = T)
  if(K==1){
    return(as.matrix(bdiag(h1,h3)))
  } else{
    h2 <- sigma_hessian(sigma_e_inv, y_u - x_covs%*%coeffs, calcu_diag = F)
    h12 <- beta_sigma_deriv(sigma_e_inv, y_u - x_covs%*%coeffs, x_covs)
    hes_mat <- as.matrix(bdiag(rbind(cbind(h1, h12),cbind(t(h12), h2)), h3))
    return(hes_mat)
  }
}
#' @export hiera_1fac_fisher_info
hiera_1fac_fisher_info <- function(Y, x_covs, i_ind, stem_res, burnin, max_steps) {
  K <- ncol(Y)
  rcd_num <- nrow(Y)
  ind_num <- max(i_ind)
  cov_num <- ncol(x_covs)

  t_len <- aggregate(i_ind, by=list(i_ind), length)$x

  ## set estimated parameters
  params_hat <- res_summary(stem_res, burnin)
  coeffs <- params_hat$coeffs
  sigma_e <- params_hat$Sigma_e
  sigma_u <- params_hat$Sigma_u
  sigma_e_inv <- solve(sigma_e)
  sigma_u_inv <- solve(sigma_u)

  ## initialize random variables
  # Y_star <- stem_res$Y_star
  # U_all <- stem_res$U_all
  U_all <- rmvnorm(ind_num, sigma = sigma_u)
  Y_star <- x_covs %*% coeffs + U_all[i_ind,] + rmvnorm(rcd_num, sigma = sigma_e)

  params <- c(coeffs,
              sigma_e[lower.tri(sigma_e, diag = F)],
              sigma_u[lower.tri(sigma_u, diag = T)])
  params_num <- length(params)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  Xbeta <- x_covs %*% coeffs

  info_mat1 <- info_mat2 <- matrix(0, max_steps, params_num * params_num)
  info_mat3 <- matrix(0, max_steps, params_num)
  # store_res <- matrix(0, max_steps / 100, params_num * params_num)
  for(iter in 1:max_steps){
    cat('\r iter: ', iter, 'mu=',coeffs[1,])
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, sigma_e,
                                 Xbeta + U_all[i_ind,])
    cat(' sample Y done |')
    # sample U
    U_all <- sample_X_all1_3fac(U_all, sigma_e_inv, sigma_u_inv, t_len,
                                rowsum(Y_star - Xbeta, i_ind, reorder = T))
    cat('sample U done |')

    # caculate information parts
    g_vals <- log_lik_1fac_deriv(params, x_covs, Y_star - U_all[i_ind, ], U_all)
    # portion_old = (iter - 1.0) / iter
    # portion_new = 1.0 / iter
    # info_mat1 <- portion_old * info_mat1 +
    #   portion_new * jacobian(func=log_lik_1fac_deriv, x=params, x_covs = x_covs,
    #                          y_u = Y_star - U_all[i_ind, ], u_all = U_all, method = "simple")
    # info_mat1 <- portion_old * info_mat1 +
    #   portion_new * log_lik_1fac_hessian(coeffs, sigma_e_inv, sigma_u_inv, x_covs,
    #                                      Y_star - U_all[i_ind, ], U_all)
    # info_mat2 <- portion_old * info_mat2 + portion_new * outer(g_vals, g_vals)
    # info_mat3 <- portion_old * info_mat3 + portion_new * g_vals
    info_mat1[iter,] <- c(log_lik_1fac_hessian(coeffs, sigma_e_inv, sigma_u_inv, x_covs,
                                             Y_star - U_all[i_ind, ], U_all))
    info_mat2[iter,] <- c(outer(g_vals, g_vals))
    info_mat3[iter,] <- g_vals
    # if(!(iter %% 100)){
    #   ## calcuate fisher infomation matrix
    #   temp <- -info_mat1 - info_mat2 + outer(info_mat3, info_mat3)
    #   store_res[iter / 100,] <- c(temp)
    #   cat(round(sqrt(diag(solve(temp))), 3),'\n')
    # }
  }
  return(list('info_mat1' = info_mat1,
              'info_mat2' = info_mat2,
              'info_mat3' = info_mat3))
}

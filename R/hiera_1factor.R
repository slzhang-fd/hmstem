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
              'B_u_all'=B_u_all))
}
#' @export test
test <- function(){
  cat("hello world12!")
}

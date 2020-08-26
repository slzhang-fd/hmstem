## Stochastic EM depend functions

#' @export hiera_2factor
hiera_2factor <- function(Y, x_covs, j_ind, i_ind, max_steps){
  NN <- nrow(Y)
  N <- max(i_ind)
  J <- max(j_ind)
  K <- ncol(Y)
  p <- ncol(x_covs)
  t_len <- aggregate(i_ind, by=list(i_ind), length)$x
  it_len <- aggregate(j_ind, by=list(j_ind), length)$x
  ## initialize random variables
  Y_star <- matrix(0, NN, K)
  U_all <- matrix(rnorm(N*K),N,K)
  V_all <- matrix(rnorm(J*K),J,K)

  ## initialize parameters
  coeffs <- matrix(0, p, K)
  B_e <- B_v <- B_u <- diag(rep(1,K))

  coeffs_all <- matrix(0, max_steps, p*K)
  B_e_all <- B_v_all <- B_u_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  for(iter in 1:max_steps){
    cat('\riter: ', iter, 'mu=',coeffs[1,])
    Sigma_e <- B_e %*% t(B_e)
    Sigma_u <- B_u %*% t(B_u)
    Sigma_v <- B_v %*% t(B_v)
    Sigma_e_inv <- chol2inv(t(B_e))
    Sigma_u_inv <- chol2inv(t(B_u))
    Sigma_v_inv <- chol2inv(t(B_v))
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star(Y_pos_ind, Y_zero_ind, Y_star, U_all, V_all, Xbeta, Sigma_e, j_ind, i_ind)
    cat('\tsample Y done |')
    # sample U
    # U_all <- sample_U_all(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_u_inv, t_cums, j_ind)
    U_all <- sample_U_all1(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_u_inv, t_len, j_ind, i_ind)
    cat('sample U done |')
    # sample V
    # V_all <- sample_V_all(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_v_inv, it_cums, i_ind)
    V_all <- sample_V_all1(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_v_inv, it_len, j_ind, i_ind)
    cat('sample V done \n')
    # update parameters
    params <- Optim_step(x_covs, Y_star, U_all, V_all, coeffs, B_e, Sigma_e_inv, B_u, B_v, j_ind, i_ind)
    coeffs <- params$coeffs
    B_e <- params$B_e
    B_u <- params$B_u
    B_v <- params$B_v
    # store results
    scale_mat <- diag(1/sqrt(rowSums(B_e^2)))
    coeffs_all[iter,] <- coeffs %*% scale_mat
    B_e_all[iter,] <- scale_mat %*% B_e
    B_u_all[iter,] <- scale_mat %*% B_u
    B_v_all[iter,] <- scale_mat %*% B_v
  }
  return(list('coeffs_all'=coeffs_all,
              'B_e_all'=B_e_all,
              'B_u_all'=B_u_all,
              'B_v_all'=B_v_all))
}
#' @export my_in_op
my_in_op <- function(x, tab){
  res <- x %in% tab
  res[which(is.na(x))] <- NA
  return(res)
}


hiera_stem <- function(Y, x_covs, j_ind, i_ind, t_cums, t_len, it_cums, it_len, max_steps){
  NN <- nrow(Y)
  N <- max(i_ind)
  J <- max(j_ind)
  K <- ncol(Y)
  p <- ncol(x_covs)
  ## initialize random variables
  Y_star <- matrix(0, NN, K)
  U_all <- matrix(rnorm(N*K),N,K)
  V_all <- matrix(rnorm(J*K),J,K)

  ## initialize parameters
  coeffs <- matrix(0, p, K)
  B_e <- B_v <- B_u <- diag(rep(1,K))

  coeffs_all <- matrix(0, max_steps, p*K)
  B_e_all <- B_v_all <- B_u_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  for(iter in 1:max_steps){
    cat('iter: ', iter, 'mu=',coeffs[1,])
    Sigma_e <- B_e %*% t(B_e)
    Sigma_u <- B_u %*% t(B_u)
    Sigma_v <- B_v %*% t(B_v)
    Sigma_e_inv <- chol2inv(t(B_e))
    Sigma_u_inv <- chol2inv(t(B_u))
    Sigma_v_inv <- chol2inv(t(B_v))
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star(Y_pos_ind, Y_zero_ind, Y_star, U_all, V_all, Xbeta, Sigma_e, j_ind, i_ind)
    cat('sample Y done |')
    # sample U
    U_all <- sample_U_all(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_u_inv, t_cums, j_ind)
    # U_all <- sample_U_all1(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_u_inv, t_len, j_ind, i_ind)
    cat('sample U done |')
    # sample V
    V_all <- sample_V_all(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_v_inv, it_cums, i_ind)
    # V_all <- sample_V_all1(Y_star, U_all, V_all, Xbeta, Sigma_e_inv, Sigma_v_inv, it_len, j_ind, i_ind)
    cat('sample V done \n')
    # update parameters
    params <- Optim_step(x_covs, Y_star, U_all, V_all, coeffs, B_e, Sigma_e_inv, B_u, B_v, j_ind, i_ind)
    coeffs <- params$coeffs
    B_e <- params$B_e
    B_u <- params$B_u
    B_v <- params$B_v
    # store results
    scale_mat <- diag(1/sqrt(rowSums(B_e^2)))
    coeffs_all[iter,] <- coeffs %*% scale_mat
    B_e_all[iter,] <- scale_mat %*% B_e
    B_u_all[iter,] <- scale_mat %*% B_u
    B_v_all[iter,] <- scale_mat %*% B_v
  }
  return(list('coeffs_all'=coeffs_all,
              'B_e_all'=B_e_all,
              'B_u_all'=B_u_all,
              'B_v_all'=B_v_all))
}




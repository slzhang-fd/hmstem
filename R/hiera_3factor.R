#' @export hiera_3factor
hiera_3factor <- function(Y, x_covs, j_ind, i_ind, jt_ind, max_steps){
  K <- ncol(Y)
  rcd_num <- nrow(Y)
  ind_num <- max(i_ind)
  grp_num <- max(j_ind)
  jt_num <- max(jt_ind)
  cov_num <- ncol(x_covs)

  t_len <- aggregate(i_ind, by=list(i_ind), length)$x
  it_len <- aggregate(j_ind, by=list(j_ind), length)$x
  i_len <- aggregate(jt_ind, by=list(jt_ind), length)$x
  ## initialize random variables
  Y_star <- matrix(0, rcd_num, K)
  U_all <- matrix(rnorm(ind_num*K),ind_num,K)
  V_all <- matrix(rnorm(grp_num*K),grp_num,K)
  W_all <- matrix(rnorm(jt_num*K),jt_num,K)
  ## initialize parameters
  coeffs <- matrix(0, cov_num, K)
  B_e <- B_v <- B_u <- B_w <- diag(rep(1,K))

  coeffs_all <- matrix(0, max_steps, cov_num*K)
  B_e_all <- B_v_all <- B_u_all <- B_w_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  for(iter in 1:max_steps){
    cat('\r iter: ', iter, 'mu=',coeffs[1,])
    Sigma_e <- B_e %*% t(B_e)
    Sigma_u <- B_u %*% t(B_u)
    Sigma_v <- B_v %*% t(B_v)
    Sigma_w <- B_w %*% t(B_w)
    Sigma_e_inv <- chol2inv(t(B_e))
    Sigma_u_inv <- chol2inv(t(B_u))
    Sigma_v_inv <- chol2inv(t(B_v))
    Sigma_w_inv <- chol2inv(t(B_w))
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e,
                                 Xbeta + V_all[j_ind,] + U_all[i_ind,] + W_all[jt_ind,])
    cat(' sample Y done |')
    # sample U
    U_all <- sample_X_all1_3fac(U_all, Sigma_e_inv, Sigma_u_inv, t_len,
                                rowsum(Y_star - Xbeta - V_all[j_ind,] - W_all[jt_ind,], i_ind, reorder = T))
    cat('sample U done |')
    # sample V
    V_all <- sample_X_all1_3fac(V_all, Sigma_e_inv, Sigma_v_inv, it_len,
                                rowsum(Y_star - Xbeta - U_all[i_ind,] - W_all[jt_ind,], j_ind, reorder = T))
    cat('sample V done |')
    # sample W
    W_all <- sample_X_all1_3fac(W_all, Sigma_e_inv, Sigma_w_inv, i_len,
                                rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T))
    cat('sample W done \n')
    # update parameters
    params <- Optim_step_3fac(x_covs, Y_star, U_all, V_all, W_all, coeffs,
                              B_e, Sigma_e_inv, B_u, B_v, B_w, j_ind, i_ind, jt_ind)
    coeffs <- params$coeffs
    B_e <- params$B_e
    B_u <- params$B_u
    B_v <- params$B_v
    B_w <- params$B_w
    # store results
    scale_mat <- diag(1/sqrt(rowSums(B_e^2)))
    coeffs_all[iter,] <- coeffs %*% scale_mat
    B_e_all[iter,] <- scale_mat %*% B_e
    B_u_all[iter,] <- scale_mat %*% B_u
    B_v_all[iter,] <- scale_mat %*% B_v
    B_w_all[iter,] <- scale_mat %*% B_w
  }
  return(list('coeffs_all'=coeffs_all,
              'B_e_all'=B_e_all,
              'B_u_all'=B_u_all,
              'B_v_all'=B_v_all,
              'B_w_all'=B_w_all))
}
#' @export
res_summary <- function(hiera_res, burnin){
  K <- sqrt(ncol(hiera_res[[2]]))
  p <- ncol(hiera_res[[1]]) / K
  stem_len <- nrow(hiera_res[[2]])

  coeffs_hat = matrix(colMeans(matrix(hiera_res[[1]][-(1:burnin),],ncol = p*K)),p,K)
  coeffs_sd = matrix(apply(matrix(hiera_res[[1]][-(1:burnin),],ncol = p*K), 2, sd),p,K)
  Sigma_e_all <- Sigma_u_all <- Sigma_v_all <- Sigma_w_all <- matrix(0, stem_len - burnin, K*K)
  for(i in (burnin+1):stem_len){
    B_e_tmp <- matrix(hiera_res[[2]][i,],K,K)
    Sigma_e_all[i-burnin,] <- c(B_e_tmp %*% t(B_e_tmp))
    B_u_tmp <- matrix(hiera_res[[3]][i,],K,K)
    Sigma_u_all[i-burnin,] <- c(B_u_tmp %*% t(B_u_tmp))
    if(length(hiera_res)>=4){
      B_v_tmp <- matrix(hiera_res[[4]][i,],K,K)
      Sigma_v_all[i-burnin,] <- c(B_v_tmp %*% t(B_v_tmp))
    }
    if(length(hiera_res)>=5){
      B_w_tmp <- matrix(hiera_res[[5]][i,], K,K)
      Sigma_w_all[i-burnin,] <- c(B_w_tmp %*% t(B_w_tmp))
    }
  }
  if(length(hiera_res)==3){
    return(list('coeffs'= coeffs_hat,
                'coeffs_sd' = coeffs_sd,
                'Sigma_e' = matrix(apply(Sigma_e_all,2,mean),K,K),
                'Sigma_e_sd' = matrix(apply(Sigma_e_all,2,sd),K,K),
                'Sigma_u' = matrix(apply(Sigma_u_all,2,mean),K,K),
                'Sigma_u_sd' = matrix(apply(Sigma_u_all,2,sd),K,K)))
  } else if(length(hiera_res)==4){
    return(list('coeffs'= coeffs_hat,
                'coeffs_sd' = coeffs_sd,
                'Sigma_e' = matrix(apply(Sigma_e_all,2,mean),K,K),
                'Sigma_e_sd' = matrix(apply(Sigma_e_all,2,sd),K,K),
                'Sigma_u' = matrix(apply(Sigma_u_all,2,mean),K,K),
                'Sigma_u_sd' = matrix(apply(Sigma_u_all,2,sd),K,K),
                'Sigma_v' = matrix(apply(Sigma_v_all,2,mean),K,K),
                'Sigma_v_sd' = matrix(apply(Sigma_v_all,2,sd),K,K)))
  } else{
    return(list('coeffs'= coeffs_hat,
                'coeffs_sd' = coeffs_sd,
                'Sigma_e' = matrix(apply(Sigma_e_all,2,mean),K,K),
                'Sigma_e_sd' = matrix(apply(Sigma_e_all,2,sd),K,K),
                'Sigma_u' = matrix(apply(Sigma_u_all,2,mean),K,K),
                'Sigma_u_sd' = matrix(apply(Sigma_u_all,2,sd),K,K),
                'Sigma_v' = matrix(apply(Sigma_v_all,2,mean),K,K),
                'Sigma_v_sd' = matrix(apply(Sigma_v_all,2,sd),K,K),
                'Sigma_w' = matrix(apply(Sigma_w_all,2,mean),K,K),
                'Sigma_w_sd' = matrix(apply(Sigma_w_all,2,sd),K,K)))
  }
}

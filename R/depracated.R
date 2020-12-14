# sigma_hessian_expect <- function(sigma_inv, residu, calcu_diag){
#   K <- nrow(sigma_inv)
#   off_num <- K*(K-1)/2 + calcu_diag * K
#   residu_num <- nrow(residu)
#   h_res <- matrix(0, off_num, off_num)
#   for(k1 in 1:off_num){
#     I_tmp1 <- matrix(0, K, K)
#     I_tmp1[lower.tri(I_tmp1, diag = calcu_diag)][k1] <- 1
#     I_tmp1[upper.tri(I_tmp1)] <- t(I_tmp1)[upper.tri(I_tmp1)]

#     for(k2 in k1:off_num){
#       I_tmp2 <- matrix(0, K, K)
#       I_tmp2[lower.tri(I_tmp2, diag = calcu_diag)][k2] <- 1
#       I_tmp2[upper.tri(I_tmp2)] <- t(I_tmp2)[upper.tri(I_tmp2)]

#       h_res[k1,k2] <- h_res[k2,k1] <-
#         - sum(diag(residu_num * sigma_inv %*% I_tmp1 %*% sigma_inv %*% I_tmp2)) / 2
#     }
#   }
#   h_res
# }
# #' @importFrom numDeriv jacobian
# #' @export log_lik_3fac
# log_lik_3fac <- function(params, x_covs, y_u_v_w) {
#     cov_num <- ncol(x_covs)
#     k_dims <- ncol(y_u_v_w)
#     y_nums <- nrow(y_u_v_w)
#     # v_nums <- nrow(v_all)
#     u_nums <- nrow(u_all)
#     # w_nums <- nrow(w_all)
#     # decompose params
#     coeffs <- matrix(params[1:(cov_num * k_dims)], cov_num, k_dims)
#     sigma_e <- sigma_v <- sigma_u <- sigma_w <- diag(rep(1, k_dims))
#     lower_num <- k_dims * (k_dims - 1) / 2
#     sigma_e[lower.tri(sigma_e, diag = F)] <-
#         params[cov_num * k_dims + 1:lower_num]
#     # sigma_v[lower.tri(sigma_v, diag = T)] <-
#     #     params[cov_num * k_dims + lower_num + 1:(lower_num + k_dims)]
#     sigma_u[lower.tri(sigma_u, diag = T)] <-
#         params[cov_num * k_dims + 2 * lower_num + k_dims + 1:(lower_num + k_dims)]
#     # sigma_w[lower.tri(sigma_w, diag = T)] <-
#     #     params[cov_num * k_dims + 3 * lower_num + 2* k_dims + 1:(lower_num + k_dims)]
#     sigma_e[upper.tri(sigma_e)] <- t(sigma_e)[upper.tri(sigma_e)]
#     # sigma_v[upper.tri(sigma_v)] <- t(sigma_v)[upper.tri(sigma_v)]
#     sigma_u[upper.tri(sigma_u)] <- t(sigma_u)[upper.tri(sigma_u)]
#     # sigma_w[upper.tri(sigma_w)] <- t(sigma_w)[upper.tri(sigma_w)]
#     temp <- y_u_v_w - x_covs %*% coeffs
#     sigma_e_inv <- solve(sigma_e)
#     # sigma_v_inv <- solve(sigma_v)
#     sigma_u_inv <- solve(sigma_u)
#     # sigma_w_inv <- solve(sigma_w)
#     res <- 0.5 * log(det(sigma_e_inv)) * y_nums -
#         0.5 * sum(diag(sigma_e_inv %*% t(temp) %*% temp)) +
#         0.5 * log(det(sigma_u_inv)) * u_nums -
#         0.5 * sum(diag(sigma_u_inv %*% t(u_all) %*% u_all))
#         # (-0.5 * k_dims * log(2 * pi) +
#         # 0.5 * log(det(sigma_v_inv))) * v_nums -
#         # 0.5 * sum(diag(sigma_v_inv %*% t(v_all) %*% v_all)) +
#         # (-0.5 * k_dims * log(2 * pi) +
#         # 0.5 * log(det(sigma_w_inv))) * w_nums -
#         # 0.5 * sum(diag(sigma_w_inv %*% t(w_all) %*% w_all))
#     res
# }

# #' @export log_lik_3fac_deriv
# log_lik_3fac_deriv <- function(params, x_covs, y_u_v_w) {
#     cov_num <- ncol(x_covs)
#     k_dims <- ncol(y_u_v_w)
#     y_nums <- nrow(y_u_v_w)
#     # decompose params
#     coeffs <- matrix(params[1:(cov_num * k_dims)], cov_num, k_dims)
#     sigma_e <- diag(rep(1, k_dims))
#     lower_num <- k_dims * (k_dims - 1) / 2
#     sigma_e[lower.tri(sigma_e, diag = F)] <-
#         params[cov_num * k_dims + 1:lower_num]
#     sigma_e[upper.tri(sigma_e)] <- t(sigma_e)[upper.tri(sigma_e)]
#     sigma_e_inv <- solve(sigma_e)

#     temp <- y_u_v_w - x_covs %*% coeffs
#     sigma_e_deriv <- - y_nums * sigma_e_inv +
#                     sigma_e_inv %*%  t(temp) %*% temp %*% sigma_e_inv
#     c(t(x_covs) %*% temp %*% sigma_e_inv,
#      sigma_e_deriv[lower.tri(sigma_e_deriv, diag = F)])
# }
#' #' @export hiera_3fac_fisher_info
#' hiera_3fac_fisher_info <- function(Y, x_covs, j_ind, i_ind, jt_ind, stem_res, burnin, max_steps) {
#'   K <- ncol(Y)
#'   rcd_num <- nrow(Y)
#'   ind_num <- max(i_ind)
#'   grp_num <- max(j_ind)
#'   jt_num <- max(jt_ind)
#'   cov_num <- ncol(x_covs)
#'
#'   t_len <- aggregate(i_ind, by=list(i_ind), length)$x
#'   it_len <- aggregate(j_ind, by=list(j_ind), length)$x
#'   i_len <- aggregate(jt_ind, by=list(jt_ind), length)$x
#'
#'   ## initialize random variables
#'   Y_star <- stem_res$Y_star
#'   U_all <- stem_res$U_all
#'   V_all <- stem_res$V_all
#'   W_all <- stem_res$W_all
#'
#'   ## set estimated parameters
#'   params_hat <- res_summary(stem_res, burnin)
#'   coeffs <- params_hat$coeffs
#'   sigma_e <- params_hat$Sigma_e
#'   sigma_u <- params_hat$Sigma_u
#'   sigma_v <- params_hat$Sigma_v
#'   sigma_w <- params_hat$Sigma_w
#'   sigma_e_inv <- solve(sigma_e)
#'   sigma_u_inv <- solve(sigma_u)
#'   sigma_v_inv <- solve(sigma_v)
#'   sigma_w_inv <- solve(sigma_w)
#'   params <- c(coeffs, sigma_e[lower.tri(sigma_e, diag = F)])
#'   params_num <- length(params)
#'   Y_pos_ind <- Y==1
#'   Y_zero_ind <- Y==0
#'   Xbeta <- x_covs %*% coeffs
#'
#'   info_mat1 <- info_mat2 <- matrix(0, params_num, params_num)
#'   info_mat3 <- rep(0, params_num)
#'   for(iter in 1:max_steps){
#'     cat('\r iter: ', iter, 'mu=',coeffs[1,])
#'     ## stochastic E step
#'     # sample Y_star
#'     Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, sigma_e,
#'                                  Xbeta + V_all[j_ind,] + U_all[i_ind,] + W_all[jt_ind,])
#'     cat(' sample Y done |')
#'     # sample U
#'     U_all <- sample_X_all1_3fac(U_all, sigma_e_inv, sigma_u_inv, t_len,
#'                                 rowsum(Y_star - Xbeta - V_all[j_ind,] - W_all[jt_ind,], i_ind, reorder = T))
#'     cat('sample U done |')
#'     # sample V
#'     V_all <- sample_X_all1_3fac(V_all, sigma_e_inv, sigma_v_inv, it_len,
#'                                 rowsum(Y_star - Xbeta - U_all[i_ind,] - W_all[jt_ind,], j_ind, reorder = T))
#'     cat('sample V done |')
#'     # sample W
#'     W_all <- sample_X_all1_3fac(W_all, sigma_e_inv, sigma_w_inv, i_len,
#'                                 rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T))
#'     cat('sample W done \n')
#'
#'     # caculate information parts
#'     g_vals <- log_lik_3fac_deriv(params, x_covs, Y_star - U_all[i_ind, ] - V_all[j_ind, ] - W_all[jt_ind, ])
#'     info_mat1 <- info_mat1 - jacobian(func=log_lik_3fac_deriv, x=params, x_covs = x_covs,
#'                                       y_u_v_w = Y_star - U_all[i_ind, ] - V_all[j_ind, ] - W_all[jt_ind, ], method = "simple") / max_steps
#'     info_mat2 <- info_mat2 + outer(g_vals, g_vals) / max_steps
#'     info_mat3 <- info_mat3 + g_vals / max_steps
#'   }
#'   return(list('full_info' = info_mat1,
#'               'cond_info' = info_mat2 - outer(info_mat3, info_mat3)))
#' }


######
#' #' @export log_lik_3fac
#' log_lik_3fac <- function(params, x_covs, y_u_v_w, u_all, v_all, w_all) {
#'   cov_num <- ncol(x_covs)
#'   k_dims <- ncol(y_u_v_w)
#'   y_nums <- nrow(y_u_v_w)
#'   u_nums <- nrow(u_all)
#'   v_nums <- nrow(v_all)
#'   w_nums <- nrow(w_all)
#'
#'   # decompose params
#'   coeffs <- matrix(params[1:(cov_num * k_dims)], cov_num, k_dims)
#'   sigma_e <- sigma_u <- sigma_v <- sigma_w <- diag(rep(1, k_dims))
#'   lower_num <- k_dims * (k_dims - 1) / 2
#'   sigma_e[lower.tri(sigma_e, diag = F)] <-
#'     params[cov_num * k_dims + 1:lower_num]
#'   sigma_u[lower.tri(sigma_u, diag = T)] <-
#'     params[cov_num * k_dims + lower_num + 1:(lower_num + k_dims)]
#'   sigma_v[lower.tri(sigma_v, diag = T)] <-
#'     params[cov_num * k_dims + 2*lower_num+k_dims + 1:(lower_num + k_dims)]
#'   sigma_w[lower.tri(sigma_w, diag = T)] <-
#'     params[cov_num * k_dims + 3*lower_num+2*k_dims + 1:(lower_num+k_dims)]
#'
#'   sigma_e[upper.tri(sigma_e)] <- t(sigma_e)[upper.tri(sigma_e)]
#'   sigma_u[upper.tri(sigma_u)] <- t(sigma_u)[upper.tri(sigma_u)]
#'   sigma_v[upper.tri(sigma_v)] <- t(sigma_v)[upper.tri(sigma_v)]
#'   sigma_w[upper.tri(sigma_w)] <- t(sigma_w)[upper.tri(sigma_w)]
#'   temp <- y_u_v_w - x_covs %*% coeffs
#'   sigma_e_inv <- solve(sigma_e)
#'   sigma_u_inv <- solve(sigma_u)
#'   sigma_v_inv <- solve(sigma_v)
#'   sigma_w_inv <- solve(sigma_w)
#'   res <- 0.5 * log(det(sigma_e_inv)) * y_nums -
#'     0.5 * sum(diag(sigma_e_inv %*% t(temp) %*% temp)) +
#'     0.5 * log(det(sigma_u_inv)) * u_nums -
#'     0.5 * sum(diag(sigma_u_inv %*% t(u_all) %*% u_all)) +
#'     0.5 * log(det(sigma_v_inv)) * v_nums -
#'     0.5 * sum(diag(sigma_v_inv %*% t(v_all) %*% v_all)) +
#'     0.5 * log(det(sigma_w_inv)) * w_nums -
#'     0.5 * sum(diag(sigma_w_inv %*% t(w_all) %*% w_all))
#'   res
#' }
#' #' @export log_lik_3fac_deriv
#' log_lik_3fac_deriv <- function(params, x_covs, y_u_v_w, u_all, v_all, w_all) {
#'   cov_num <- ncol(x_covs)
#'   k_dims <- ncol(y_u_v_w)
#'   y_nums <- nrow(y_u_v_w)
#'   u_nums <- nrow(u_all)
#'   v_nums <- nrow(v_all)
#'   w_nums <- nrow(w_all)
#'   # decompose params
#'   coeffs <- matrix(params[1:(cov_num * k_dims)], cov_num, k_dims)
#'   sigma_e <- sigma_u <- sigma_v <- sigma_w <- diag(rep(1, k_dims))
#'   lower_num <- k_dims * (k_dims - 1) / 2
#'   sigma_e[lower.tri(sigma_e, diag = F)] <-
#'     params[cov_num * k_dims + 1:lower_num]
#'   sigma_u[lower.tri(sigma_u, diag = T)] <-
#'     params[cov_num * k_dims + lower_num + 1:(lower_num+k_dims)]
#'   sigma_v[lower.tri(sigma_v, diag = T)] <-
#'     params[cov_num * k_dims + 2*lower_num+k_dims + 1:(lower_num + k_dims)]
#'   sigma_w[lower.tri(sigma_w, diag = T)] <-
#'     params[cov_num * k_dims + 3*lower_num+2*k_dims + 1:(lower_num+k_dims)]
#'
#'   sigma_e[upper.tri(sigma_e)] <- t(sigma_e)[upper.tri(sigma_e)]
#'   sigma_u[upper.tri(sigma_u)] <- t(sigma_u)[upper.tri(sigma_u)]
#'   sigma_v[upper.tri(sigma_v)] <- t(sigma_v)[upper.tri(sigma_v)]
#'   sigma_w[upper.tri(sigma_w)] <- t(sigma_w)[upper.tri(sigma_w)]
#'   sigma_e_inv <- solve(sigma_e)
#'   sigma_u_inv <- solve(sigma_u)
#'   sigma_v_inv <- solve(sigma_v)
#'   sigma_w_inv <- solve(sigma_w)
#'
#'   temp <- y_u_v_w - x_covs %*% coeffs
#'   sigma_e_deriv <- - y_nums * sigma_e_inv +
#'     sigma_e_inv %*%  t(temp) %*% temp %*% sigma_e_inv
#'   sigma_u_deriv <- - u_nums * sigma_u_inv +
#'     sigma_u_inv %*% t(u_all) %*% u_all %*% sigma_u_inv
#'   sigma_v_deriv <- - v_nums * sigma_v_inv +
#'     sigma_v_inv %*% t(v_all) %*% v_all %*% sigma_v_inv
#'   sigma_w_deriv <- - w_nums * sigma_w_inv +
#'     sigma_w_inv %*% t(w_all) %*% w_all %*% sigma_w_inv
#'   diag(sigma_u_deriv) <- diag(sigma_u_deriv) / 2
#'   diag(sigma_v_deriv) <- diag(sigma_v_deriv) / 2
#'   diag(sigma_w_deriv) <- diag(sigma_w_deriv) / 2
#'   c(t(x_covs) %*% temp %*% sigma_e_inv,
#'     sigma_e_deriv[lower.tri(sigma_e_deriv, diag = F)],
#'     sigma_u_deriv[lower.tri(sigma_u_deriv, diag = T)],
#'     sigma_v_deriv[lower.tri(sigma_v_deriv, diag = T)],
#'     sigma_w_deriv[lower.tri(sigma_w_deriv, diag = T)])
#' }
#'
#' #' @export log_lik_3fac_deriv1
#' log_lik_3fac_deriv1 <- function(sigma_e_inv, sigma_u_inv, sigma_v_inv, sigma_w_inv,
#'                                 x_covs, e_all, u_all, v_all, w_all) {
#'   cov_num <- ncol(x_covs)
#'   k_dims <- ncol(e_all)
#'   e_nums <- nrow(e_all)
#'   u_nums <- nrow(u_all)
#'   v_nums <- nrow(v_all)
#'   w_nums <- nrow(w_all)
#'
#'   sigma_e_deriv <- - e_nums * sigma_e_inv +
#'     sigma_e_inv %*%  t(e_all) %*% e_all %*% sigma_e_inv
#'   sigma_u_deriv <- - u_nums * sigma_u_inv +
#'     sigma_u_inv %*% t(u_all) %*% u_all %*% sigma_u_inv
#'   sigma_v_deriv <- - v_nums * sigma_v_inv +
#'     sigma_v_inv %*% t(v_all) %*% v_all %*% sigma_v_inv
#'   sigma_w_deriv <- - w_nums * sigma_w_inv +
#'     sigma_w_inv %*% t(w_all) %*% w_all %*% sigma_w_inv
#'   diag(sigma_u_deriv) <- diag(sigma_u_deriv) / 2
#'   diag(sigma_v_deriv) <- diag(sigma_v_deriv) / 2
#'   diag(sigma_w_deriv) <- diag(sigma_w_deriv) / 2
#'   c(t(x_covs) %*% e_all %*% sigma_e_inv,
#'     sigma_e_deriv[lower.tri(sigma_e_deriv, diag = F)],
#'     sigma_u_deriv[lower.tri(sigma_u_deriv, diag = T)],
#'     sigma_v_deriv[lower.tri(sigma_v_deriv, diag = T)],
#'     sigma_w_deriv[lower.tri(sigma_w_deriv, diag = T)])
#' }
#' #' @export log_lik_3fac_hessian
#' log_lik_3fac_hessian <- function(sigma_e_inv, sigma_u_inv, sigma_v_inv, sigma_w_inv,
#'                                  x_covs, e_all, u_all, v_all, w_all) {
#'   K <- nrow(sigma_e_inv)
#'   h1 <- - kronecker(sigma_e_inv, t(x_covs) %*% x_covs)
#'   h3 <- sigma_hessian(sigma_u_inv, u_all, calcu_diag = T)
#'   h4 <- sigma_hessian(sigma_v_inv, v_all, calcu_diag = T)
#'   h5 <- sigma_hessian(sigma_w_inv, w_all, calcu_diag = T)
#'   if(K==1){
#'     return(as.matrix(bdiag(h1,h3,h4,h5)))
#'   } else{
#'     h2 <- sigma_hessian(sigma_e_inv, e_all, calcu_diag = F)
#'     h12 <- beta_sigma_deriv(sigma_e_inv, e_all, x_covs)
#'     hes_mat <- as.matrix(bdiag(rbind(cbind(h1, h12),cbind(t(h12), h2)), h3, h4, h5))
#'     return(hes_mat)
#'   }
#' }
#' #' @export hiera_3fac_fisher_info
#' hiera_3fac_fisher_info <- function(Y, x_covs, i_ind, j_ind, jt_ind, stem_res, max_steps) {
#'   rcd_num <- nrow(Y)
#'   ind_num <- max(i_ind)
#'   grp_num <- max(j_ind)
#'   jt_num <- max(jt_ind)
#'   cov_num <- ncol(x_covs)
#'   k_dims <- ncol(Y)
#'
#'   t_len <- aggregate(i_ind, by=list(i_ind), length)$x
#'   it_len <- aggregate(j_ind, by=list(j_ind), length)$x
#'   i_len <- aggregate(jt_ind, by=list(jt_ind), length)$x
#'
#'   ## set estimated parameters
#'   coeffs <- stem_res$coeffs
#'   sigma_e <- stem_res$Sigma_e
#'   sigma_u <- stem_res$Sigma_u
#'   sigma_v <- stem_res$Sigma_v
#'   sigma_w <- stem_res$Sigma_w
#'   sigma_e_inv <- solve(sigma_e)
#'   sigma_u_inv <- solve(sigma_u)
#'   sigma_v_inv <- solve(sigma_v)
#'   sigma_w_inv <- solve(sigma_w)
#'
#'   ## initialize random variables
#'   # Y_star <- stem_res$Y_star
#'   # U_all <- stem_res$U_all
#'   U_all <- rmvnorm(ind_num, sigma = sigma_u)
#'   V_all <- rmvnorm(grp_num, sigma = sigma_v)
#'   W_all <- rmvnorm(jt_num, sigma = sigma_w)
#'   Y_star <- x_covs %*% coeffs + U_all[i_ind,] + V_all[j_ind,] + W_all[jt_ind,] + rmvnorm(rcd_num, sigma = sigma_e)
#'
#'   Y_pos_ind <- Y==1
#'   Y_zero_ind <- Y==0
#'   Xbeta <- x_covs %*% coeffs
#'
#'   params_num <- cov_num * k_dims + 4 * k_dims * (k_dims-1) / 2 + 3 *k_dims
#'   h1 <- - kronecker(sigma_e_inv, t(x_covs) %*% x_covs)
#'   h3 <- h4 <- h5 <- matrix(0, k_dims*(k_dims+1)/2, k_dims*(k_dims+1)/2)
#'   h2 <- matrix(0, k_dims*(k_dims-1)/2,k_dims*(k_dims-1)/2)
#'   h12 <- matrix(0, k_dims * cov_num, k_dims*(k_dims-1)/2)
#'   info_mat1 <- info_mat2 <- matrix(0, params_num, params_num)
#'   info_mat3 <- rep(0, params_num)
#'   store_res <- matrix(0, max_steps/1000, params_num * params_num)
#'   for(iter in 1:(max_steps+1000)){
#'     cat('\r iter: ', iter, 'mu=',coeffs[1,])
#'     ## stochastic E step
#'     # sample Y_star
#'     Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, sigma_e,
#'                                  Xbeta + V_all[j_ind,] + U_all[i_ind,] + W_all[jt_ind,])
#'     cat(' sample Y done |')
#'     # sample U
#'     U_all <- sample_X_all1_3fac(U_all, sigma_e_inv, sigma_u_inv, t_len,
#'                                 rowsum(Y_star - Xbeta - V_all[j_ind,] - W_all[jt_ind,], i_ind, reorder = T))
#'     cat('sample U done |')
#'     # sample V
#'     V_all <- sample_X_all1_3fac(V_all, sigma_e_inv, sigma_v_inv, it_len,
#'                                 rowsum(Y_star - Xbeta - U_all[i_ind,] - W_all[jt_ind,], j_ind, reorder = T))
#'     cat('sample V done |')
#'     # sample W
#'     W_all <- sample_X_all1_3fac(W_all, sigma_e_inv, sigma_w_inv, i_len,
#'                                 rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T))
#'     cat('sample W done')
#'     # caculate information parts
#'     if(iter > 1000){
#'       iter_count <- iter - 1000
#'       E_all <- Y_star - U_all[i_ind,] - V_all[j_ind,] - W_all[jt_ind,]
#'       g_vals <- log_lik_3fac_deriv1(sigma_e_inv, sigma_u_inv, sigma_v_inv, sigma_w_inv,
#'                                     x_covs, E_all, U_all, V_all, W_all)
#'       portion_old = (iter_count - 1.0) / iter_count
#'       portion_new = 1.0 / iter_count
#'       # info_mat1 <- portion_old * info_mat1 +
#'       #   portion_new * log_lik_3fac_hessian(sigma_e_inv, sigma_u_inv, sigma_v_inv, sigma_w_inv,
#'       #                                      x_covs, E_all, U_all, V_all, W_all)
#'       h3 <- portion_old * h3 + portion_new * sigma_hessian(sigma_u_inv, U_all, calcu_diag = T)
#'       h4 <- portion_old * h4 + portion_new * sigma_hessian(sigma_v_inv, V_all, calcu_diag = T)
#'       h5 <- portion_old * h5 + portion_new * sigma_hessian(sigma_w_inv, W_all, calcu_diag = T)
#'       if(k_dims > 1){
#'         h2 <- portion_old * h2 + portion_new *  sigma_hessian(sigma_e_inv, E_all, calcu_diag = F)
#'         h12 <- portion_old * h12 + portion_new * beta_sigma_deriv(sigma_e_inv, E_all, x_covs)
#'       }
#'       info_mat2 <- portion_old * info_mat2 + portion_new * outer(g_vals, g_vals)
#'       info_mat3 <- portion_old * info_mat3 + portion_new * g_vals
#'       if(!(iter_count %% 1000)){
#'         ## calcuate fisher infomation matrix
#'         if(k_dims==1){
#'           info_mat1 <- as.matrix(bdiag(h1,h3,h4,h5))
#'         }else{
#'           info_mat1 <- as.matrix(bdiag(rbind(cbind(h1, h12),cbind(t(h12), h2)), h3, h4, h5))
#'         }
#'         temp <- -info_mat1 - info_mat2 + outer(info_mat3, info_mat3)
#'         store_res[iter_count / 1000,] <- c(temp)
#'         cat(round(sqrt(diag(solve(temp))), 4),'\n')
#'       }
#'     }
#'   }
#'   return(list('info_mat1' = info_mat1,
#'               'info_mat2' = info_mat2,
#'               'info_mat3' = info_mat3,
#'               'store_res' = store_res,
#'               'Y_star' = Y_star,
#'               'U_all' = U_all,
#'               'V_all' = V_all,
#'               'W_all' = W_all))
#' }
#'
#' #' @export hiera_3fac_fisher_addon
#' hiera_3fac_fisher_addon <- function(Y, x_covs, i_ind, j_ind, jt_ind, stem_res, fisher_res, max_steps, iter_done) {
#'   rcd_num <- nrow(Y)
#'   ind_num <- max(i_ind)
#'   grp_num <- max(j_ind)
#'   jt_num <- max(jt_ind)
#'   cov_num <- ncol(x_covs)
#'   k_dims <- ncol(Y)
#'
#'   t_len <- aggregate(i_ind, by=list(i_ind), length)$x
#'   it_len <- aggregate(j_ind, by=list(j_ind), length)$x
#'   i_len <- aggregate(jt_ind, by=list(jt_ind), length)$x
#'
#'   ## set estimated parameters
#'   coeffs <- stem_res$coeffs
#'   sigma_e <- stem_res$Sigma_e
#'   sigma_u <- stem_res$Sigma_u
#'   sigma_v <- stem_res$Sigma_v
#'   sigma_w <- stem_res$Sigma_w
#'   sigma_e_inv <- solve(sigma_e)
#'   sigma_u_inv <- solve(sigma_u)
#'   sigma_v_inv <- solve(sigma_v)
#'   sigma_w_inv <- solve(sigma_w)
#'
#'   ## initialize random variables
#'   U_all <- fisher_res$U_all
#'   V_all <- fisher_res$V_all
#'   W_all <- fisher_res$W_all
#'   Y_star <- fisher_res$Y_star
#'
#'   info_mat1 <- fisher_res$info_mat1
#'   info_mat2 <- fisher_res$info_mat2
#'   info_mat3 <- fisher_res$info_mat3
#'   h1 <- info_mat1[1:(cov_num*k_dims),1:(cov_num*k_dims)]
#'   h2 <- info_mat1[cov_num*k_dims + 1:(k_dims*(k_dims-1)/2), cov_num*k_dims + 1:(k_dims*(k_dims-1)/2)]
#'   h3 <- info_mat1[cov_num*k_dims + k_dims*(k_dims-1)/2 + 1:(k_dims*(k_dims+1)/2),
#'                   cov_num*k_dims + k_dims*(k_dims-1)/2 + 1:(k_dims*(k_dims+1)/2)]
#'   h4 <- info_mat1[cov_num*k_dims + k_dims*(k_dims-1)/2 + k_dims*(k_dims+1)/2 + 1:(k_dims*(k_dims+1)/2),
#'                   cov_num*k_dims + k_dims*(k_dims-1)/2 + k_dims*(k_dims+1)/2 + 1:(k_dims*(k_dims+1)/2)]
#'   h5 <- info_mat1[cov_num*k_dims + k_dims*(k_dims-1)/2 + k_dims*(k_dims+1) + 1:(k_dims*(k_dims+1)/2),
#'                   cov_num*k_dims + k_dims*(k_dims-1)/2 + k_dims*(k_dims+1) + 1:(k_dims*(k_dims+1)/2)]
#'
#'   h12 <- info_mat1[1:(cov_num * k_dims), cov_num * k_dims + 1:(k_dims*(k_dims-1)/2)]
#'
#'
#'   Y_pos_ind <- Y==1
#'   Y_zero_ind <- Y==0
#'   Xbeta <- x_covs %*% coeffs
#'
#'   params_num <- cov_num * k_dims + 4 * k_dims * (k_dims-1) / 2 + 3 *k_dims
#'   store_res <- matrix(0, max_steps/1000, params_num * params_num)
#'   for(iter in 1:max_steps){
#'     cat('\r iter: ', iter, 'mu=',coeffs[1,])
#'     ## stochastic E step
#'     # sample Y_star
#'     Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, sigma_e,
#'                                  Xbeta + V_all[j_ind,] + U_all[i_ind,] + W_all[jt_ind,])
#'     cat(' sample Y done |')
#'     U_all <- sample_X_all1_3fac(U_all, sigma_e_inv, sigma_u_inv, t_len,
#'                                 rowsum(Y_star - Xbeta - V_all[j_ind,] - W_all[jt_ind,], i_ind, reorder = T))
#'     cat('sample U done |')
#'     V_all <- sample_X_all1_3fac(V_all, sigma_e_inv, sigma_v_inv, it_len,
#'                                 rowsum(Y_star - Xbeta - U_all[i_ind,] - W_all[jt_ind,], j_ind, reorder = T))
#'     cat('sample V done |')
#'     W_all <- sample_X_all1_3fac(W_all, sigma_e_inv, sigma_w_inv, i_len,
#'                                 rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T))
#'     cat('sample W done')
#'
#'     # caculate information parts
#'     iter_count <- iter_done + iter
#'     E_all <- Y_star - U_all[i_ind,] - V_all[j_ind,] - W_all[jt_ind,]
#'     g_vals <- log_lik_3fac_deriv1(sigma_e_inv, sigma_u_inv, sigma_v_inv, sigma_w_inv,
#'                                   x_covs, E_all, U_all, V_all, W_all)
#'     portion_old = (iter_count - 1.0) / iter_count
#'     portion_new = 1.0 / iter_count
#'     info_mat1 <- portion_old * info_mat1 +
#'       portion_new * log_lik_3fac_hessian(sigma_e_inv, sigma_u_inv, sigma_v_inv, sigma_w_inv,
#'                                          x_covs, E_all, U_all, V_all, W_all)
#'     info_mat2 <- portion_old * info_mat2 + portion_new * outer(g_vals, g_vals)
#'     info_mat3 <- portion_old * info_mat3 + portion_new * g_vals
#'     if(!((iter-1000) %% 1000)){
#'       ## calcuate fisher infomation matrix
#'       temp <- -info_mat1 - info_mat2 + outer(info_mat3, info_mat3)
#'       store_res[iter / 1000,] <- c(temp)
#'       cat(round(sqrt(diag(solve(temp))), 4),'\n')
#'     }
#'   }
#'   return(list('info_mat1' = info_mat1,
#'               'info_mat2' = info_mat2,
#'               'info_mat3' = info_mat3,
#'               'store_res' = store_res,
#'               'Y_star' = Y_star,
#'               'U_all' = U_all,
#'               'V_all' = V_all,
#'               'W_all' = W_all))
#' }

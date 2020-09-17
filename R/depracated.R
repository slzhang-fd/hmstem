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

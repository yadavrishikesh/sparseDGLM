#' Update Fixed Effects Coefficients With Space-Time Covariates
#' This function updates the space-time beta coefficients in the model using posterior mean and covariance based on the input parameters.
#'
#' @param lambda A reshaped matrix of dimension (nt * ns) x 1, originally of dimension nt x ns.
#' @param m.lambda A vector of length (nt * ns), originally a matrix of dimension nt x ns.
#' @param X A reshaped matrix of covariates of dimension (nt * ns) x p, originally an array of dimension nt x ns x p.
#' @param quad.X A matrix representing X'X (quadratic form of X).
#' @param tau2 A scalar representing the variance of the nugget.
#'
#' @return A vector representing the updated beta coefficients (posterior mean).
#' @export
#'
#' @examples
#' \dontrun{
#' lambda <- matrix(rnorm(100), nrow = 100)
#' m.lambda <- matrix(rnorm(100), nrow = 100)
#' X <- matrix(rnorm(500), nrow = 100, ncol = 5)
#' quad.X <- t(X) %*% X
#' tau2 <- 0.1
#' beta <- update_space_time_beta(lambda, m, X, quad.X, tau2)
#' }
update_space_time_beta <- function(lambda, m.lambda, X, quad.X, tau2) {
#  browser()
  p <- ncol(X)
  prior_cov <- diag(1, p)
  prior_prec <- solve(prior_cov)
  post_prec <- (1 / tau2) * quad.X + prior_prec
  post_cov <- solve(post_prec)
  post_mean <- post_cov %*% t(X) %*% c(lambda - m.lambda) / tau2
  return(post_mean)
}

## test this function
# q.t <- 2
# q.s=2
# q.st =2
# nt = 100
# ns=200
#
# Xt.o <- matrix(X[,1,1:2, drop =FALSE], nrow = nt, ncol = q.t)
# Xs.o<- matrix(X[1, , 3:4, drop = FALSE], nrow = ns, ncol = q.s)
# Xst.o <- matrix(X[, ,5:6, drop = FALSE], nrow = nt * ns, ncol = q.st)
# quad.Xst.o <- t(Xst.o) %*% Xst.o
#
# Ft.o<- Ft[1:nt, ]
# Ft.mean<-   matrix(rowSums(Ft.o * sim_M0$theta), nrow = nt, ncol=ns)
#
# comp.cov.spat<- c(Xs.o%*% beta[3:4])
# comp.cov.spat<- rep(comp.cov.spat, each=nt)
# comp.cov.temp<-  c(Xt.o%*%beta[1:2])
# comp.cov.st<- matrix(c(Xst.o%*% beta[5:6]), nrow=nt, ncol = ns)
#
# update_space_time_beta(lambda = sim_M0$lambda,
#                        m.lambda = sim_M0$mu.st + Ft.mean + comp.cov.spat +  comp.cov.temp,
#                        X = Xst.o,
#                        quad = quad.Xst.o,
#                        tau2 = hyper$tau2
#                        )



#' Update Fixed Effects Coefficients With Purely Temporal Covariates
#'
#' This function updates the purely temporal beta coefficients using the posterior mean and covariance matrix.
#'
#' @param p Integer, the number of covariates.
#' @param ns Integer, the number of spatial locations.
#' @param X A matrix of covariates.
#' @param sigma2 A scalar, the variance parameter.
#' @param Sigma.inv A matrix, the inverse of the covariance matrix.
#' @param mean_lat A vector, the mean of latent variables.
#' @param cov_quad A matrix, the quadratic form of covariates.
#' @param var_hyprior A scalar, the hyperprior variance for the beta coefficients.
#'
#' @return A matrix of updated beta coefficients.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- 5
#' ns <- 50
#' X <- matrix(rnorm(500), nrow = ns, ncol = p)
#' sigma2 <- 0.1
#' Sigma.inv <- diag(ns)
#' mean_lat <- rnorm(ns)
#' cov_quad <- t(X) %*% X
#' var_hyprior <- 0.01
#' beta <- update_purely_temp_beta(p, ns, X, sigma2, Sigma.inv, mean_lat, cov_quad, var_hyprior)
#' }
update_purely_temp_beta<- function(
                                   ns,
                                   X,
                                   sigma2,
                                   Sigma.inv,
                                   mean_lat,
                                   cov_quad,
                                   var_hyprior){
  #browser()
  p<- ncol(X)
  sigma.s.inv<- Sigma.inv / sigma2
  quadr_cov<- sum(sigma.s.inv) # t(rep(1,ns)) %*% sigma.s.inv %*% rep(1,ns)
  mult_mean<- (t(rep(1,ns)) %*% sigma.s.inv) %*% t(mean_lat)

  cov.betas<- Matrix::chol2inv(chol(diag(1/var_hyprior, p) + quadr_cov *  cov_quad))
  mean.betas<- cov.betas %*% colSums(c(mult_mean) * X)

  proposals <- mvtnorm::rmvnorm(1, mean=mean.betas, sigma = cov.betas)
  return(t(proposals))
}


# cross_prod_cov_X_t.o <- matrix(0, nrow = q.t, ncol = q.t)
# if(q.t > 0){
#   for(t in 1:nt){
#     cross_prod_cov_X_t.o <- cross_prod_cov_X_t.o +
#       Xt.o[t,] %*% t(Xt.o[t,])
#   }
# }
# update_purely_temp_beta(ns = ns,
#                         X = matrix(X[,1,1:2, drop =FALSE], nrow = nt, ncol = q.t),
#                         sigma2 =  hyper$sigma2,
#                         Sigma.inv = diag(1,ns),
#                         mean_lat = sim_M0$lambda - Ft.mean - comp.cov.spat - comp.cov.st - sim_M0$mu.st,
#                         cov_quad =cross_prod_cov_X_t.o,
#                         var_hyprior = 2
#                         )
#
# beta[1:2]


#' Update Fixed Effects Coefficients With Purely Spatial Covariates
#'
#' This function updates the purely spatial beta coefficients using the posterior mean and covariance matrix for a spatio-temporal model.
#'
#' @param ns Integer, the number of spatial locations.
#' @param nt Integer, the number of temporal points.
#' @param X A matrix of covariates of dimension `ns x p`.
#' @param sigma2 A scalar, the variance parameter.
#' @param Sigma.inv A matrix, the inverse of the covariance matrix.
#' @param mean_lat A matrix, the mean of latent variables.
#' @param var_hyprior A scalar, the hyperprior variance for the beta coefficients.
#'
#' @return A matrix of updated beta coefficients.
#' @export
#'
#' @examples
#' \dontrun{
#' ns <- 50
#' nt <- 100
#' X <- matrix(rnorm(ns * 5), nrow = ns, ncol = 5)
#' sigma2 <- 0.1
#' Sigma.inv <- diag(ns)
#' mean_lat <- matrix(rnorm(ns * nt), nrow = ns, ncol = nt)
#' var_hyprior <- 0.01
#' beta <- update_purely_spat_beta(ns, nt, X, sigma2, Sigma.inv, mean_lat, var_hyprior)
#' }
update_purely_spat_beta<- function(ns,
                                   nt,
                                   X,
                                   sigma2,
                                   Sigma.inv,
                                   mean_lat,
                                   var_hyprior){
  #browser()
  p<- ncol(X)
  sigma.s.inv<- Sigma.inv / sigma2
  quadr_cov<- nt * (t(X)%*% sigma.s.inv %*% X) # t(rep(1,ns)) %*% sigma.s.inv %*% rep(1,ns)

  mult_mean<- t(X) %*% rowSums(sigma.s.inv %*% t(mean_lat))

  cov.betas<- Matrix::chol2inv(chol(diag(1/var_hyprior, p) + quadr_cov))
  mean.betas<- cov.betas %*% mult_mean

  proposals <- mvtnorm::rmvnorm(1, mean=mean.betas, sigma = cov.betas)
  return(t(proposals))
}

## Check for spatial beta
# update_purely_spat_beta(ns = ns,
#                         nt = nt,
#                         X = Xs.o,
#                         sigma2 = hyper$sigma2,
#                         Sigma.inv = diag(1, ns),
#                         mean_lat = sim_M0$lambda - Ft.mean - comp.cov.temp - comp.cov.st - sim_M0$mu.st,
#                         var_hyprior = 2
#                                    )
# beta[3:4]



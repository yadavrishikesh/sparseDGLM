#' Update Fixed Effects Coefficients With Space-Time Covariates
#'
#' This function updates the space-time beta coefficients in the model using posterior mean and covariance based on the input parameters.
#'
#' @param lambda A reshaped matrix of dimension (nt * ns) x 1, originally of dimension nt x ns.
#' @param m A vector of length (nt * ns), originally a matrix of dimension nt x ns.
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
#' m <- matrix(rnorm(100), nrow = 100)
#' X <- matrix(rnorm(500), nrow = 100, ncol = 5)
#' quad.X <- t(X) %*% X
#' tau2 <- 0.1
#' beta <- update_space_time_beta(lambda, m, X, quad.X, tau2)
#' }
update_space_time_beta <- function(lambda, m, X, quad.X, tau2) {
  p <- ncol(X)
  prior_cov <- diag(10, p)
  prior_prec <- solve(prior_cov)
  post_prec <- (1 / tau2) * quad.X + prior_prec
  post_cov <- solve(post_prec)
  post_mean <- post_cov %*% t(X) %*% c(lambda - m) / tau2
  return(post_mean)
}



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
update_purely_temp_beta<- function(p,
                                   ns,
                                   X,
                                   sigma2,
                                   Sigma.inv,
                                   mean_lat,
                                   cov_quad,
                                   var_hyprior){
  # browser()
  sigma.s.inv<- Sigma.inv / sigma2
  quadr_cov<- sum(sigma.s.inv) # t(rep(1,ns)) %*% sigma.s.inv %*% rep(1,ns)
  mult_mean<- (t(rep(1,ns)) %*% sigma.s.inv) %*% t(mean_lat)

  cov.betas<- Matrix::chol2inv(chol(diag(1/var_hyprior, p) + quadr_cov *  cov_quad))
  mean.betas<- cov.betas %*% colSums(c(mult_mean) * X)

  proposals <- mvtnorm::rmvnorm(1, mean=mean.betas, sigma = cov.betas)
  return(t(proposals))
}




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


#' Update the Intercept for the Regression Model
#'
#' This function updates the intercept in a regression model using the posterior distribution, incorporating the likelihood from the observed data and a prior distribution.
#'
#' @param nt Integer, number of temporal points.
#' @param ns Integer, number of spatial locations.
#' @param lambda.mean.diff A matrix of the difference between lambda and the mean values.
#' @param tau2 A scalar, variance parameter of the likelihood.
#' @param prior.var A scalar, variance parameter of the prior distribution for the intercept.
#'
#' @return A scalar representing the updated intercept.
#' @export
#'
#' @examples
#' \dontrun{
#' nt <- 10
#' ns <- 5
#' lambda.mean.diff <- matrix(rnorm(nt * ns), nrow = nt, ncol = ns)
#' tau2 <- 0.1
#' prior.var <- 1
#' intercept <- update_intercept_reg(nt, ns, lambda.mean.diff, tau2, prior.var)
#' }
update_intercept_reg<- function(nt, ns,
                               lambda.mean.diff,
                                tau2, prior.var){
  var<- 1/(nt * ns / tau2 + 1/prior.var)
  mean<- (1/tau2) * sum(lambda.mean.diff) * var
  return(rnorm(1, mean=mean, sd=sqrt(var)))
}


#' Update the Intercept for the Dense Correlation Model
#'
#' This function updates the intercept in a dense correlation model using spatio-temporal data. The update is based on posterior sampling from a Gaussian distribution, accounting for the current and neighboring intercept values.
#'
#' @param nt Integer, number of temporal points.
#' @param ns Integer, number of spatial locations.
#' @param mu_current Matrix of current intercept values with dimensions `nt x ns`.
#' @param mu0 Vector, initial intercept for the first time step.
#' @param lambdas_mean Matrix of mean lambda values with dimensions `nt x ns`.
#' @param sigma.kappa.inv Inverse of the spatial covariance matrix.
#' @param sigma2 Scalar, variance of the spatial covariance function.
#' @param tau2 Scalar, variance parameter for the noise in the model.
#'
#' @return A matrix of updated intercept values of dimensions `nt x ns`.
#' @export
#'
#' @examples
#' \dontrun{
#' nt <- 10
#' ns <- 5
#' mu_current <- matrix(rnorm(nt * ns), nrow = nt, ncol = ns)
#' mu0 <- rep(0, ns)
#' lambdas_mean <- matrix(rnorm(nt * ns), nrow = nt, ncol = ns)
#' sigma.kappa.inv <- diag(ns)
#' sigma2 <- 1
#' tau2 <- 0.5
#' updated_mu <- update_intercept.dense(nt, ns, mu_current, mu0, lambdas_mean, sigma.kappa.inv, sigma2, tau2)
#' }
update_intercept.dense<- function(nt,
                                      ns,
                                      mu_current,
                                      mu0,
                                      lambdas_mean,
                                      sigma.kappa.inv,
                                      sigma2,
                                      tau2){
  #browser()
  A_inv.nt<- tau2^(-1) * diag(ns) + sigma.kappa.inv/sigma2
  chol.latent.cov.inv.nt<- spam::chol(A_inv.nt)
  tchol.latent.cov.inv.nt<-  t(chol.latent.cov.inv.nt)

  A_inv.others <- tau2^(-1) * diag(ns) + 2 * sigma.kappa.inv/sigma2
  chol.latent.cov.inv.other<- spam::chol(A_inv.others)
  tchol.latent.cov.inv.other<- t(chol.latent.cov.inv.other)

  lambdas.part<- tau2^(-1) * lambdas_mean

  #browser()
  for (t in 1:nt) {
    if (t == 1) {
      mu_prev <- rep(0, ns)
      mu_next <- mu_current[t+1, ]
    } else if (t == nt) {
      mu_prev <- mu_current[t - 1, ]
      mu_next <- rep(0, ns)
    } else {
      mu_prev <- mu_current[t - 1, ]
      mu_next <- mu_current[t + 1, ]
    }

    # Posterior parameters
    if (t == 1) {
      chol.latent.cov.inv <- chol.latent.cov.inv.other
      tchol.latent.cov.inv <- tchol.latent.cov.inv.other
      latent.mean.part <- lambdas.part[t,] + (sigma.kappa.inv/sigma2) %*% (mu_next + mu0)
    } else if (t == nt) {
      chol.latent.cov.inv <- chol.latent.cov.inv.nt
      tchol.latent.cov.inv <- tchol.latent.cov.inv.nt
      latent.mean.part <- lambdas.part[t,] + (sigma.kappa.inv/sigma2) %*% mu_prev
    } else {
      chol.latent.cov.inv <- chol.latent.cov.inv.other
      tchol.latent.cov.inv <- tchol.latent.cov.inv.other
      latent.mean.part <- lambdas.part[t,] + (sigma.kappa.inv/sigma2) %*% (mu_prev + mu_next)
    }
    # Sample from the posterior using canonical representations of Gaussian distributions
    omega <- spam::forwardsolve(tchol.latent.cov.inv, (latent.mean.part))
    mm <- spam::backsolve(chol.latent.cov.inv, omega)
    zz <- matrix(rnorm(ns, sd=1), ns, 1)
    vv <- spam::backsolve(chol.latent.cov.inv, zz)
    proposals <- mm + vv

    mu_current[t, ] <- proposals
  }
  return(mu_current)
}



#' Update Intercept for the Sparse Correlation Model
#'
#' This function updates the intercept in a sparse correlation model using spatio-temporal data. The update is based on posterior sampling from a Gaussian distribution, incorporating both current and neighboring intercept values in the sparse matrix structure.
#'
#' @param nt Integer, number of temporal points.
#' @param ns Integer, number of spatial locations.
#' @param nb Integer, number of basis functions used in the sparse representation.
#' @param mu_current Matrix of current intercept values with dimensions `nt x nb`.
#' @param R0 Vector, initial intercept for the first time step.
#' @param lambdas_mean Matrix of mean lambda values with dimensions `nt x ns`.
#' @param Q Matrix, precision matrix for the sparse basis functions.
#' @param t.A.A Matrix, transpose of projection matrix times projection matrix.
#' @param t.A Matrix, transpose of projection matrix.
#' @param sigma2 Scalar, variance of the spatial covariance function.
#' @param tau2 Scalar, variance parameter for the noise in the model.
#'
#' @return A matrix of updated intercept values with dimensions `nt x nb`.
#' @export
#'
#' @examples
#' \dontrun{
#' nt <- 10
#' ns <- 5
#' nb <- 3
#' mu_current <- matrix(rnorm(nt * nb), nrow = nt, ncol = nb)
#' R0 <- rep(0, nb)
#' lambdas_mean <- matrix(rnorm(nt * ns), nrow = nt, ncol = ns)
#' Q <- diag(nb)
#' t.A.A <- diag(nb)
#' t.A <- diag(nb)
#' sigma2 <- 1
#' tau2 <- 0.5
#' updated_mu <- update_intercept.sparse(nt, ns, nb, mu_current, R0, lambdas_mean, Q, t.A.A, t.A, sigma2, tau2)
#' }
update_intercept.sparse<- function(nt,
                                      ns,
                                      nb,
                                      mu_current,
                                      R0,
                                      lambdas_mean,
                                      Q,
                                      t.A.A,
                                      t.A,
                                      sigma2,
                                      tau2
){
  #browser()
  A_inv.nt<- tau2^(-1) * t.A.A + Q/sigma2
  chol.latent.cov.inv.nt<- spam::chol(A_inv.nt)
  tchol.latent.cov.inv.nt<-   Matrix::t(chol.latent.cov.inv.nt)

  A_inv.others <- tau2^(-1) * t.A.A + 2 * Q/sigma2
  chol.latent.cov.inv.other<- spam::chol(A_inv.others)
  tchol.latent.cov.inv.other<-  Matrix::t(chol.latent.cov.inv.other)

  lambdas.part<- tau2^(-1) * (lambdas_mean %*% t.A)

  for (t in 1:nt) {
    if (t == 1) {
      mu_prev <- rep(0, nb)
      mu_next <- mu_current[t+1, ]
    } else if (t == nt) {
      mu_prev <- mu_current[t - 1, ]
      mu_next <- rep(0, nb)
    } else {
      mu_prev <- mu_current[t - 1, ]
      mu_next <- mu_current[t + 1, ]
    }

    # Posterior parameters
    if (t == 1) {
      chol.latent.cov.inv <- chol.latent.cov.inv.other
      tchol.latent.cov.inv <- tchol.latent.cov.inv.other
      latent.mean.part <- lambdas.part[t,] + (Q/sigma2) %*% (mu_next + R0)
    } else if (t == nt) {
      chol.latent.cov.inv <- chol.latent.cov.inv.nt
      tchol.latent.cov.inv <- tchol.latent.cov.inv.nt
      latent.mean.part <- lambdas.part[t,] + (Q/sigma2) %*% mu_prev
    } else {
      chol.latent.cov.inv <- chol.latent.cov.inv.other
      tchol.latent.cov.inv <- tchol.latent.cov.inv.other
      latent.mean.part <- lambdas.part[t,] + (Q/sigma2) %*% (mu_prev + mu_next)
    }
    # Sample from the posterior using canonical representations of Gaussian distributions
    omega <- spam::forwardsolve(tchol.latent.cov.inv, (latent.mean.part))
    mm <- spam::backsolve(chol.latent.cov.inv, omega)
    zz <- matrix(rnorm(nb, sd=1), nb, 1)
    vv <- spam::backsolve(chol.latent.cov.inv, zz)
    proposals <- mm + vv

    mu_current[t, ] <- proposals
  }
  return(mu_current)
}




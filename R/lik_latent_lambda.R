#' Log Likelihood for Lambda with Gaussian Prior
#'
#' This function computes the log likelihood of the latent parameter `lambda` under a Gaussian prior, with mean `mean.lambda` and variance `tau2`.
#'
#' @param lambda A matrix or vector of the latent parameter values.
#' @param mean.lambda A matrix or vector representing the mean of the Gaussian prior.
#' @param tau2 A numeric value representing the variance of the Gaussian prior.
#'
#' @return A matrix or vector of log likelihood values with respect to `lambda`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' lambda <- matrix(rnorm(10), nrow = 5, ncol = 2)
#' mean.lambda <- matrix(0, nrow = 5, ncol = 2)
#' tau2 <- 1
#' log_lik <- lambda_lik_fun(lambda, mean.lambda, tau2)
#' print(log_lik)
#' }
lambda_lik_fun<- function(lambda, mean.lambda, tau2){
  logdens<-  - 0.5/(tau2) * (lambda - mean.lambda)^2
  return(logdens)
}


#' Gradient of the Log Likelihood for Lambda with Gaussian Prior
#'
#' This function computes the gradient of the log likelihood for the latent parameter `lambda` under a Gaussian prior with mean `mean.lambda` and variance `tau2`.
#'
#' @param lambda A matrix or vector of the latent parameter values.
#' @param mean.lambda A matrix or vector representing the mean of the Gaussian prior.
#' @param tau2 A numeric value representing the variance of the Gaussian prior.
#'
#' @return A matrix or vector representing the gradient of the log likelihood with respect to `lambda`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' lambda <- matrix(rnorm(10), nrow = 5, ncol = 2)
#' mean.lambda <- matrix(0, nrow = 5, ncol = 2)
#' tau2 <- 1
#' grad_lik <- lambda_grad_fun(lambda, mean.lambda, tau2)
#' print(grad_lik)
#' }
lambda_grad_fun<- function(lambda, mean.lambda, tau2){
  logdens<- - (1/tau2) * (lambda - mean.lambda)
  return(logdens)
}


#' Hessian of the Log Likelihood for Lambda with Gaussian Prior
#'
#' This function computes the Hessian (second derivative) of the log likelihood for the latent parameter `lambda` under a Gaussian prior with mean `mean.lambda` and variance `tau2`.
#'
#' @param mean.lambda A numeric value representing the mean of the Gaussian prior.
#' @param tau2 A numeric value representing the variance of the Gaussian prior.
#'
#' @return A numeric value representing the Hessian of the log likelihood with respect to `lambda`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' mean.lambda <- 0
#' tau2 <- 1
#' hess_lik <- hess_lambda(mean.lambda, tau2)
#' print(hess_lik)
#' }
hess_lambda<- function(mean.lambda, tau2){
  hess<- -1/tau2
  return(-hess)
}

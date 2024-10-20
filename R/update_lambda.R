#' Update Lambda Using MMALA (Manifold Metropolis-Adjusted Langevin Algorithm)
#'
#' This function updates the lambda parameter using the MMALA method in a spatio-temporal model. It computes the proposal for lambda based on the gradient and Hessian information for different likelihoods (Poisson, Negative Binomial, and Log-Normal).
#'
#' @param ns Integer, number of spatial locations.
#' @param Y Matrix of observed data with dimensions `nt x ns`.
#' @param nt Integer, number of time points.
#' @param cur.lambda Matrix, current values of lambda with dimensions `nt x ns`.
#' @param mean.lambda Matrix, mean of lambda with dimensions `nt x ns`.
#' @param tau2 Scalar, variance of the nugget effect.
#' @param r Scalar, shape parameter for the Negative Binomial distribution.
#' @param k Scalar, precision parameter for the Log-Normal distribution.
#' @param data_lik Character, likelihood type. Options are `"Poisson"`, `"NegB"`, or `"lNormal"`.
#' @param tun_lambda Scalar, tuning parameter for the lambda proposal variance.
#'
#' @return A list with:
#'   - `cur_lambda`: Updated lambda values.
#'   - `ind.acc`: Acceptance indicator for each location and time point.
#' @export
#'
#' @examples
#' \dontrun{
#' ns <- 50
#' nt <- 100
#' Y <- matrix(rpois(ns * nt, lambda = 1), ncol = ns)
#' cur.lambda <- matrix(rnorm(ns * nt), ncol = ns)
#' mean.lambda <- matrix(runif(ns * nt), ncol = ns)
#' tau2 <- 0.5
#' r <- 1
#' k <- 1
#' data_lik <- "Poisson"
#' tun_lambda <- 0.1
#' result <- update_lambda_MMALA(ns, nt, Y, cur.lambda, mean.lambda, tau2, r, k, data_lik, tun_lambda)
#' }
update_lambda_MMALA<- function(ns, nt, Y, cur.lambda, mean.lambda, tau2, r, k, data_lik, tun_lambda){

    neghess<- hess_lambda_data(Y=Y, nt=nt, ns=ns, mean.lambda = mean.lambda, tau2 = tau2, r = r, k = k, data_lik = data_lik) +
    hess_lambda(mean.lambda=mean.lambda, tau2=tau2)

  hess<- 1/neghess
  inv.hess<- neghess

  grad_cur_lambda<-  grad_lambda_data(Y = Y, lambda = cur.lambda, r=r, k=k, data_lik= data_lik) +
    lambda_grad_fun(lambda = cur.lambda, mean.lambda = mean.lambda, tau2=tau2)

  mean_cur<- cur.lambda + 0.5 * tun_lambda * grad_cur_lambda * hess
  prop_lambda<-  mean_cur + matrix(rnorm(n= ns * nt, sd =  sqrt(tun_lambda * hess)), nrow = nt, ncol=ns)

  grad_prop_lambda<- grad_lambda_data(Y = Y, lambda = prop_lambda, r=r, k=k, data_lik= data_lik) +
    lambda_grad_fun(lambda = prop_lambda, mean.lambda = mean.lambda, tau2=tau2)
  mean_prop<- prop_lambda + 0.5 * tun_lambda * grad_prop_lambda * hess

  diff<- (data_lik_fun(Y = Y, lambda = prop_lambda, r = r, k = k, data_lik = data_lik) +
    lambda_lik_fun(lambda  = prop_lambda, mean.lambda = mean.lambda, tau2 = tau2) -
    data_lik_fun(Y = Y, lambda = cur.lambda, r = r, k = k, data_lik = data_lik) -
    lambda_lik_fun(lambda  = cur.lambda, mean.lambda = mean.lambda, tau2 = tau2) -
    (0.5/tun_lambda) * ((cur.lambda - mean_prop) * sqrt(inv.hess))^2 +
    (0.5/tun_lambda) * ((prop_lambda - mean_cur) * sqrt(inv.hess))^2)

  ind.acc<- matrix(log(runif(ns * nt)), nrow = nt, ncol = ns) < diff
  cur.lambda<- ifelse(ind.acc, prop_lambda,  cur.lambda)

  return(list("cur_lambda"= cur.lambda, "ind.acc"=ind.acc))
}

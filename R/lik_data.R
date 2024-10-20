#' Log Likelihood Calculation for Different GLM Families
#'
#' This function computes the log likelihood of the data `Y` given the linear predictor `lambda` for various generalized linear model (GLM) families of distributions. It supports Poisson, Negative Binomial, and Log-Normal distributions.
#'
#' @param Y A matrix representing the observed data, where rows correspond to time points and columns to spatial locations.
#' @param lambda A matrix representing the linear predictor for the GLM, with the same dimensions as `Y`.
#' @param r A numeric parameter for the Negative Binomial distribution, representing the dispersion (shape parameter).
#' @param k A numeric precision parameter for the Log-Normal distribution.
#' @param data_lik A character string specifying the likelihood family. Options include `"Poisson"`, `"NegB"` (Negative Binomial), and `"lNormal"` (Log-Normal).
#'
#' @return A matrix of log likelihood values with the same dimensions as `Y`.
#' @export
#'
#' @examples
#' # Example with Poisson likelihood
#' Y <- matrix(rpois(20, lambda = 3), nrow = 4, ncol = 5)
#' lambda <- log(Y + 1)
#' log_lik_poisson <- data_lik_fun(Y, lambda, r = NULL, k = NULL, data_lik = "Poisson")
#' print(log_lik_poisson)
#'
#' # Example with Negative Binomial likelihood
#' r <- 5
#' log_lik_negb <- data_lik_fun(Y, lambda, r, k = NULL, data_lik = "NegB")
#' print(log_lik_negb)
#'
#' # Example with Log-Normal likelihood
#' k <- 2
#' log_lik_lnormal <- data_lik_fun(Y, lambda, r = NULL, k, data_lik = "lNormal")
#' print(log_lik_lnormal)

data_lik_fun<- function(Y, lambda, r, k, data_lik){
  #browser()
  nt<- nrow(Y)
  ns<- ncol(Y)
  if(data_lik=="Poisson"){
      logdens<- Y * lambda - exp(lambda) # matrix(dpois(Y, lambda=exp(lambda),log=TRUE), nrow = nt, ncol=ns)
  } else if(data_lik=="NegB"){
      logdens<- Y * lambda  - (r + Y) * log(r + exp(lambda)) #matrix(dnbinom(Y, mu=exp(lambda), size= r, log=TRUE), nrow = nt, ncol=ns)
  } else if(data_lik=="lNormal"){
    logdens<- -(0.5/k) * (log(1+Y) - lambda)^2  #matrix(dnorm(log(1+Y), mean= lambda, sd=sqrt(1/k), log = TRUE), nrow=nt, ncol=ns) ### Jacobean 1/(1+Y) not written here as this does not make any difference when using  MH
  }
  return(logdens)
}



#' Gradient of the Log Likelihood with Respect to Lambda
#'
#' This function computes the gradient of the log likelihood with respect to the linear predictor `lambda` for various generalized linear model (GLM) families. It supports Poisson, Negative Binomial, and Log-Normal distributions.
#'
#' @param Y A matrix representing the observed data, where rows correspond to time points and columns to spatial locations.
#' @param lambda A matrix representing the linear predictor for the GLM, with the same dimensions as `Y`.
#' @param r A numeric parameter for the Negative Binomial distribution, representing the dispersion (shape parameter).
#' @param k A numeric precision parameter for the Log-Normal distribution.
#' @param data_lik A character string specifying the likelihood family. Options include `"Poisson"`, `"NegB"` (Negative Binomial), and `"lNormal"` (Log-Normal).
#'
#' @return A matrix of gradient values with respect to `lambda`, with the same dimensions as `Y`.
#' @export
#'
#' @examples
#' # Example with Poisson likelihood
#' Y <- matrix(rpois(20, lambda = 3), nrow = 4, ncol = 5)
#' lambda <- log(Y + 1)
#' grad_poisson <- grad_lambda_data(Y, lambda, r = NULL, k = NULL, data_lik = "Poisson")
#' print(grad_poisson)
#'
#' # Example with Negative Binomial likelihood
#' r <- 5
#' grad_negb <- grad_lambda_data(Y, lambda, r, k = NULL, data_lik = "NegB")
#' print(grad_negb)
#'
#' # Example with Log-Normal likelihood
#' k <- 2
#' grad_lnormal <- grad_lambda_data(Y, lambda, r = NULL, k, data_lik = "lNormal")
#' print(grad_lnormal)
grad_lambda_data<- function(Y, lambda, r, k, data_lik){
  if(data_lik=="Poisson"){
    logdens<- Y - exp(lambda)
  } else if(data_lik=="NegB"){
    logdens<- Y - (Y+r) * exp(lambda)/(r + exp(lambda))
  } else if(data_lik=="lNormal"){
    logdens<-   (1/k) * (log(1+Y) - lambda)
  }
  return(logdens)
}




#' Hessian of the Log-Likelihood with Respect to Lambda for Different GLM Families
#'
#' This function computes the Hessian of the log-likelihood with respect to lambda for different types of generalized linear model (GLM) families such as Poisson, Negative Binomial (NegB), and log-Normal (lNormal).
#'
#' @param Y A matrix of observed data.
#' @param nt Integer, number of temporal points.
#' @param ns Integer, number of spatial locations.
#' @param mean.lambda A matrix of the mean of lambda values.
#' @param tau2 A scalar, variance parameter.
#' @param r A scalar, shape parameter for the Negative Binomial distribution.
#' @param k A scalar, precision parameter for the log-Normal distribution.
#' @param data_lik Character, the type of likelihood. Options include `"Poisson"`, `"NegB"`, and `"lNormal"`.
#'
#' @return A matrix representing the Hessian of the log-likelihood.
#' @export
#'
#' @examples
#' \dontrun{
#' Y <- matrix(rpois(100, lambda=2), nrow=10, ncol=10)
#' nt <- 10
#' ns <- 10
#' mean.lambda <- matrix(log(2), nrow=10, ncol=10)
#' tau2 <- 0.1
#' r <- 1
#' k <- 0.5
#' hessian <- hess_lambda_data(Y, nt, ns, mean.lambda, tau2, r, k, "Poisson")
#' }
hess_lambda_data<- function(Y, nt, ns, mean.lambda, tau2, r, k, data_lik){
  if(data_lik=="Poisson"){
    hess<- - exp(mean.lambda + tau2/2)
  } else if(data_lik=="NegB"){
    hess<- - (Y+r) * r * exp(mean.lambda + tau2/2) / ((r + exp(mean.lambda + tau2/2))^2)
  } else if(data_lik=="lNormal"){
    hess<- -matrix(1/k, ncol = ns, nrow = nt)
  }
  return(-hess)
}










#' Update Variance Hyperparameter Using Inverse Gamma Distribution
#'
#' This function updates the variance parameter by drawing from an inverse gamma posterior distribution, typically used for variance components in Bayesian hierarchical models.
#'
#' @param ns Integer, number of spatial locations.
#' @param nt Integer, number of time points.
#' @param param.hyperprior A vector of length 2, hyperprior parameters for the inverse gamma distribution (shape and rate).
#' @param qud_sum A numeric value, the sum of squared differences or residuals.
#'
#' @return A simulated variance value from the inverse gamma distribution.
#' @export
#'
#' @examples
#' \dontrun{
#' ns <- 5
#' nt <- 100
#' param.hyperprior <- c(0.01, 0.01)
#' qud_sum <- 10.5
#' updated_variance <- update_var(ns, nt, param.hyperprior, qud_sum)
#' }
update_var<- function(ns, nt, param.hyperprior, qud_sum){
  sim<- invgamma::rinvgamma(n=1, shape = 0.5 * ns * nt + param.hyperprior[1],
                            rate = param.hyperprior[2] + 0.5 *  qud_sum)
  return(sim)
}



#' Update Variance of State Vectors (Wt) Using Inverse Gamma Distribution
#'
#' This function updates the variance of state vectors (Wt) in a dynamic linear model by drawing from an inverse gamma posterior distribution, based on the differences between the current and previous state vectors.
#'
#' @param nt Integer, number of time points.
#' @param p Integer, the dimension of the state vector.
#' @param theta0 A vector of length `p`, the initial state vector.
#' @param theta A matrix of dimension `nt x p`, the current state vectors at each time point.
#' @param G A matrix of dimension `p x p`, the state transition matrix.
#' @param param.hyperprior A vector of length 2, hyperprior parameters for the inverse gamma distribution (shape and rate).
#'
#' @return A vector of length `p`, simulated values of the variance for each state vector component from the inverse gamma distribution.
#' @export
#'
#' @examples
#' \dontrun{
#' nt <- 100
#' p <- 5
#' theta0 <- runif(p)
#' theta <- matrix(runif(nt * p), nt, p)
#' G <- diag(p)
#' param.hyperprior <- c(0.01, 0.01)
#' Wt <- update.Wt(nt, p, theta0, theta, G, param.hyperprior)
#' }
update.Wt<- function(
    nt,
    p,
    theta0,
    theta,
    G,
    param.hyperprior
)
{
  #browser()
  quads<- matrix(nrow = nt, ncol = p)
  quads[1,]<- c((theta[1,] - G %*% theta0)^2)
  quads[2:nt,]<-  (theta[2:(nt),] - t(tcrossprod(G, as.matrix(theta[1:(nt-1),]))))^2

  sim<- invgamma::rinvgamma(n=p, shape = 0.5 * nt + param.hyperprior[1],
                            rate = param.hyperprior[2] + 0.5 * colSums(quads))
  return(sim)
}

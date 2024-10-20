#' Update Shape Parameter (r) for Negative Binomial Distribution
#'
#' This function updates the shape parameter (r) for the Negative Binomial distribution using the Metropolis-Hastings algorithm.
#'
#' @param Y Observed data, a matrix.
#' @param lambda Mean parameter for the Negative Binomial distribution.
#' @param r Current value of the shape parameter.
#' @param tun.r Tuning parameter for the proposal distribution of r.
#' @param ind Logical, indicates whether the proposal was accepted.
#'
#' @return A list with:
#'   - `r`: Updated value of the shape parameter.
#'   - `ind`: Acceptance indicator.
#' @export
#'
#' @examples
#' \dontrun{
#' Y <- matrix(rnbinom(100, mu = 10, size = 5), ncol = 10)
#' lambda <- matrix(log(10), ncol = 10)
#' r <- 5
#' tun.r <- 0.1
#' result <- update_r(Y, lambda, r, tun.r)
#' }
update_r<- function(Y, lambda, r, tun.r, ind=FALSE){

  tran.r<- transfo(par = r, lb= 0 , ub=Inf)
  r.prop.tran<- rnorm(n=1, mean = tran.r, sd=sqrt(tun.r))
  r.prop<- inv_transfo(tpar=r.prop.tran, lb= 0, ub=Inf)

  prio.diff<- dexp(r.prop, rate = 1/100, log = TRUE) + jac_inv_transfo(r.prop.tran, lb = 0, ub = Inf, log = TRUE) -
    jac_inv_transfo(tran.r, lb = 0, ub= Inf, log = TRUE) - dexp(r, rate = 1/100, log = TRUE)

 log.lik.diff<-  sum(dnbinom(Y, mu=exp(lambda), size= r.prop, log=TRUE))  -
   sum(dnbinom(Y, mu=exp(lambda), size= r, log=TRUE))

 log.diff<- log.lik.diff + prio.diff

 if(log(runif(1)) < log.diff){
   r = r.prop
   ind=TRUE
 }
 return(list(r=r,
             ind=ind))

}



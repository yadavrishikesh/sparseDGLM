#' Updating SPDE Basis Weight Parameters Using Gibbs Sampling
#'
#' This function updates the SPDE basis weight parameters using Gibbs sampling. It relies on the precision matrix `Q`, the projection matrix `A` (through its crossproduct `t.A.A`), and the latent mean terms.
#'
#' @param nb Integer, the number of basis functions.
#' @param nt Integer, the number of time points.
#' @param tau2 A scalar, variance of the nugget term.
#' @param Q A matrix of dimension `nb x nb`, the precision matrix of the SPDE basis.
#' @param t.A.A A matrix, crossproduct of the projection matrix that projects the SPDE basis nodes to the observed spatial locations.
#' @param latent.mean.parts A matrix, terms that appear in the mean of the SPDE basis weight parameters.
#'
#' @return A matrix of updated SPDE basis weight parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' nb <- 50
#' nt <- 100
#' tau2 <- 0.1
#' Q <- diag(nb)
#' t.A.A <- matrix(runif(nb * nb), nb, nb)
#' latent.mean.parts <- matrix(runif(nt * nb), nt, nb)
#' spde_weights <- update_spde_weights.M1_sp(nb, nt, tau2, Q, t.A.A, latent.mean.parts)
#' }
update_spde_weights.M1_sp<- function(nb,
                                        nt,
                                        tau2,
                                        Q,
                                        t.A.A,
                                        latent.mean.parts
){
  #browser()
  latent.cov.inv<- Q +  nt * t.A.A /  tau2
  latent.mean.part<- latent.mean.parts
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, (latent.mean.part))
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- matrix(rnorm(nb, sd=1), nb, 1)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(t(proposals))
}

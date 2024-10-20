#' Update the Range Parameter (Kappa) for the Dense Correlation Model
#'
#' This function updates the range parameter (kappa) in a dense correlation model using a Metropolis-Hastings step. The update is based on the current kappa, the spatio-temporal data, and covariance matrix.
#'
#' @param cur.kappa Scalar, current value of the range parameter (kappa).
#' @param nt Integer, number of temporal points.
#' @param sigma2 Scalar, variance of the spatial covariance function.
#' @param log.det Scalar, log determinant of the covariance matrix.
#' @param delta Scalar, tuning parameter for the proposal distribution.
#' @param tun.kappa Scalar, tuning parameter for adjusting the proposal variance.
#' @param dist.mat Matrix, distance matrix between spatial locations.
#' @param mu Matrix of current intercept values with dimensions `nt x ns`.
#' @param mu0 Vector, prior mean for the intercept.
#' @param cov_inv Matrix, inverse of the covariance matrix.
#' @param cor.type Character, type of correlation function to use (e.g., "exponential").
#' @param ind Logical, whether to use independent updating.
#'
#' @return Updated value of the range parameter (kappa).
#' @export
#'
#' @examples
#' \dontrun{
#' cur.kappa <- 0.5
#' nt <- 100
#' sigma2 <- 1
#' log.det <- 0
#' delta <- 0.01
#' tun.kappa <- 0.1
#' dist.mat <- matrix(runif(25), nrow = 5)
#' mu <- matrix(rnorm(500), nrow = nt, ncol = 5)
#' mu0 <- rep(0, 5)
#' cov_inv <- diag(5)
#' cor.type <- "exponential"
#' updated_kappa <- update_kappa.dense(nt, cur.kappa, mu, mu0, sigma2, log.det, cov_inv, delta, tun.kappa, dist.mat, cor.type)
#' }
update_kappa.dense<- function(nt, cur.kappa, mu, mu0, sigma2, log.det, cov_inv, delta, tun.kappa, dist.mat, cor.type, ind=FALSE){
  cur.kappa.logit<- transfo(par = cur.kappa, lb= 0.05 * delta, ub = delta)
  prop.kappa.logit<- rnorm(n=1, mean = cur.kappa.logit, sd=sqrt(tun.kappa))
  prop.kappa<- inv_transfo(tpar = prop.kappa.logit, lb = 0.05 * delta, ub = delta)

  jaco_diff<- jac_inv_transfo(tpar = prop.kappa.logit, lb =0.05 * delta, ub = delta, log = TRUE) -
    jac_inv_transfo(tpar = cur.kappa.logit, lb = 0.05 * delta, ub = delta, log = TRUE)

  inv_covlog.det<- det_inv_cov_matrix(cov_mat =   cor.fun(cor.type = cor.type, dist.mat = dist.mat, kappa = prop.kappa),
                                      det = TRUE)
  prop_log_det<- inv_covlog.det$detr
  cov_inv_prop<- inv_covlog.det$inv

  diff_mus<- mu[2:nt,] - mu[1:(nt-1),]
  # diff_lik<- -0.5 * nt * (prop_log_det  - log.det) -
  #   (0.5 /sigma2) * ((t(mu[1,] - mu0) %*% cov_inv_prop %*% (mu[1,] - mu0)) - (t(mu[1,]- mu0) %*% cov_inv %*% (mu[1,] - mu0))) -
  #   (0.5 / sigma2) * (sum((diff_mus %*% cov_inv_prop) * diff_mus) - sum((diff_mus %*% cov_inv) * diff_mus))  +
  #   jaco_diff
  #
  diff_lik<- (-0.5 * (nt-1) * (prop_log_det  - log.det) -
    (0.5 / sigma2) * (sum((diff_mus %*% cov_inv_prop) * diff_mus) - sum((diff_mus %*% cov_inv) * diff_mus))  +
    jaco_diff)

  if(log(runif(1)) < diff_lik){
    cur.kappa <- prop.kappa
    log.det<- prop_log_det
    cov_inv <- cov_inv_prop
    ind<- TRUE
  }
  return(list("cur.kappa" = cur.kappa,
              "log.det" = log.det,
              "cov_inv" = cov_inv,
              "ind" = ind
  )
  )
}



#' Update Kappa using Metropolis-Hastings for Sparse Model
#'
#' This function updates the range parameter (kappa) for a sparse spatial model using a Metropolis-Hastings algorithm. It considers the sparse precision matrix structure and spatial covariance.
#'
#' @param cur.kappa Scalar, current value of the range parameter (kappa).
#' @param sigma2 Scalar, variance of the spatial covariance function.
#' @param c.mat Matrix, constant matrix in the SPDE precision formulation.
#' @param g1.mat Matrix, part of the first-order spatial derivative matrix.
#' @param g2.mat Matrix, part of the second-order spatial derivative matrix.
#' @param cor.mat.inv Matrix, inverse of the current correlation matrix.
#' @param log.det.cor.mat Scalar, log determinant of the correlation matrix.
#' @param delta Scalar, tuning parameter for the proposal distribution.
#' @param tun.kappa Scalar, tuning parameter for adjusting the proposal variance.
#' @param ind Logical, whether to use independent updating.
#' @param nt Integer, number of temporal points.
#' @param R Matrix, current values of the latent process.
#' @param R0 Matrix, prior values of the latent process at time `t = 0`.
#'
#' @return Updated value of the range parameter (kappa).
#' @export
#'
#' @examples
#' \dontrun{
#' cur.kappa <- 0.5
#' sigma2 <- 1
#' c.mat <- matrix(runif(25), nrow = 5)
#' g1.mat <- matrix(runif(25), nrow = 5)
#' g2.mat <- matrix(runif(25), nrow = 5)
#' cor.mat.inv <- diag(5)
#' log.det.cor.mat <- 0
#' delta <- 0.01
#' tun.kappa <- 0.1
#' nt <- 100
#' R <- matrix(rnorm(500), nrow = nt, ncol = 5)
#' R0 <- rep(0, 5)
#' updated_kappa <- update_kappa.sparse(nt, cur.kappa, sigma2, c.mat, g1.mat, g2.mat, R, R0, cor.mat.inv, log.det.cor.mat, delta, tun.kappa)
#' }
update_kappa.sparse<- function(nt, cur.kappa, sigma2, c.mat, g1.mat, g2.mat,
                                R, R0, cor.mat.inv, log.det.cor.mat, delta, tun.kappa, ind=FALSE){
  # browser()
  tran.kappa<- transfo(par = cur.kappa, lb=0.05 * delta, ub=delta)
  prop.tran.kappa<- rnorm(n=1, mean = tran.kappa, sd=sqrt(tun.kappa))
  prop.kappa<- inv_transfo(tpar=prop.tran.kappa, lb=0.05 * delta, ub=delta)

  cor.mat.inv.and.log.det.prop<- cormat.inv.update.inla(rho = prop.kappa, c.mat=c.mat,
                                                        g1.mat=g1.mat, g2.mat=g2.mat, alpha = 2)

  # cor.mat.inv.and.log.det.prop<- covmat.inv.update.inla(rho = prop.kappa,  mar.std = sqrt(sigma2), c.mat=c.mat,
  #                                                       g1.mat=g1.mat, g2.mat=g2.mat, alpha = 2)
  cor.mat.inv.prop<- cor.mat.inv.and.log.det.prop$cormat.inv
  log.det.cor.mat.prop<- cor.mat.inv.and.log.det.prop$cormat.logdet

  diff_mus<- R[2:nt,] - R[1:(nt-1),]
  # log.diff.lik<- - 0.5 * nt * (log.det.cor.mat.prop - log.det.cor.mat) -
  #   (0.5) * sum((t(R[1,] - R0) %*% cor.mat.inv.prop %*% (R[1,] - R0)) - (t(R[1,]- R0) %*% cor.mat.inv %*% (R[1,] - R0))) -
  #   (0.5) * (sum((diff_mus %*% cor.mat.inv.prop) * diff_mus) - sum((diff_mus %*% cor.mat.inv) * diff_mus))


  log.diff.lik<- - 0.5 * nt * (log.det.cor.mat.prop - log.det.cor.mat) -
    (0.5 /sigma2) * sum((t(R[1,] - R0) %*% cor.mat.inv.prop %*% (R[1,] - R0)) - (t(R[1,]- R0) %*% cor.mat.inv %*% (R[1,] - R0))) -
    (0.5/sigma2) * (sum((diff_mus %*% cor.mat.inv.prop) * diff_mus) - sum((diff_mus %*% cor.mat.inv) * diff_mus))


  prio.diff<- jac_inv_transfo(prop.tran.kappa, lb = 0.05 * delta, ub = delta, log = TRUE) -
    jac_inv_transfo(tran.kappa, lb = 0.05 * delta, ub= delta, log = TRUE)

  log.diff<- log.diff.lik + prio.diff
  if(log(runif(1)) < log.diff){
    cur.kappa = prop.kappa
    log.det.cor.mat = log.det.cor.mat.prop
    cor.mat.inv = cor.mat.inv.prop
    ind=TRUE
  }

  return(list(cur.kappa=cur.kappa,
              log.det.cor.mat=log.det.cor.mat,
              cor.mat.inv=cor.mat.inv,
              ind=ind))
}






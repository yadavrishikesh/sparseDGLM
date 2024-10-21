#' Inverse of a Matrix using the Sherman-Morrison-Woodbury Formula
#'
#' This function computes the inverse of a matrix using the Sherman-Morrison-Woodbury (SMW) formula.
#'
#' @param alpha A numeric scalar that represents a parameter involved in the inversion calculation.
#' @param tau A numeric scalar that represents the diagonal elements of the matrix to be inverted.
#' @param ns An integer representing the dimension of the square matrix (n x n).
#'
#' @return A matrix representing the inverse of the original matrix calculated using the SMW formula.
#'
#' @details
#' The function calculates the inverse of a matrix of the form \( A + u v^top \) efficiently using the SMW formula.
#' Here, the matrix \( A \) is diagonal with elements \(tau\), and \( u \) and \( v \) are vectors of 1s.
#'
#' @export
#'
#' @examples
#' alpha <- 0.5
#' tau <- 1
#' ns <- 5
#' inverse_SMW(alpha, tau, ns)
inverse_SMW <- function(alpha, tau, ns) {
  A.inv<- diag(1/tau, ns)
  inv = A.inv  - ((1/tau^2)/(1 / alpha + ns/tau)) * matrix(1,nrow = ns, ncol=ns)
  return(inv)
}


#' Compute the Inverse and Log-Determinant of the Correlation Matrix
#'
#' This function computes the inverse and log-determinant of a correlation matrix based on a given range parameter (`kappa`), a distance matrix (`dist.mat`), and a specified correlation type.
#'
#' @param kappa A numeric value representing the range parameter.
#' @param dist.mat A matrix representing the distance matrix between spatial locations.
#' @param cor.type A character string specifying the type of correlation. This is passed to the `cor.fun` function to determine the correlation structure.
#'
#' @return A list containing:
#' \describe{
#'   \item{cormat.inv}{The inverse of the correlation matrix.}
#'   \item{cormat.logdet}{The log-determinant of the correlation matrix.}
#' }
#'
#' @details
#' The function uses the Cholesky decomposition to efficiently compute the inverse and log-determinant of the correlation matrix. The correlation matrix is constructed using the `cor.fun` function, which takes into account the type of correlation and the distance matrix.
#'
#' @export
#'
#' @examples
#' kappa <- 1
#' dist.mat <- matrix(runif(9), 3, 3)
#' cor.type <- "Matern1.5"
#' result <- cormat.inv.fun(kappa, dist.mat, cor.type)
#' result$cormat.inv  # Inverse correlation matrix
#' result$cormat.logdet  # Log-determinant of the correlation matrix
cormat.inv.fun <- function(kappa, dist.mat, cor.type){
 # browser()
  Sigma<- cor.fun(cor.type = cor.type, dist.mat = dist.mat, kappa = kappa)
  chol.mat<- Matrix::chol(Sigma)
  cormat.inv = Matrix::chol2inv(chol.mat)
  cormat.logdet <- 2 * sum(log(diag(chol.mat))) ## log(det(C)) based on cormat.inv, where C is covariance matrix
  return(list(cormat.inv = cormat.inv, cormat.logdet = cormat.logdet))
}



#' Create Precision Matrix Based on the SPDE Triangulation
#'
#' This function creates the precision matrix for a spatial process based on the Stochastic Partial Differential Equation (SPDE) triangulation. The matrix is constructed using the range parameter (`rho`), marginal standard deviation (`mar.std`), and several matrices (`c.mat`, `g1.mat`, `g2.mat`).
#'
#' @param rho A numeric value representing the range parameter, parametrized such that it represents the distance at which the correlations drop approximately to 0.13.
#' @param mar.std A numeric value for the marginal standard deviation. Default is 1.
#' @param c.mat A matrix involved in the update of the inverse covariance matrix.
#' @param g1.mat A matrix involved in the update of the inverse covariance matrix.
#' @param g2.mat A matrix involved in the update of the inverse covariance matrix.
#' @param alpha A numeric value representing the smoothness parameter, default is 2.
#' @param logDet Logical, if `TRUE`, the log-determinant of the updated covariance matrix will be computed and returned.
#'
#' @return A list containing:
#' \describe{
#'   \item{cormat.inv}{The updated inverse covariance matrix.}
#'   \item{cormat.logdet}{The log-determinant of the covariance matrix if `logDet = TRUE`, otherwise `NULL`.}
#' }
#'
#' @details
#' The precision matrix is calculated using the given range parameter (`rho`), marginal standard deviation (`mar.std`), and the matrices (`c.mat`, `g1.mat`, `g2.mat`) representing the SPDE structure. The function also provides the option to compute the log-determinant of the covariance matrix.
#'
#' @export
#'
#' @examples
#' rho <- 1
#' mar.std <- 1
#' c.mat <- matrix(runif(9), 3, 3)
#' g1.mat <- matrix(runif(9), 3, 3)
#' g2.mat <- matrix(runif(9), 3, 3)
#' result <- covmat.inv.update.inla(rho, mar.std, c.mat, g1.mat, g2.mat)
#' result$cormat.inv  # Inverse covariance matrix
#' result$cormat.logdet  # Log-determinant (if logDet = TRUE)
covmat.inv.update.inla <- function(rho, mar.std = 1, c.mat, g1.mat, g2.mat, alpha = 2, logDet=TRUE){
  kappa= sqrt(8)/rho
  tau2s= 1 / (4 * pi * mar.std^2 * kappa^2)
  cormat.inv <- tau2s * (kappa^4 * c.mat + 2 * kappa^2 * g1.mat + g2.mat)
  if(logDet) {cormat.logdet <- -2 * sum(log(diag(spam::chol(cormat.inv)))) ## log(det(C)) based on cormat.inv, where C is covariance matrix
  return(list(cormat.inv = cormat.inv, cormat.logdet = cormat.logdet))
  } else {
    return(list(cormat.inv = cormat.inv))
  }
}



#' Update the Inverse Correlation Matrix in INLA Framework
#'
#' This function updates the inverse of the correlation matrix used in the INLA (Integrated Nested Laplace Approximation) framework based on given parameters and matrices.
#'
#' @param rho A numeric value representing the range parameter.
#' @param c.mat A matrix involved in the update of the inverse correlation matrix.
#' @param g1.mat A matrix involved in the update of the inverse correlation matrix.
#' @param g2.mat A matrix involved in the update of the inverse correlation matrix.
#' @param alpha A numeric value representing the smoothness parameter, default is 2.
#' @param logDet Logical, if `TRUE`, the log-determinant of the updated correlation matrix will be computed and returned.
#'
#' @return A list containing:
#' \describe{
#'   \item{cormat.inv}{The updated inverse correlation matrix.}
#'   \item{cormat.logdet}{The log-determinant of the correlation matrix if `logDet = TRUE`, otherwise `NULL`.}
#' }
#'
#' @details
#' The function calculates the inverse correlation matrix using the given range parameter (`rho`) and matrices (`c.mat`, `g1.mat`, `g2.mat`).
#' The log-determinant of the matrix is optionally calculated and returned when `logDet = TRUE`.
#'
#' @export
#'
#' @examples
#' rho <- 1
#' c.mat <- matrix(runif(9), 3, 3)
#' g1.mat <- matrix(runif(9), 3, 3)
#' g2.mat <- matrix(runif(9), 3, 3)
#' result <- cormat.inv.update.inla(rho, c.mat, g1.mat, g2.mat)
#' result$cormat.inv  # Inverse correlation matrix
#' result$cormat.logdet  # Log-determinant (if logDet = TRUE)
cormat.inv.update.inla <- function(rho, c.mat, g1.mat, g2.mat, alpha = 2, logDet=TRUE){
 # browser()
  kappa= sqrt(8)/rho
  tau2= 1 / (4 * pi * kappa^2)
  cormat.inv <- tau2 * (kappa^4 * c.mat + 2 * kappa^2 * g1.mat + g2.mat)
  if(logDet) {cormat.logdet <- -2 * sum(log(Matrix::diag(spam::chol(cormat.inv)))) ## log(det(C)) based on cormat.inv, where C is covariance matrix
  return(list(cormat.inv = cormat.inv, cormat.logdet = cormat.logdet))
  } else {
    return(list(cormat.inv = cormat.inv))
  }
}

#' Generate Samples from a Multivariate Normal Distribution with a Given Precision Matrix
#'
#' This function generates samples from a multivariate normal distribution given a precision matrix.
#'
#' @param n Integer, the number of samples to draw.
#' @param mu A numeric vector of length `p`, representing the mean of the multivariate normal distribution.
#' @param Omega A numeric matrix of dimension `p x p`, representing the precision matrix (the inverse of the covariance matrix).
#'
#' @return A numeric matrix of dimension `p x n`, where each column represents a sample from the multivariate normal distribution.
#'
#' @details
#' The function generates `n` samples from a multivariate normal distribution using the Cholesky decomposition of the precision matrix (`Omega`).
#' It first generates standard normal random variables, then transforms them using the Cholesky factor of `Omega`, and finally adds the mean vector (`mu`) to obtain the desired samples.
#'
#' @export
#'
#' @examples
#' # Example with a 3-dimensional multivariate normal distribution
#' mu <- c(0, 0, 0)
#' Omega <- diag(3)  # Precision matrix is the identity matrix
#' samples <- rMVNormP(n = 5, mu = mu, Omega = Omega)
#' print(samples)
rMVNormP <- function(n, mu, Omega) {
  p <- length(mu)
  Z <- matrix(rnorm(p * n), p, n)
  U <- chol(Omega)  # Cholesky decomposition of the precision matrix
  X <- backsolve(U, Z)  # Transform standard normal samples to multivariate normal
  X <- sweep(X, 1, mu, FUN = `+`)  # Add the mean to each sample
  return(X)
}


#' Inverse and Determinant of a Matrix using Cholesky Decomposition
#'
#' This function computes the inverse and determinant of a covariance matrix using Cholesky decomposition for efficiency.
#'
#' @param cov_mat A numeric matrix representing the covariance matrix to be inverted.
#' @param det Logical, if `TRUE`, the determinant of the matrix will be computed and returned; if `FALSE`, only the inverse will be returned.
#'
#' @return A list containing:
#' \describe{
#'   \item{inv}{The inverse of the covariance matrix.}
#'   \item{detr}{The determinant of the matrix if `det = TRUE`, otherwise `NULL`.}
#' }
#'
#' @details
#' The Cholesky decomposition is used to compute both the inverse and the log determinant efficiently.
#' The determinant is computed as the sum of the logarithms of the diagonal elements of the Cholesky factor, multiplied by 2.
#'
#' @export
#'
#' @examples
#' cov_mat <- matrix(c(4, 2, 2, 3), nrow = 2)
#' result <- det_inv_cov_matrix(cov_mat)
#' result$inv  # Inverse of the matrix
#' result$detr  # Determinant of the matrix
#'
#' # If you only need the inverse
#' det_inv_cov_matrix(cov_mat, det = FALSE)
det_inv_cov_matrix <- function(cov_mat, det = TRUE) {
  chol <-  Matrix::chol(cov_mat)
  inv <- Matrix::chol2inv(chol)
  if (det) {
    detr <-  2 * sum(log(diag(chol)))
  } else {
    detr = NULL
  }
  return(list("inv" = inv, "detr" = detr))
}


#' Transform Parameter to Unconstrained Scale
#'
#' This function transforms a scalar parameter from a constrained scale (bounded by `lb` and `ub`) to an unconstrained scale, which is useful in optimization and MCMC algorithms.
#'
#' @param par A numeric scalar representing the parameter to be transformed.
#' @param lb A numeric scalar representing the lower bound of the parameter.
#' @param ub A numeric scalar representing the upper bound of the parameter.
#'
#' @return A numeric value representing the transformed parameter on the unconstrained scale.
#'
#' @details
#' The transformation depends on the bounds (`lb`, `ub`):
#' - If both `lb` and `ub` are infinite, the parameter is returned unchanged.
#' - If `lb` is finite and `ub` is infinite, the log transformation is applied: `log(par - lb)`.
#' - If `lb` is infinite and `ub` is finite, the transformation is `log(ub - par)`.
#' - If both bounds are finite, the transformation is the logit function: `qlogis((par - lb) / (ub - lb))`.
#'
#' @references
#' Section 56 of the Stan Reference manual version 2.9 at \url{https://github.com/stan-dev/stan/releases/download/v2.9.0/stan-reference-2.9.0.pdf}
#'
#' @export
#'
#' @examples
#' transfo(par = 0.5, lb = 0, ub = 1)
#' transfo(par = 2, lb = -Inf, ub = Inf)
#' transfo(par = 1, lb = 0, ub = Inf)
#' transfo(par = 1, lb = -Inf, ub = 5)
transfo <- Vectorize(function(par, lb, ub) {
  stopifnot(length(par) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if (lb == -Inf & ub == Inf) {
    return(par)
  } else if (lb > -Inf & ub == Inf) {
    return(log(par - lb))
  } else if (lb == -Inf & ub < Inf) {
    return(log(ub - par))
  } else if (lb > -Inf & ub < Inf) {
    return(qlogis((par - lb) / (ub - lb)))
  }
}, vectorize.args = c("par", "lb", "ub"))

#' Jacobian of the Transformations Made in \code{transfo}
#'
#' This function computes the Jacobian (or log-Jacobian) of the parameter transformations based on specified bounds (`lb`, `ub`) and returns either the Jacobian or its logarithm depending on the `log` argument.
#'
#' @param tpar The transformed parameter.
#' @param lb The lower bound of the parameter.
#' @param ub The upper bound of the parameter.
#' @param log Logical, if `TRUE`, the log-Jacobian is returned; if `FALSE`, the Jacobian is returned.
#'
#' @return A numeric value representing the Jacobian or log-Jacobian for the transformed parameter.
#'
#' @details
#' The Jacobian depends on the bounds:
#' - If both `lb` and `ub` are infinite, the Jacobian is `log(tpar)`.
#' - If `lb` is finite and `ub` is infinite, the Jacobian is `tpar`.
#' - If `lb` is infinite and `ub` is finite, the Jacobian is `tpar`.
#' - If both bounds are finite, the Jacobian is `log(ub - lb) + logit(tpar) + log(1 - logit(tpar))`.
#'
#' @export
#'
#' @examples
#' jac_inv_transfo(tpar = 0.5, lb = 0, ub = 1, log = TRUE)
#' jac_inv_transfo(tpar = 2, lb = -Inf, ub = Inf, log = FALSE)
#' jac_inv_transfo(tpar = 1, lb = 0, ub = Inf, log = FALSE)
#' jac_inv_transfo(tpar = 1, lb = -Inf, ub = 5, log = TRUE)
jac_inv_transfo <- Vectorize(function(tpar, lb, ub, log = FALSE) {
  if (lb == -Inf & ub == Inf) {
    ljac <- log(tpar)
  }
  if (lb > -Inf & ub == Inf) {
    ljac <- tpar
  } else if (lb == -Inf & ub < Inf) {
    ljac <- tpar
  } else if (lb > -Inf & ub < Inf) {
    ljac <- log(ub - lb) + plogis(tpar, log.p = TRUE) + plogis(tpar, log.p = TRUE, lower.tail = FALSE)
  }
  if (log) {
    return(ljac)
  } else{
    return(exp(ljac))
  }
}, vectorize.args = c("tpar", "lb", "ub"))


#' Gradient of the Jacobian
#'
#' This function computes the gradient of the Jacobian of the transformation, depending on the parameter bounds (`lb` and `ub`).
#'
#' @param tpar The transformed parameter.
#' @param lb The lower bound of the parameter.
#' @param ub The upper bound of the parameter.
#'
#' @return A numeric value representing the gradient of the Jacobian for the transformed parameter.
#'
#' @details
#' The gradient depends on the values of `lb` and `ub`:
#' - If both `lb` and `ub` are infinite, the gradient is 0.
#' - If `lb` is finite and `ub` is infinite, the gradient is 1.
#' - If `lb` is infinite and `ub` is finite, the gradient is 1.
#' - If both bounds are finite, the gradient is calculated as \(-1 + 2 x {plogis}(-tpar)\).
#'
#' @export
#'
#' @examples
#' dlogjac_inv_transfo(tpar = 0.5, lb = 0, ub = 1)
#' dlogjac_inv_transfo(tpar = 2, lb = -Inf, ub = Inf)
#' dlogjac_inv_transfo(tpar = 1, lb = 0, ub = Inf)
#' dlogjac_inv_transfo(tpar = 1, lb = -Inf, ub = 5)
dlogjac_inv_transfo <- Vectorize(function(tpar, lb, ub) {
  if (lb == -Inf & ub == Inf) {
    return(0)
  }
  if (lb > -Inf & ub == Inf) {
    return(1)
  } else if (lb == -Inf & ub < Inf) {
    return(1)
  } else{
    -1 + 2 * plogis(-tpar)
  }
}, vectorize.args = c("tpar", "lb", "ub"))


#' Transform Parameters Back to Original Scales
#'
#' This function transforms a parameter from an unconstrained scale (used in optimization or MCMC) back to its original scale, based on specified lower (`lb`) and upper (`ub`) bounds.
#'
#' @param tpar A numeric value representing the parameter on the unconstrained scale.
#' @param lb A numeric value representing the lower bound of the parameter's original scale.
#' @param ub A numeric value representing the upper bound of the parameter's original scale.
#'
#' @return A numeric value representing the transformed parameter back on the original scale.
#'
#' @details
#' The transformation logic depends on the specified bounds:
#' - If both `lb` and `ub` are infinite, the parameter is unchanged.
#' - If `lb` is finite and `ub` is infinite, an exponential transformation is applied and shifted by `lb`.
#' - If `lb` is infinite and `ub` is finite, the parameter is transformed by subtracting an exponential term from `ub`.
#' - If both `lb` and `ub` are finite, the inverse logit transformation is applied to map the parameter back to the original scale.
#'
#' @export
#'
#' @examples
#' # Case with finite bounds
#' inv_transfo(tpar = 0.5, lb = 0, ub = 1)
#'
#' # Case with lb = -Inf and ub = Inf (no transformation)
#' inv_transfo(tpar = 2, lb = -Inf, ub = Inf)
#'
#' # Case with lb = 0 and ub = Inf (exponential transformation)
#' inv_transfo(tpar = 1, lb = 0, ub = Inf)
#'
#' # Case with lb = -Inf and ub = 5
#' inv_transfo(tpar = 1, lb = -Inf, ub = 5)
inv_transfo <- Vectorize(function(tpar, lb, ub) {
  stopifnot(length(tpar) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if (lb == -Inf & ub == Inf) {
    return(tpar)
  } else if (lb > -Inf & ub == Inf) {
    return(exp(tpar) + lb)
  } else if (lb == -Inf & ub < Inf) {
    return(ub - exp(tpar))
  } else{
    return(lb + (ub - lb) * plogis(tpar))
  }
}, vectorize.args = c("tpar", "lb", "ub"))




#' Adaptive Function for Updating Variance Components in MH and MALA Algorithms
#'
#' This function adapts the tuning parameter for variance updates in the Metropolis-Hastings (MH) and Metropolis-Adjusted Langevin Algorithm (MALA) during MCMC iterations.
#' The adaptation is based on the acceptance rate of the MCMC algorithm.
#'
#' @param index_MCMC_iter Integer, the current iteration index of the MCMC algorithm.
#' @param sigma2_adapt Numeric, the tuning parameter for variance updates in MH and MALA algorithms.
#' @param target_accept Numeric, the target acceptance probability. Typically 0.224 for random walk MH and 0.57 for MALA.
#' @param rate_adapt Numeric, the observed acceptance rate over a fixed number of MCMC iterations.
#' @param burn_in1 Integer, the first burn-in period for the MCMC algorithm.
#' @param burn_in2 Integer, the second burn-in period for the MCMC algorithm.
#' @param adapt Integer, the sequence interval at which `sigma2_adapt` is updated.
#' @param adpat_param Numeric, the adaptation parameter for scaling updates to `sigma2_adapt`.
#' @param adapt_seq Integer vector, specifies the iterations at which adaptation occurs.
#' @param lower.acc Numeric, the lower bound of the acceptance rate range. For MH, it is 0.15; for MALA, it is 0.50.
#' @param upper.acc Numeric, the upper bound of the acceptance rate range. For MH, it is 0.30; for MALA, it is 0.65.
#'
#' @return A numeric value representing the updated `sigma2_adapt` tuning parameter.
#'
#' @details
#' This adaptive mechanism works by adjusting the `sigma2_adapt` parameter to keep the acceptance rate of MCMC proposals within a desirable range. It adjusts more frequently during the burn-in period and less frequently afterward. If the acceptance rate deviates too far from the target range, the function adjusts the variance to improve mixing.
#'
#' @export
#'
#' @examples
#' index_MCMC_iter <- 100
#' sigma2_adapt <- 1.0
#' target_accept <- 0.57
#' rate_adapt <- 0.50
#' burn_in1 <- 500
#' burn_in2 <- 200
#' adapt <- 50
#' adpat_param <- 0.1
#' adapt_seq <- seq(50, 1000, by = 50)
#' lower.acc <- 0.50
#' upper.acc <- 0.65
#' adpative_function(index_MCMC_iter, sigma2_adapt, target_accept, rate_adapt,
#'                   burn_in1, burn_in2, adapt, adpat_param, adapt_seq,
#'                   lower.acc, upper.acc)
adpative_function <- function(index_MCMC_iter,
                              sigma2_adapt,
                              target_accept,
                              rate_adapt,
                              burn_in1,
                              burn_in2,
                              adapt,
                              adpat_param,
                              adapt_seq,
                              lower.acc,
                              upper.acc) {
  if (index_MCMC_iter < burn_in1 + burn_in2) {
    if (index_MCMC_iter %in% adapt_seq) {
      if (index_MCMC_iter < burn_in1) {
        sigma2_adapt <- exp(((rate_adapt / adapt) - target_accept) / adpat_param) * sigma2_adapt
      } else {
        sigma2_adapt <- ifelse((((
          rate_adapt / adapt
        ) > upper.acc) | ((
          rate_adapt / adapt
        ) < lower.acc)), exp(((rate_adapt / adapt) - target_accept
        ) / adpat_param) * sigma2_adapt, sigma2_adapt)
      }
    }
  } else {
    sigma2_adapt <- sigma2_adapt
  }
  return(sigma2_adapt)
}

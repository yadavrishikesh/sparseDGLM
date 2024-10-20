#' Initial Values for Model Parameters
#'
#' This function generates initial values for all the parameters of the models specified. The initial values are generated based on the provided data, model type, and likelihood. It supports three model types: `bayes.reg`, `dense`, and `sparse`.
#'
#' @param p Integer, total number of dynamic temporal covariates.
#' @param q Integer, total number of fixed covariates.
#' @param ns Integer, number of spatial locations.
#' @param nt Integer, number of time points.
#' @param nb Integer, number of mesh nodes in the case of the SPDE approach.
#' @param Y A matrix representing the original data.
#' @param model A character string specifying the model type. The possible values are `"bayes.reg"`, `"dense"`, or `"sparse"`.
#' @param data_lik A character string specifying the data likelihood. Options include `"NegB"` for negative binomial or `"lNormal"` for log-normal.
#' @param delta A numeric value representing the maximum distance between the observed locations.
#' @param seed An integer seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{init}{A list of initial parameter values based on the model type.}
#'   \item{param.name}{A vector of parameter names for the model.}
#' }
#'
#' @details
#' This function initializes the parameters for three different models: `bayes.reg` (Bayesian regression), `dense` (dense model with space-time varying intercepts), and `sparse` (sparse model with space-time varying intercepts). The initialization process adjusts for the likelihood type and the model specifications, setting up parameters like spatial and temporal covariates, intercepts, and other components based on the data.
#'
#' @export
#'
#' @examples
#' p <- 3
#' q <- 2
#' ns <- 5
#' nt <- 10
#' nb <- 4
#' Y <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' model <- "bayes.reg"
#' data_lik <- "lNormal"
#' delta <- 1
#' seed <- 123
#' result <- init_fun(p, q, ns, nt, nb, Y, model, data_lik, delta, seed)
#' print(result$init)
#' print(result$param.name)
init_fun <- function(p, q, ns, nt, nb, Y, model, data_lik, delta, seed) {
  na.ind <- is.na(Y)
  mean.temp <- apply(Y,
                     MARGIN = 1,
                     FUN = mean,
                     na.rm = T)
  Y.init <- Y
  for (i in 1:nt) {
    Y.init[i, ][na.ind[i, ]] <-   mean.temp[i]
  }
  set.seed(seed)
  if (model == "bayes.reg") {
    init <- list(
      "r" = if(data_lik=="NegB") {runif(1,0,1000)} else{NULL},
      "k" = if(data_lik=="lNormal") {runif(1,1,10)} else{NULL},
      "beta0" =  mean(log(1 + Y), na.rm = TRUE),
      "tau2" = runif(1, 0, 1),
      "beta" = runif(q, -1, 1),
      "lambda_init" =  matrix(runif(ns * nt), nrow = nt, ncol = ns) + log(1 + Y.init)
    )
    param_names <-  c(
      if(data_lik=="NegB"){expression(r)} else{NULL},
      if(data_lik=="lNormal") {expression(k)} else{NULL},
      expression(beta[0]),
      expression(tau ^ 2),
      paste0("beta", 1:q),
      expression(lambda[1]),
      expression(lambda[2]),
      expression(lambda[3]),
      expression(lambda[4])
    )
  } else if (model =="dense") {  ## dense model with space-time varying intercepts
    init <- list(
      "r" = if(data_lik=="NegB") {runif(1,0,1000)} else{NULL},
      "k" = if(data_lik=="lNormal") {runif(1,1,10)} else{NULL},
      "kappa" = runif(1, 0.1 * delta, delta),
      "sigma2" = runif(1, 0, 1),
      "tau2" = runif(1, 0, 1),
      "Wt" = runif(p, 0, 1),
      "mu.st" = log(1 + Y.init),
      "beta" = runif(q, -1, 1),
      "theta" = matrix(runif(nt * p, -1, 1), ncol = p, nrow = nt),
      "lambda_init" = 0.5 * matrix(runif(ns * nt), nrow = nt, ncol = ns) + log(1 +
                                                                                 Y.init)
    )
    param_names <-  c(
      if(data_lik=="NegB"){expression(r)} else{NULL},
      if(data_lik=="lNormal") {expression(k)} else{NULL},
      expression(kappa),
      expression(sigma ^ 2),
      expression(tau ^ 2),
      paste0("sigma_omega", 1:p),
      "intercept[1]",
      "intercept[2]",
      "intercept[3]",
      paste0("beta", 1:q),
      expression(theta[1]),
      expression(theta[2]),
      expression(theta[3]),
      expression(theta[4]),
      expression(lambda[1]),
      expression(lambda[2]),
      expression(lambda[3]),
      expression(lambda[4])
    )
  } else if (model =="sparse") {
    ## sparse model with space-time varying intercepts
    init <- list(
      "r" = if(data_lik=="NegB") {runif(1,0,1000)} else{NULL},
      "k" = if(data_lik=="lNormal") {runif(1,1,10)} else{NULL},
      "kappa" = runif(1, 0.1 * delta, delta),
      "sigma2" = runif(1, 0, 1),
      "tau2" = runif(1, 0, 1),
      "Wt" = runif(p, 0, 1),
      "R.st" = matrix(
        rnorm(nb * nt, mean = 0, sd = sqrt(3)),
        ncol = nb,
        nrow = nt
      ),
      "beta" = runif(q, -1, 1),
      "theta" = matrix(runif(nt * p, -1, 1), ncol = p, nrow = nt),
      "lambda_init" = 0.5 * matrix(runif(ns * nt), nrow = nt, ncol = ns) + log(1 +
                                                                                 Y.init)
    )
    param_names <-  c(
      if(data_lik=="NegB"){expression(r)} else{NULL},
      if(data_lik=="lNormal") {expression(k)} else{NULL},
      expression(kappa),
      expression(sigma ^ 2),
      expression(tau ^ 2),
      paste0("sigma_omega", 1:p),
      "R[1]",
      "R[1]",
      "R[2]",
      paste0("beta", 1:q),
      expression(theta[1]),
      expression(theta[2]),
      expression(theta[3]),
      expression(theta[4]),
      expression(lambda[1]),
      expression(lambda[2]),
      expression(lambda[3]),
      expression(lambda[4])
    )
  }
  return(list(init = init, param.name = param_names))
}

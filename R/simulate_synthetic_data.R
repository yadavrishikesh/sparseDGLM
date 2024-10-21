#' Simulate Data from the DGLM Model
#'
#' This function simulates data from a Dynamic Generalized Linear Model (DGLM) with different response variables and space-time structures in the predictors. The function supports multiple model types (`bayes.reg`, `dense`, and `sparse`) and allows for simulation from various likelihoods (`Poisson`, `NegB`, `lNormal`).
#'
#' @param nt Integer, number of time points.
#' @param model Character, specifying the model type. Options include `"bayes.reg"`, `"dense"`, and `"sparse"`.
#' @param data_lik Character, specifying the likelihood type. Options include `"Poisson"`, `"NegB"`, and `"lNormal"`.
#' @param loc A matrix of spatial coordinates for each spatial location.
#' @param mesh.inf A list containing mesh parameters (`max.edge`, `cutoff`, and `offset`). If set to `NULL`, the mesh is created automatically inside the function with default values.
#' @param hyper A list containing model-specific hyperparameters such as `kappa`, `sigma2`, `tau2`, and `r`.
#' @param beta0 Numeric, the intercept for the regression model.
#' @param beta A numeric vector of regression coefficients.
#' @param X A 3D array of fixed covariates for the regression model.
#' @param Ft A matrix of dynamic temporal covariates.
#' @param G A matrix for updating dynamic temporal covariates over time.
#' @param cor.type Character, specifying the correlation type for the dense model.
#' @param theta0 A numeric vector of initial values for the dynamic covariates.
#' @param mu0 Numeric, the initial value for the spatial mean.
#' @param C0 Numeric, the initial variance for the spatial model.
#' @param R0 Numeric, the initial value for the spatial random field.
#'
#' @return A list containing simulated data, including the response `y`, linear predictor `lambda`, dynamic covariate `theta`, covariate effects, spatial random effects `mu.st`, and other model-specific outputs.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' nt <- 100
#' ns <- 50
#' loc <- matrix(runif(ns * 2), ncol = 2)
#' hyper <- list(kappa = 0.5, sigma2 = 1, tau2 = 0.1)
#' beta <- rnorm(5)
#' X <- array(rnorm(nt * ns * 5), dim = c(nt, ns, 5))
#' Ft <- matrix(rnorm(nt * 3), nt, 3)
#' G <- diag(3)
#' result <- sim_synthetic_data(nt, "dense", "Poisson", loc, hyper = hyper, beta = beta, X = X, Ft = Ft, G = G)
#' }
sim_synthetic_data<- function(nt, model, data_lik, loc, mesh.inf=NULL,
                              mesh.hyper = NULL,
                              hyper, beta0, beta, Wt,
                              X, Ft, G, cor.type,
                              theta0, mu0, C0, R0
){

 # browser()
  ns<- nrow(loc)
  q<- dim(X)[3]
  X.cov<- matrix(aperm(X, c(1, 2, 3)),  nrow = nt * ns, ncol = q)

  loc_matrix <- as.matrix(loc)
  dist.mat.all <- as.matrix(fields::rdist.earth(loc_matrix, loc_matrix, miles = FALSE))
  min_dist <- min(dist.mat.all[dist.mat.all > 0])
  max_dist <- max(dist.mat.all)

  dist.mat.all<- dist.mat.all/max_dist ## standardizing the  distnace matrix

  range_x <- diff(range(loc_matrix[, 1]))  # Range of longitude
  range_y <- diff(range(loc_matrix[, 2]))  # Range of latitude
  if(is.null( mesh.inf)){
  max.edge <- c(0.05 * max_dist, 0.2 * max_dist)
  cutoff <- 0.01 * min_dist
  offset <- c(0.05 * range_x, 0.1 * range_x)
  } else{
    max.edge = sparse.info.sim$max.edge
    cutoff =  sparse.info.sim$cutoff
    offset = sparse.info.sim$offset
}
  INLA.mesh <- INLA::inla.mesh.2d(loc = loc_matrix,
                             max.edge = max.edge,
                             cutoff = cutoff,
                             offset = offset)
  A.proj<- as.matrix(INLA::inla.spde.make.A(mesh = INLA.mesh, loc = as.matrix(loc)))
  fem.mesh <- INLA::inla.mesh.fem(INLA.mesh, order = 2)

  sparse.info.sim<- list(INLA.mesh=INLA.mesh,
                         A.proj = A.proj, c.mat =  fem.mesh$c0,
                         g1.mat = fem.mesh$g1, g2.mat = fem.mesh$g2
  )


  if(is.null(mesh.hyper)){
    mesh.hyper<- list(max.edge =  max.edge, cutoff =cutoff, offset=offset)
  }

  nb<- dim(A.proj)[2]

  ### structures into the linear predictors lambda according to different modeling structures
  if(model == "bayes.reg"){
    tau2<- hyper$tau2
    beta0<- beta0
    comp.cov<- matrix(c(X.cov %*% beta), nrow = nt, ncol=ns)

    lambda<- beta0 + comp.cov +
      matrix(rnorm(nt * ns, sd=sqrt(tau2)), nrow = nt, ncol=ns)

    theta = NULL
    Ft.mean = NULL
    comp.cov = NULL
    mu.st = NULL
    R.st = NULL
    corr = NULL

  } else if(model == "dense"){

    kappa<- hyper$kappa
    sigma2<- hyper$sigma2
    Wt<- hyper$Wt
    tau2<- hyper$tau2

    corr<-  cor.fun(cor.type = cor.type, dist.mat = dist.mat.all, kappa = kappa)
    mu.st<- matrix(NA, nrow=nt, ncol=ns)
    mu.st[1,]<-  mvtnorm::rmvnorm(n = 1, mean = rep(mu0, ns), sigma = sigma2 * corr)
    for (t in 2:nt) {
      mu.st[t,]<-  mu.st[t-1,] + mvtnorm::rmvnorm(n=1, sigma = sigma2 * corr)
    }

    #browser()

    theta<- matrix(NA, nrow=nt, ncol=p)
    theta[1,]<- c(G %*% theta0) + c(mvtnorm::rmvnorm(n=1, sigma = Wt))
    for (t in 2:nt) {
      theta[t,]<- c(G %*% theta[t-1,]) + c(mvtnorm::rmvnorm(n=1, sigma = Wt))
    }

    Ft.mean<- matrix(rowSums(Ft * theta), nrow = nt, ncol=ns)
    comp.cov<- matrix(c(X.cov %*% beta), nrow = nt, ncol=ns)

    lambda<- mu.st +  Ft.mean + comp.cov +
      matrix(rnorm(nt * ns, sd=sqrt(tau2)), nrow = nt, ncol=ns)

    R.st = NULL


  } else if(model =="sparse"){

    A.proj<- sparse.info.sim$A.proj
    nb<- dim(A.proj)[2]

    kappa<- hyper$kappa
    sigma2<- hyper$sigma2
    Wt<- hyper$Wt
    tau2<- hyper$tau2

    Q.rho.det <- covmat.inv.update.inla(rho = kappa,  mar.std = sqrt(sigma2), c.mat = sparse.info.sim$c.mat,
                                        g1.mat = sparse.info.sim$g1.mat, g2.mat = sparse.info.sim$g2.mat,
                                        alpha = 2, logDet = FALSE)
    Q.Z0 <- as.matrix(Q.rho.det$cormat.inv)

    R.st<- matrix(NA, nrow=nt, ncol=nb)
    mu.st<- matrix(NA, nrow=nt, ncol=ns)
    R.st[1,]<- c(rMVNormP(n = 1, mu = rep(R0, nb), Omega = Q.Z0))
    mu.st[1,]<- as.numeric(A.proj %*% R.st[1,])
    for (t in 2:nt) {
      R.st[t,]<-  R.st[t-1,] + c(rMVNormP(n = 1, mu = rep(0, nb), Omega = Q.Z0))
      mu.st[t,]<- as.numeric(A.proj %*% R.st[t,])
    }

    theta<- matrix(NA, nrow=nt, ncol=p)
    theta[1,]<- c(G %*% theta0) + c(mvtnorm::rmvnorm(n=1, sigma = Wt))
    for (t in 2:nt) {
      theta[t,]<- c(G %*% theta[t-1,]) + c(mvtnorm::rmvnorm(n=1, sigma = Wt))
    }

    Ft.mean<- matrix(rowSums(Ft * theta), nrow = nt, ncol=ns)
    comp.cov<- matrix(c(X.cov %*% beta), nrow = nt, ncol=ns)

    lambda<- mu.st +  Ft.mean + comp.cov +
      matrix(rnorm(nt * ns, sd=sqrt(tau2)), nrow = nt, ncol=ns)

    corr = NULL
  }

  ### simulate the data according to a particular likelihood
  if(data_lik=="Poisson"){
    y<- matrix(rpois(n=nt * ns, lambda=exp(lambda)), nrow = nt, ncol = ns)
  } else if(data_lik=="NegB"){
    r<- hyper$r
    y<- matrix(rnbinom(n=nt * ns, size  = r, mu=exp(lambda)), nrow = nt, ncol = ns)
  } else if(data_lik=="lNormal"){
    k<- hyper$k
    y<- matrix(exp(rnorm(n=nt * ns, mean = lambda, sd = sqrt(k)))-1, nrow = nt, ncol = ns) ## -1 here in the data, because, we will be defining log-Normal for Y+1
  }

  return(list("y"=y,
              "sparse.info.sim" = sparse.info.sim,
               "mesh.hyper" = mesh.hyper,
              "dist.mat.all" = dist.mat.all,
              "lambda"= lambda,
              "theta" = theta,
              "Ft.mean" = Ft.mean,
              "comp.cov" = comp.cov,
              "mu.st" = mu.st,
              "R.st" = R.st,
              "delta" = max_dist,
              "corr" = corr)
         )
}


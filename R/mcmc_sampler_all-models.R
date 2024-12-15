#' MCMC Sampler for Spatio-Temporal Dynamic Generalized Linear Models (DGLM)
#'
#' This function performs MCMC based inference  and spatiotemporal predictions for spatiotemporal dynamic generalized linear models (st-DGLM), including dense and its sparse model counterparts.
#'
#' @param Y Matrix of response variables (observed and/or predicted).
#' @param X Array of covariates (design matrix).
#' @param loc Matrix of spatial locations.
#' @param num_harmonics Integer, the number of harmonics.
#' @param seas.period Integer, specifying the periodicity of seasonality. For example, for hourly data, daily seasonality should be set to 24; for daily data, weekly seasonality should be set to 7.
#' @param model Character, specifying the model type. Possible values are:
#'   \itemize{
#'     \item `"bayes.reg"`: Bayesian regression model.
#'     \item `"dense"`: Dense spatio-temporal model.
#'     \item `"sparse"`: Sparse spatio-temporal model.
#'   }
#' @param data_lik Character, specifying the likelihood type. Options include:
#'   \itemize{
#'     \item `"Poisson"`: Poisson likelihood.
#'     \item `"NegB"`: Negative Binomial likelihood.
#'     \item `"lNormal"`: Log-normal likelihood.
#'   }
#' @param cor.type Character, specifying the correlation type for the dense model. Options include:
#'   \itemize{
#'     \item `"Matern0.5"`: Matern correlation with smoothness parameter 0.5.
#'     \item `"Matern1"`: Matern correlation with smoothness parameter 1.
#'     \item `"Matern1.5"`: Matern correlation with smoothness parameter 1.5.
#'     \item `"Matern2.5"`: Matern correlation with smoothness parameter 2.5.
#'     \item `"MaternInf"`: Matern correlation with infinite smoothness (Gaussian correlation).
#'   }
#' @param spatInt.ind Vector of indices for spatial interpolation.
#' @param forcast.ind Vector of indices for forecasting.
#' @param N.MCMC Integer, specifying the number of MCMC iterations.
#' @param samples.store Integer, specifying the number of samples to store.
#' @param sparse.info.sim List containing precomputed sparse model information (optional).
#' @param mesh.hyper List containing the information of hyperparameters when creating mesh nodes, specifically on `max.edge`, `cutoff`, and `offset`.
#' @param thin Integer, specifying the thinning interval for MCMC samples.
#' @param print.result Logical, if `TRUE`, results will be printed.
#' @param traceplot Logical, if `TRUE`, trace plots of MCMC samples will be generated.
#' @param true.values List of true parameter values (optional).
#' @param simulation Logical, if `TRUE`, the function runs in simulation mode.
#' @param no_parallel_chain Integer, specifying the number of parallel chains (1 to 4).
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{"summery.hyper"}{A list containing summary statistics of hyperparameters and variance of state variables.}
#'   \item{"summary.fixed.effects.coeff"}{A list containing the overall intercept and coefficient estimates for covariates.}
#'   \item{"predictions"}{A list of predictions including:
#'     \itemize{
#'       \item{"within.sample"}{Within-sample predictions and imputations.}
#'       \item{"spatial"}{Spatial interpolation predictions.}
#'       \item{"forecast"}{Forecasted values.}
#'       \item{"space.time"}{Space-time predictions.}
#'     }
#'   }
#'   \item{"summary.dynamic.temp.coeff"}{Summary of dynamic temporal coefficients.}
#'   \item{"summary.st.intercepot"}{A list containing the SPDE basis weights and spatial-temporal intercept estimates.}
#'   \item{"summary.linear.predictor"}{Summary of the linear predictor values.}
#'   \item{"traceplots.samples"}{Trace plots of MCMC samples.}
#'   \item{"run.time.costly.comt."}{A list with the runtime details for kappa and intercept parameters.}
#'   \item{"N.MCMC"}{The number of MCMC iterations.}
#'   \item{"original.data"}{A list containing:
#'     \itemize{
#'       \item{"data.fit"}{Original data used for fitting.}
#'       \item{"data.all"}{All available data.}
#'     }
#'   }
#'   \item{"forcast.ind"}{Indices for forecasting.}
#'   \item{"spatInt.ind"}{Indices for spatial interpolation.}
#'   \item{"sparse.info"}{A list with sparse model information, including INLA mesh details and projection matrices.}
#'   \item{"nt"}{The number of time points in the model.}
#'   \item{"ns"}{The number of spatial locations in the model.}
#'   \item{"model"}{The type of model used.}
#'   \item{"data_lik"}{The likelihood type used in the model.}
#'   \item{"cor.type"}{The correlation type used for the dense model.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' N.MCMC <- 1000
#' model <- "dense"
#' cor.type <- "Matern1.5"
#' no_parallel_chain <- 2
#' samples.store <- 100
#' thin <- 10
#' adapt <- 100
#' tun_kappa <- 1e-4
#' tun_lambda <- 1
#' Y <- matrix(rnorm(100), nrow = 10)
#' Ft <- matrix(rnorm(30), nrow = 10)
#' G <- diag(3)
#' X <- array(rnorm(1000), dim = c(10, 10, 3))
#' loc <- matrix(runif(20), ncol = 2)
#' spatInt.ind <- 1:2
#' forcast.ind <- 1:2
#' result <- MCMC.sampler.st.DGLM(N.MCMC, model, cor.type, no_parallel_chain,
#'                                samples.store, thin, adapt, tun_kappa, tun_lambda,
#'                                Y = Y, Ft = Ft, G = G, X = X, loc = loc,
#'                                spatInt.ind = spatInt.ind, forcast.ind = forcast.ind)
#' }
MCMC.sampler.st.DGLM<- function(Y,
                                X,
                                loc,
                                num_harmonics,
                                seas.period,
                                model = "dense",
                                data_lik = "Poisson",
                                cor.type = "Matern0.5",
                                spatInt.ind,
                                forcast.ind,
                                N.MCMC = 5000,
                                samples.store = 250,
                                sparse.info.sim = NULL,
                                mesh.hyper = NULL,
                                thin = 10,
                                print.result=TRUE,
                                traceplot=FALSE,
                                true.values=NULL,
                                simulation = FALSE,
                                no_parallel_chain = 1){

  #browser()
 ## some fixed MCMC hyperparameters
  adapt = 100
  tun_kappa = 1e-4
  tun_lambda = 1
  tun_r = 1e-2

  ## harmonics  extraction
  harmo_info<- generate_harmonic_matrices(nt = nrow(Y),
                             num_harmonics = num_harmonics,
                             period = seas.period)
  G<- harmo_info$G.all
  Ft<- harmo_info$Ft.all

  nt<- nrow(Y) - length(forcast.ind)
  ns<- ncol(Y) - length(spatInt.ind)
  q<- dim(X)[3]
  p<- ncol(Ft)

  ### data
  Y.o<- Y[1:nt, -spatInt.ind]
  ind_NA_data<- is.na(Y.o)
  data.log.mean<- mean(log(1+Y.o), na.rm=TRUE)
  data.mean<- mean(Y.o, na.rm=TRUE)
  Y.intpl<- Y[1:nt, spatInt.ind]
  Y.frcast<- Y[-(1:nt), -spatInt.ind]
  Y.st<- Y[-(1:nt), spatInt.ind]

  ### design matrix for fixed coeffiecients
  X.o<- X[1:nt, -spatInt.ind, ]
  X.o<- matrix(aperm(X.o, c(1, 2, 3)), nrow = nt * ns, ncol = q)
  quad.X.o<- t(X.o) %*% X.o
  X.intpl<- X[1:nt, spatInt.ind, ]
  X.intpl<- matrix(aperm(X.intpl, c(1, 2, 3)),
                   nrow = nt * length(spatInt.ind), ncol = q)
  X.frcast<- X[-(1:nt), -spatInt.ind, ]
  X.frcast<- matrix(aperm(X.frcast, c(1, 2, 3)),
                    nrow = length(forcast.ind) * ns, ncol = q)
  X.st<- X[-(1:nt), spatInt.ind, ]
  X.st<- matrix(aperm(X.st, c(1, 2, 3)),
                nrow = length(forcast.ind) * length(spatInt.ind), ncol = q)

  ### purely temporal covariate with dynamic temporal components
  Ft.o<- Ft[1:nt, ]
  F.mat.o<- array(dim=c(nt, ns, p))
  for (i in 1:p) {
    F.mat.o[,,i]<- Ft.o[,i]
  }
  m0.theta = rep(0, p)
  C0.theta = diag(5, p)

  Ft.frcast<- Ft[-(1:nt), ]

  #### components related to space-time intercept
  m0.mu<- ifelse(is.na(Y[1,]),
                 apply(log(1+Y), MARGIN = 2, FUN = function(x) mean(x, na.rm=TRUE)),
                 log(1+Y[1,]))
  C0.mu<- 1

  m0.R = mean(log(1+Y.o), na.rm=TRUE)
  C0.R = 5

  ### spatial location related components
  loc.o<- loc[-spatInt.ind, ]
  loc.p<- loc[spatInt.ind, ]

  ## distance matrices
  dist.mat.all<-  fields::rdist.earth(loc, loc, miles = FALSE)
  min.dist <- min(dist.mat.all[dist.mat.all > 0])
  max.dist <- max(dist.mat.all)

  ### setting the hyper parameters of sparsity in sparse model
  range_x <- diff(range(loc[, 1]))  # Range of longitude
  range_y <- diff(range(loc[, 2]))  # Range of latitude
  if(is.null(mesh.hyper)){
    max.edge <- c(0.05 * max.dist, 0.2 * max.dist)
    cutoff <- 0.01 * min.dist
    offset <- c(0.05 * range_x, 0.1 * range_x)

    mesh.hyper<- list(max.edge = max.edge, cutoff= cutoff, offset=offset)
  } else{
    max.edge = mesh.hyper$max.edge
    cutoff =  mesh.hyper$cutoff
    offset = mesh.hyper$offset
  }

  dist.mat.o<- fields::rdist.earth(loc.o, loc.o, miles = FALSE)
  max.dist<- max(dist.mat.o)
  dist.mat.o <- dist.mat.o / max.dist
  delta<- max(dist.mat.o)

  dist.mat.all <- dist.mat.all / max(dist.mat.all)


  if(model =="sparse" & is.null(sparse.info.sim)){
    INLA.mesh <- INLA::inla.mesh.2d(loc = as.matrix(loc.o),
                                    max.edge = max.edge,
                                    cutoff = cutoff,
                                    offset = offset
                                    )
    A.proj.o<- INLA::inla.spde.make.A(mesh = INLA.mesh, loc = as.matrix(loc.o)) # as.matrix(INLA::inla.spde.make.A(mesh = INLA.mesh, loc = as.matrix(loc.o)))
    A.proj.p<- as.matrix(INLA::inla.spde.make.A(mesh = INLA.mesh, loc = as.matrix(loc.p)))
    fem.mesh <- INLA::inla.mesh.fem(INLA.mesh, order = 2)
    c.mat.o<-  fem.mesh$c0
    g1.mat.o<- fem.mesh$g1
    g2.mat.o<- fem.mesh$g2

    ### a different ways of defining the mesh nodes
    # INLA.mesh <- INLA::inla.mesh.2d(loc = as.matrix(loc),
    #                                 max.edge = max.edge,
    #                                 cutoff = cutoff,
    #                                 offset = offset
    # )
    # A.proj.all.loc<- as.matrix(INLA::inla.spde.make.A(mesh = INLA.mesh, loc = as.matrix(loc)))
    # A.proj.o<- as.matrix(A.proj.all.loc[-spatInt.ind, ])
    # A.proj.p<- as.matrix(A.proj.all.loc[spatInt.ind, ])
    #
    # fem.mesh <- INLA::inla.mesh.fem(INLA.mesh, order = 2)
    # c.mat.o<-  fem.mesh$c0
    # g1.mat.o<- fem.mesh$g1
    # g2.mat.o<- fem.mesh$g2

  }else if(model =="sparse" & simulation & !is.null(sparse.info.sim)) {

    INLA.mesh<- sparse.info.sim$INLA.mesh
    A.proj<- sparse.info.sim$A.proj
    A.proj.o<- A.proj[-spatInt.ind, ]
    A.proj.p<- A.proj[spatInt.ind, ]
    c.mat.o<-  sparse.info.sim$c.mat
    g1.mat.o<- sparse.info.sim$g1.mat
    g2.mat.o<- sparse.info.sim$g2.mat

  } else{
    INLA.mesh<- NULL
    c.mat.o<- NULL
    g1.mat.o<- NULL
    g2.mat.o<- NULL
    A.proj.o<- NULL
    A.proj.p<- NULL
  }

  ### MCMC sampler related components
  burn_in1 = floor(N.MCMC/4)
  burn_in2 = floor(N.MCMC/2)
  adapt_seq<- seq(from=adapt, to=N.MCMC, by=adapt)



  ### run a Bayesian GLM to decide about the beter inital values for covarites ceofficents and linear predictors

  init_values_cov_lin.pred<- parallel::mclapply(1:no_parallel_chain, FUN = function(ii){
    MCMC.sampler_init(model = "bayes.reg",
                           data_lik = data_lik, nt = nt, ns = ns, p = p, q = q, Y.o = Y.o, data.log.mean = data.log.mean,
                           data.mean = data.mean, ind_NA_Y = ind_NA_data, X.o = X.o, quad.X.o = quad.X.o,
                           spatInt.ind = spatInt.ind, forcast.ind = forcast.ind,
                           N.MCMC = 5000, burn_in1 = 3000, burn_in2 = 1000, adapt_seq = adapt_seq, thin = thin,
                           adapt = adapt, tun_r = tun_r, tun_lambda = tun_lambda, print.result = print.result, traceplot = traceplot,
                           true.values = true.values, simulation = simulation, init.seed = 123 + ii
    )
  }, mc.cores = no_parallel_chain)

 # browser()

  if(no_parallel_chain==1){
    init_beta<- init_values_cov_lin.pred[[1]]$sumry.beta
    init_lambda<- init_values_cov_lin.pred[[1]]$pred.lambda.mean
  } else if(no_parallel_chain==2){
    init_beta<- (init_values_cov_lin.pred[[1]]$sumry.beta + init_values_cov_lin.pred[[2]]$sumry.beta)/2
    init_lambda<- (init_values_cov_lin.pred[[1]]$pred.lambda.mean + init_values_cov_lin.pred[[2]]$pred.lambda.mean)/2
  } else if(no_parallel_chain==3){
    init_beta<- (init_values_cov_lin.pred[[1]]$sumry.beta + init_values_cov_lin.pred[[2]]$sumry.beta + init_values_cov_lin.pred[[3]]$sumry.beta)/3
    init_lambda<- (init_values_cov_lin.pred[[1]]$pred.lambda.mean + init_values_cov_lin.pred[[2]]$pred.lambda.mean + init_values_cov_lin.pred[[3]]$pred.lambda.mean)/3
  } else if(no_parallel_chain==4){
    init_beta<- (init_values_cov_lin.pred[[1]]$sumry.beta + init_values_cov_lin.pred[[2]]$sumry.beta + init_values_cov_lin.pred[[3]]$sumry.beta + init_values_cov_lin.pred[[4]]$sumry.beta)/4
    init_lambda<- (init_values_cov_lin.pred[[1]]$pred.lambda.mean + init_values_cov_lin.pred[[2]]$pred.lambda.mean + init_values_cov_lin.pred[[3]]$pred.lambda.mean + init_values_cov_lin.pred[[4]]$pred.lambda.mean)/4
  }
  #### MCMC sampler for different models
if(model=="bayes.reg"){ ## baseline model (Model with only fixed covariates
  parallel_chain<- parallel::mclapply(1:no_parallel_chain, FUN = function(ii){
    MCMC.sampler_model.reg(model = model,
      data_lik = data_lik, nt = nt, ns = ns, p = p, q = q, Y.o = Y.o, data.log.mean = data.log.mean,
      data.mean = data.mean, Y.intpl = Y.intpl, Y.frcast = Y.frcast, Y.st = Y.st, ind_NA_Y = ind_NA_data,
      X.o = X.o, quad.X.o = quad.X.o, X.intpl = X.intpl, X.frcast = X.frcast, X.st = X.st, X.p = X.p,
      spatInt.ind = spatInt.ind, forcast.ind = forcast.ind, samples.store = samples.store,
      N.MCMC = N.MCMC, burn_in1 = burn_in1, burn_in2 = burn_in2, adapt_seq = adapt_seq, thin = thin,
      adapt = adapt, tun_r = tun_r, tun_lambda = tun_lambda, print.result = print.result, traceplot = traceplot,
      true.values = true.values, simulation = simulation, init_beta, init_lambda,  init.seed = 123 + ii
    )
  }, mc.cores = no_parallel_chain)
} else if(model == "dense"){
  parallel_chain<- parallel::mclapply(1:no_parallel_chain, FUN = function(ii){
    MCMC.sampler_model.dense(model = model,
      data_lik = data_lik, nt = nt, ns = ns, p = p, q = q, Y.o = Y.o, data.log.mean = data.log.mean,
      data.mean = data.mean, Y.intpl = Y.intpl, Y.frcast = Y.frcast, Y.st = Y.st, ind_NA_Y = ind_NA_data,
      X.o = X.o, quad.X.o = quad.X.o, X.intpl = X.intpl, X.frcast = X.frcast, X.st = X.st, X.p = X.p,
      Ft.o = Ft.o, F.mat.o = F.mat.o, Ft.frcast = Ft.frcast, m0.mu = m0.mu, C0.mu = C0.mu, m0.theta = m0.theta,
      C0.theta = C0.theta, G = G, dist.mat.all = dist.mat.all, dist.mat.o = dist.mat.o, delta = delta,
      spatInt.ind = spatInt.ind, forcast.ind = forcast.ind, cor.type =  cor.type, samples.store = samples.store,
      N.MCMC = N.MCMC, burn_in1 = burn_in1, burn_in2 = burn_in2, adapt_seq = adapt_seq, thin = thin, adapt = adapt,
      tun_r = tun_r, tun_kappa = tun_kappa, tun_lambda = tun_lambda, print.result = print.result, traceplot = traceplot,
      true.values = true.values, simulation = simulation, init_beta, init_lambda, init.seed =  123 + ii
      )
  }, mc.cores = no_parallel_chain)
} else if(model == "sparse"){
  parallel_chain<- parallel::mclapply(1:no_parallel_chain, FUN = function(ii){
    MCMC.sampler_model.sparse(model = model,
      data_lik = data_lik, nt = nt, ns = ns, p = p, q = q, Y.o = Y.o, data.log.mean = data.log.mean,
      data.mean = data.mean, Y.intpl = Y.intpl, Y.frcast = Y.frcast, Y.st = Y.st, ind_NA_Y = ind_NA_data,
      X.o = X.o, quad.X.o = quad.X.o, X.intpl = X.intpl, X.frcast = X.frcast, X.st = X.st, X.p = X.p,
      Ft.o = Ft.o, F.mat.o = F.mat.o, Ft.frcast = Ft.frcast, m0.R = m0.R, C0.R = C0.R, m0.theta = m0.theta,
      C0.theta = C0.theta, G = G, delta = delta, A.proj.o = A.proj.o, A.proj.p = A.proj.p, c.mat.o = c.mat.o,
      g1.mat.o = g1.mat.o, g2.mat.o = g2.mat.o, spatInt.ind = spatInt.ind, forcast.ind = forcast.ind, samples.store = samples.store,
      N.MCMC = N.MCMC, burn_in1 = burn_in1, burn_in2 = burn_in2, adapt_seq = adapt_seq, thin = thin, adapt = adapt,
      tun_r = tun_r, tun_kappa = tun_kappa, tun_lambda = tun_lambda, print.result = print.result, traceplot = traceplot,
      true.values = true.values, simulation = simulation, init_beta, init_lambda, init.seed =  123 + ii
      )
  }, mc.cores = no_parallel_chain)
}


 # Generating the required results
  ### assessing the performances based on the multiple chains
    if(no_parallel_chain==1){
    parallel_chain[[2]]<- parallel_chain[[3]]<- parallel_chain[[4]]<- parallel_chain[[1]]
  } else if(no_parallel_chain==2){
    parallel_chain[[3]]<- parallel_chain[[1]]
    parallel_chain[[4]]<- parallel_chain[[2]]
  } else if(no_parallel_chain==3){
    parallel_chain[[4]]<- parallel_chain[[1]]
  }
  ### samples of some parameter to check for convergence and mixing properties
  traceplots.samples<- list(chain1=parallel_chain[[1]]$samples.trace.plot,
                            chain2=parallel_chain[[2]]$samples.trace.plot,
                            chain3=parallel_chain[[3]]$samples.trace.plot,
                            chain4=parallel_chain[[4]]$samples.trace.plot)


  ### Summary of model hyperparmeters
  thin_burn_in<- floor((burn_in1 + burn_in2)/thin)
  if(model == "dense" | model =="sparse"){
  if(data_lik=="Poisson"){
    kappa.samples<- sapply(1:4, function(jj) {max.dist * traceplots.samples[[jj]]$samples.kappa[-(1:thin_burn_in)]}) %>% rowMeans
    sigma2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.sigma2[-(1:thin_burn_in)]}) %>% rowMeans
    tau2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.tau2[-(1:thin_burn_in)]}) %>% rowMeans

    sumry.hyper.mat<- cbind(kappa.samples, sigma2.samples, tau2.samples)

    sumry.hyper<- cbind.data.frame(mean = apply(sumry.hyper.mat, MARGIN = 2, FUN = mean),
                                 median = apply(sumry.hyper.mat, MARGIN = 2, FUN = median),
                                 std  = apply(sumry.hyper.mat, MARGIN = 2, FUN = sd),
                                 quant2.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.025),
                                 quant97.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.975)
                                 )
    rownames(sumry.hyper) = c("kappa (spat. range)", "sigma2 (var. intcpt.)", "tau2 (var.nugget)")

  } else if(data_lik=="NegB"){
    r.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.r[-(1:thin_burn_in)]}) %>% rowMeans
    kappa.samples<- sapply(1:4, function(jj) {max.dist * traceplots.samples[[jj]]$samples.kappa[-(1:thin_burn_in)]}) %>% rowMeans
    sigma2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.sigma2[-(1:thin_burn_in)]}) %>% rowMeans
    tau2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.tau2[-(1:thin_burn_in)]}) %>% rowMeans

    sumry.hyper.mat<- cbind(r.samples, kappa.samples, sigma2.samples, tau2.samples)

    sumry.hyper<- cbind.data.frame(mean = apply(sumry.hyper.mat, MARGIN = 2, FUN = mean),
                                   median = apply(sumry.hyper.mat, MARGIN = 2, FUN = median),
                                   std  = apply(sumry.hyper.mat, MARGIN = 2, FUN = sd),
                                   quant2.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.025),
                                   quant97.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.975)
    )
    rownames(sumry.hyper) = c( "r (shape NegB)", "kappa (spat. range)", "sigma2 (var. intcpt.)", "tau2 (var.nugget)")

  } else{
    k.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.k[-(1:thin_burn_in)]}) %>% rowMeans
    kappa.samples<- sapply(1:4, function(jj) {max.dist * traceplots.samples[[jj]]$samples.kappa[-(1:thin_burn_in)]}) %>% rowMeans
    sigma2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.sigma2[-(1:thin_burn_in)]}) %>% rowMeans
    tau2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.tau2[-(1:thin_burn_in)]}) %>% rowMeans

    sumry.hyper.mat<- cbind(k.samples, kappa.samples, sigma2.samples, tau2.samples)

    sumry.hyper<- cbind.data.frame(mean = apply(sumry.hyper.mat, MARGIN = 2, FUN = mean),
                                   median = apply(sumry.hyper.mat, MARGIN = 2, FUN = median),
                                   std  = apply(sumry.hyper.mat, MARGIN = 2, FUN = sd),
                                   quant2.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.025),
                                   quant97.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.975)
    )
    rownames(sumry.hyper) = c( "k (prec. lNomral)", "kappa (spat. range)", "sigma2 (var. intcpt.)", "tau2 (var.nugget)")

  }
  } else{
    if(data_lik=="Poisson"){
      tau2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.tau2[-(1:thin_burn_in)]}) %>% rowMeans
      sumry.hyper.mat<- matrix(tau2.samples, ncol=1)

      sumry.hyper<- cbind.data.frame(mean = apply(sumry.hyper.mat, MARGIN = 2, FUN = mean),
                                     median = apply(sumry.hyper.mat, MARGIN = 2, FUN = median),
                                     std  = apply(sumry.hyper.mat, MARGIN = 2, FUN = sd),
                                     quant2.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.025),
                                     quant97.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.975)
      )
      rownames(sumry.hyper) = c( "tau2 (var.nugget)")

    } else if(data_lik=="NegB"){
      r.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.r[-(1:thin_burn_in)]}) %>% rowMeans
      tau2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.tau2[-(1:thin_burn_in)]}) %>% rowMeans
      sumry.hyper.mat<- cbind(r.samples, tau2.samples)
      sumry.hyper<- cbind.data.frame(mean = apply(sumry.hyper.mat, MARGIN = 2, FUN = mean),
                                     median = apply(sumry.hyper.mat, MARGIN = 2, FUN = median),
                                     std  = apply(sumry.hyper.mat, MARGIN = 2, FUN = sd),
                                     quant2.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.025),
                                     quant97.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.975)
      )
      rownames(sumry.hyper) = c( "r (shape NegB)", "tau2 (var.nugget)")
    } else{
      k.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.k[-(1:thin_burn_in)]}) %>% rowMeans
      tau2.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.tau2[-(1:thin_burn_in)]}) %>% rowMeans
      sumry.hyper.mat<- cbind(k.samples, tau2.samples)
      sumry.hyper<- cbind.data.frame(mean = apply(sumry.hyper.mat, MARGIN = 2, FUN = mean),
                                     median = apply(sumry.hyper.mat, MARGIN = 2, FUN = median),
                                     std  = apply(sumry.hyper.mat, MARGIN = 2, FUN = sd),
                                     quant2.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.025),
                                     quant97.5  = apply(sumry.hyper.mat, MARGIN = 2, FUN = quantile, probs =0.975)
      )
      rownames(sumry.hyper) = c( "k (prec. lNomral)", "tau2 (var.nugget)")

    }
  }


  ### summary of state variances phi's
 if(model == "dense" | model == "sparse"){
  Wt.samples<-  Reduce(`+`, lapply(1:4, function(jj) traceplots.samples[[jj]]$sample.Wt[-(1:thin_burn_in), ])) / 4
  sumry.var.state<- cbind.data.frame(mean = apply(Wt.samples, MARGIN = 2, FUN = mean),
                                     median = apply(Wt.samples, MARGIN = 2, FUN = median),
                                     std  = apply(Wt.samples, MARGIN = 2, FUN = sd),
                                     quant2.5  = apply(Wt.samples, MARGIN = 2, FUN = quantile, probs =0.025),
                                     quant97.5  = apply(Wt.samples, MARGIN = 2, FUN = quantile, probs =0.975)
  )
  rownames(sumry.var.state) = paste0("phi",1:p)
 } else{
   sumry.var.state<- NULL
 }

  ## summary of covariate coefficients
  beta.samples<-  Reduce(`+`, lapply(1:4, function(jj) traceplots.samples[[jj]]$samples.beta[-(1:thin_burn_in), ])) / 4
  sumry.beta<- cbind.data.frame(mean = apply(beta.samples, MARGIN = 2, FUN = mean),
                                median = apply(beta.samples, MARGIN = 2, FUN = median),
                                std  = apply(beta.samples, MARGIN = 2, FUN = sd),
                                quant2.5  = apply(beta.samples, MARGIN = 2, FUN = quantile, probs =0.025),
                                quant97.5  = apply(beta.samples, MARGIN = 2, FUN = quantile, probs =0.975)
  )
  if(simulation){
  rownames(sumry.beta) = paste0("beta",1:q)
  }else{
  rownames(sumry.beta) = dimnames(X)[[3]]
  }


  if(model == "bayes.reg"){
  beta0.samples<- sapply(1:4, function(jj) {traceplots.samples[[jj]]$samples.beta0[-(1:thin_burn_in)]}) %>% rowMeans
  overall.intcpt<- matrix(beta0.samples, ncol=1)

  sumry.overall.intcpt<- cbind.data.frame(mean = apply(overall.intcpt, MARGIN = 2, FUN = mean),
                                          median = apply(overall.intcpt, MARGIN = 2, FUN = median),
                                          std  = apply(overall.intcpt, MARGIN = 2, FUN = sd),
                                          quant2.5  = apply(overall.intcpt, MARGIN = 2, FUN = quantile, probs =0.025),
                                          quant97.5  = apply(overall.intcpt, MARGIN = 2, FUN = quantile, probs =0.975)
  )
  rownames(sumry.overall.intcpt) = c( "beta0")
  } else{

    sumry.overall.intcpt<- NULL
  }

  ### summary of linear predictors lambda for all the models and all the responses
  lambda.est = Reduce("+", lapply(1:4, function(jj) parallel_chain[[jj]]$within.sample.pred$pred.lambda.mean)) / 4
  lambda.std<- Reduce("+", lapply(1:4, function(jj) parallel_chain[[jj]]$within.sample.pred$pred.lambda.mean)) / 4
  est.info.lambda.st<- list(mean = lambda.est,
                            std= lambda.std,
                            quant2.5 = lambda.est - 1.96 * lambda.std,
                            quant97.5 = lambda.est + 1.96 * lambda.std
                            )

  ### the results of space-time intercepts mu, only for sparse and dense model ####
  if(data_lik == "bayes.reg"){
    est.info.mu.st<- NULL
  } else{
    means.mus.est<-  Reduce("+", lapply(1:4, function(jj){parallel_chain[[jj]]$within.sample.pred$pred.mu.st.mean}))/4
    sd.mus.est<- Reduce("+", lapply(1:4, function(jj){parallel_chain[[jj]]$within.sample.pred$pred.mu.st.sd}))/4

    est.info.mu.st<- list(mean = means.mus.est,
                          std = sd.mus.est,
                          quant2.5 = means.mus.est - 1.96 * sd.mus.est,
                          quant97.5 = means.mus.est + 1.96 * sd.mus.est)
  }

  ### results for space-time spde basis weights parameters R in case of sparse model ####
  if(model =="sparse"){
    means.R.est<-  Reduce("+", lapply(parallel_chain, function(x) x$within.sample.pred$pred.R.mean)) / length(parallel_chain)
    sd.R.est<-  Reduce("+", lapply(parallel_chain, function(x) x$within.sample.pred$pred.R.sd)) / length(parallel_chain)

    est.info.R.st<- list(mean = means.R.est,
                         std = sd.R.est,
                         quant2.5 = means.R.est - 1.96 * sd.R.est,
                         quant97.5 = means.R.est + 1.96 * sd.R.est)

  } else{
    est.info.R.st<- NULL
  }


  #### returning estimates of thetas
  if(model=="bayes.reg"){
    thetas.info<- NULL
    run_time_kappas<- NULL
    run_time_intercepts<- NULL
  } else {
    samples.thetas<- (parallel_chain[[1]]$samples.theta + parallel_chain[[2]]$samples.theta +
                        parallel_chain[[3]]$samples.theta + parallel_chain[[4]]$samples.theta)/4

    est.mean.theta<- apply(samples.thetas, MARGIN = c(1,2), FUN = mean)
    est.sd.theta<- apply(samples.thetas, MARGIN = c(1,2), FUN = sd)
    lci.theta<- apply(samples.thetas, MARGIN = c(1,2), FUN = quantile, probs=0.025)
    uci.theta<- apply(samples.thetas, MARGIN = c(1,2), FUN = quantile, probs=0.975)
    thetas.info<- list(mean = est.mean.theta,
                       std = est.sd.theta,
                       quant2.5 = lci.theta,
                       quant97.5 = uci.theta)
  ###### run times of kappas and updating intercepts
  run_time_kappas<- (parallel_chain[[1]]$runtime$spat_range + parallel_chain[[2]]$runtime$spat_range +
                       parallel_chain[[3]]$runtime$spat_range + parallel_chain[[4]]$runtime$spat_range)/4

  run_time_intercepts<- (parallel_chain[[1]]$runtime$spat_intercept + parallel_chain[[2]]$runtime$spat_intercept +
                           parallel_chain[[3]]$runtime$spat_intercept + parallel_chain[[4]]$runtime$spat_intercept)/4

  }

  ############### space-time predictions #################
  ######## spatial interpolation
  Y.spat.pred.array<- array(dim=c(4, dim(parallel_chain[[1]]$spat.temp.pred$samples.spat.intpol)))
  for (ll in 1:4) {
    Y.spat.pred.array[ll,,,]<- parallel_chain[[ll]]$spat.temp.pred$samples.spat.intpol
  }
  Y.spat.pred.array[Y.spat.pred.array > 1e10 | is.na(Y.spat.pred.array)] <- NA

  # Convert the array to a list to detect `numeric(0)` entries
  Y.list <- as.list(Y.spat.pred.array)
  # Replace `numeric(0)` with NA
  Y.list <- lapply(Y.list, function(x) {
    if (length(x) == 0) {
      return(NA)
    } else {
      return(x)
    }
  })
  Y.spat.pred.array <- array(unlist(Y.list), dim = dim(Y.spat.pred.array))   # Convert the list back to an array with the original dimensions

  Y.spat.pred<- apply(Y.spat.pred.array, MARGIN = c(2,3,4), FUN=function(x) {x[which.min(abs(x - median(x, na.rm=TRUE)))]})
  Y.list <- as.list(Y.spat.pred)
  # Replace `numeric(0)` with NA
  Y.list <- lapply(Y.list, function(x) {
    if (length(x) == 0) {
      return(NA)
    } else {
      return(x)
    }
  })
  Y.spat.pred <- array(unlist(Y.list), dim = dim(Y.spat.pred))   # Convert the list back to an array with the original dimensions

  est.mean.spat<- apply(Y.spat.pred, MARGIN = c(2,3), FUN=function(x) {x[which.min(abs(x - median(x, na.rm=TRUE)))]})
  est.sd.spat<- apply(Y.spat.pred, MARGIN = c(2,3), FUN=sd, na.rm=TRUE)
  lci.spat<- apply(Y.spat.pred, MARGIN = c(2,3), FUN=quantile, probs=0.025, na.rm=TRUE)
  uci.spat<- apply(Y.spat.pred, MARGIN = c(2,3), FUN=quantile, probs=0.975, na.rm=TRUE)

  spat.interpol<- list(true = Y[1:(nt),spatInt.ind],
                       median =  est.mean.spat,
                       std = est.sd.spat,
                       quant2.5 = lci.spat,
                       quant97.5  = uci.spat
                       )

  ############### forecasting #################
  Y.forecast.pred.array<- array(dim=c(4, dim(parallel_chain[[1]]$spat.temp.pred$samples.future.forecast)))
  for (ll in 1:4) {
    Y.forecast.pred.array[ll,,,]<- parallel_chain[[ll]]$spat.temp.pred$samples.future.forecast
  }
  Y.forecast.pred.array[Y.forecast.pred.array > 1e10 | is.na(Y.forecast.pred.array)] <- NA
  # Convert the array to a list to detect `numeric(0)` entries
  Y.list <- as.list(Y.forecast.pred.array)
  # Replace `numeric(0)` with NA
  Y.list <- lapply(Y.list, function(x) {
    if (length(x) == 0) {
      return(NA)
    } else {
      return(x)
    }
  })
  Y.forecast.pred.array <- array(unlist(Y.list), dim = dim(Y.forecast.pred.array))   # Convert the list back to an array with the original dimensions

  Y.forecast.pred<- apply(Y.forecast.pred.array, MARGIN = c(2,3,4), FUN=function(x) {x[which.min(abs(x - median(x, na.rm=TRUE)))]})
  Y.forecast.pred[Y.forecast.pred > 1e10 | is.na(Y.forecast.pred)] <- NA
  # Convert the array to a list to detect `numeric(0)` entries
  Y.list <- as.list(Y.forecast.pred)
  # Replace `numeric(0)` with NA
  Y.list <- lapply(Y.list, function(x) {
    if (length(x) == 0) {
      return(NA)
    } else {
      return(x)
    }
  })
  Y.forecast.pred <- array(unlist(Y.list), dim = dim(Y.forecast.pred))   # Convert the list back to an array with the original dimensions

  est.mean.forecast<- apply(Y.forecast.pred, MARGIN = c(2,3), FUN=function(x) {x[which.min(abs(x - median(x, na.rm=TRUE)))]})
  est.sd.forecast<- apply(Y.forecast.pred, MARGIN = c(2,3), FUN=sd, na.rm=TRUE)
  lci.forecast<- apply(Y.forecast.pred, MARGIN = c(2,3), FUN=quantile, probs=0.025, na.rm=TRUE)
  uci.forecast<- apply(Y.forecast.pred, MARGIN = c(2,3), FUN=quantile, probs=0.975, na.rm=TRUE)

  forecast<- list(true = Y[(nt + forcast.ind), -spatInt.ind],
                  median =  est.mean.forecast,
                  std = est.sd.forecast,
                  quant2.5 = lci.forecast,
                  quant97.5 = uci.forecast
                  )


  ####### spatio-temporal predictions ########
  Y.spat.forecast.pred.array<- array(dim=c(4, dim(parallel_chain[[1]]$spat.temp.pred$samples.spat.future.forecast)))
  for (ll in 1:4) {
    Y.spat.forecast.pred.array[ll,,,]<- parallel_chain[[ll]]$spat.temp.pred$samples.spat.future.forecast
  }

  Y.spat.forecast.pred.array[Y.spat.forecast.pred.array > 1e10 | is.na(Y.spat.forecast.pred.array)] <- NA

  # Convert the array to a list to detect `numeric(0)` entries
  Y.list <- as.list(Y.spat.forecast.pred.array)
  # Replace `numeric(0)` with NA
  Y.list <- lapply(Y.list, function(x) {
    if (length(x) == 0) {
      return(NA)
    } else {
      return(x)
    }
  })
  Y.spat.forecast.pred.array <- array(unlist(Y.list), dim = dim(Y.spat.forecast.pred.array))   # Convert the list back to an array with the original dimensions
  Y.spat.forecast.pred<- apply(Y.spat.forecast.pred.array, MARGIN = c(2,3,4), FUN=function(x) {x[which.min(abs(x - median(x, na.rm=TRUE)))]})

  Y.list <- as.list(Y.spat.forecast.pred)
  # Replace `numeric(0)` with NA
  Y.list <- lapply(Y.list, function(x) {
    if (length(x) == 0) {
      return(NA)
    } else {
      return(x)
    }
  })
  Y.spat.forecast.pred <- array(unlist(Y.list), dim = dim(Y.spat.forecast.pred))   # Convert the list back to an array with the original dimensions


  est.mean.spat.forecast<- apply(Y.spat.forecast.pred, MARGIN = c(2,3), FUN=function(x) {x[which.min(abs(x - median(x, na.rm=TRUE)))]})
  est.sd.spat.forecast<- apply(Y.spat.forecast.pred, MARGIN = c(2,3), FUN=sd,  na.rm=TRUE)
  lci.spat.forecast<- apply(Y.spat.forecast.pred, MARGIN = c(2,3), FUN=quantile, probs=0.025, na.rm=TRUE)
  uci.spat.forecast<- apply(Y.spat.forecast.pred, MARGIN = c(2,3), FUN=quantile, probs=0.975, na.rm=TRUE)

  spat_time_pred<- list(true = Y[(nt + forcast.ind), spatInt.ind],
                        median =  est.mean.spat.forecast,
                        std = est.sd.spat.forecast,
                        quant2.5 = lci.spat.forecast,
                        quant97.5 = uci.spat.forecast
                        )


 # browser()
  ### within sample predictions and missing value imputations using sample medians
  na.ind<-  is.na(Y.o)
  Y.within.sample.mean.array<- array(dim=c(4, length(parallel_chain[[1]]$within.sample.pred$pred.cont.mean)))
  Y.within.sample.sd.array<- array(dim=c(4, length(parallel_chain[[1]]$within.sample.pred$pred.cont.sd)))
  for (ll in 1:4) {
    Y.within.sample.mean.array[ll,]<- parallel_chain[[ll]]$within.sample.pred$pred.cont.mean
    Y.within.sample.sd.array[ll,]<- parallel_chain[[ll]]$within.sample.pred$pred.cont.sd
  }
  est.mean.within.sample<- matrix(NA, nrow = nt, ncol =ns)
  est.mean.within.sample[!na.ind]<- c(apply(Y.within.sample.mean.array, MARGIN = c(2), FUN=median))

  est.sd.within.sample<- matrix(NA, nrow = nt, nco=ns)
  est.sd.within.sample[!na.ind]<- c(apply(Y.within.sample.sd.array, MARGIN = c(2), FUN=median))

  lci.within.sample<- matrix(NA, nrow = nt, ncol=ns)
  if(data_lik=="lNormal"){
    lci.within.sample[!na.ind]<- pmax(-Inf, c(apply(Y.within.sample.mean.array, MARGIN = c(2), FUN=median)) -
                                        1.96 * c(apply(Y.within.sample.sd.array, MARGIN = c(2), FUN=median)))
  } else{
     lci.within.sample[!na.ind]<- pmax(0, c(apply(Y.within.sample.mean.array, MARGIN = c(2), FUN=median)) -
                                         1.96 * c(apply(Y.within.sample.sd.array, MARGIN = c(2), FUN=median)))
  }


  uci.within.sample<- matrix(NA, nrow = nt, ncol=ns)
  uci.within.sample[!na.ind]<- pmin(Inf, c(apply(Y.within.sample.mean.array, MARGIN = c(2), FUN=median)) +
                                      1.96 * c(apply(Y.within.sample.sd.array, MARGIN = c(2), FUN=median)))

  Y.miss.imput.array<- array(dim=c(4, dim(parallel_chain[[1]]$miss.value.imputation.result$samples.miss.imput)))
  for (ll in 1:4) {
    Y.miss.imput.array[ll,,]<- parallel_chain[[ll]]$miss.value.imputation.result$samples.miss.imput
  }
  Y.miss.imput.pred<- apply(Y.miss.imput.array, MARGIN = c(2,3), FUN=function(x) {x[which.min(abs(x - median(x)))]})

  est.mean.within.sample[na.ind]<- apply(Y.miss.imput.pred, MARGIN = c(2), FUN=function(x) {x[which.min(abs(x - median(x)))]})
  est.sd.within.sample[na.ind]<- apply(Y.miss.imput.pred, MARGIN = c(2), FUN=sd)
  lci.within.sample[na.ind]<- apply(Y.miss.imput.pred, MARGIN = c(2), FUN=quantile, probs=0.025)
  uci.within.sample[na.ind]<- apply(Y.miss.imput.pred, MARGIN = c(2), FUN=quantile, probs=0.975)

  within.sample.pred.and.imputations<- list(true = Y.o,
                                            median =  est.mean.within.sample,
                                            std = est.sd.within.sample,
                                            quant2.5 = lci.within.sample,
                                            quant97.5 = uci.within.sample
  )


  ####### Estimations of average Annual daily bike (AADB), site wise
  ### at new spatial locations
  sites.sums.new.loc<- apply(Y.spat.pred, MARGIN = c(1,3), FUN=mean, na.rm=TRUE)
  est.aadb.new.loc<- apply(sites.sums.new.loc, MARGIN = 2,  FUN=function(x) {x[which.min(abs(x - median(x)))]})
  sd.aadb.new.loc<- apply(sites.sums.new.loc, MARGIN = 2,  FUN=sd)
  lci.aadb.new.loc<- apply(sites.sums.new.loc, MARGIN = 2, FUN=quantile, probs=0.025, na.rm=TRUE)
  uci.aadb.new.loc<- apply(sites.sums.new.loc, MARGIN = 2, FUN=quantile, probs=0.975, na.rm=TRUE)
  aadb.new.loc<- list(
                       mean =  est.aadb.new.loc,
                       std = sd.aadb.new.loc,
                       quant2.5 = lci.aadb.new.loc,
                       quant97.5  = uci.aadb.new.loc
  )

 ### at already observed  locations
  aadb.obsd.loc<- list(mean = apply(est.mean.within.sample, MARGIN = 2,  FUN = mean),
                       std = apply(est.sd.within.sample, MARGIN = 2,  FUN = mean),
                       quant2.5 = apply(lci.within.sample, MARGIN = 2,  FUN = mean),
                       quant97.5 = apply(uci.within.sample, MARGIN = 2, FUN = mean)
  )

return(
  list("summery.hyper" = list(sumry.hyper = sumry.hyper, sumry.var.state= sumry.var.state),
       "summary.fixed.effects.coeff" = list(overall.intcp = sumry.overall.intcpt, cov.coeff = sumry.beta),
       "predictions" = list(within.sample = within.sample.pred.and.imputations,
                            spatial = spat.interpol,
                            forecast = forecast,
                            space.time = spat_time_pred),

       "summary.dynamic.temp.coeff" = thetas.info,
       "summary.st.intercepot" = list(spde.basis.weights = est.info.R.st,
                                      mu.st = est.info.mu.st),
       "summary.linear.predictor" =  est.info.lambda.st,
       "traceplots.samples" = traceplots.samples,
       "run.time.costly.comt." = list(kappa = run_time_kappas,
                                     intercept = run_time_intercepts),
       "N.MCMC" =  N.MCMC,
       "original.data" = list(data.fit= Y.o, data.all = Y),
       "forcast.ind" = forcast.ind,
       "spatInt.ind" = spatInt.ind,
       "sparse.info" =  list(INLA.mesh = INLA.mesh,
                             A.proj.o = A.proj.o,
                             A.proj.p =  A.proj.p,
                             c.mat.o =  A.proj.p,
                             g1.mat.o = g1.mat.o,
                             g2.mat.o = g2.mat.o,
                             mesh.hyper = mesh.hyper
                             ),
       "AADB.est" =  list("new.locations" = aadb.new.loc,
                          "obsd.locations" = aadb.obsd.loc
                          ),
       "nt" = nt,
       "ns" = ns,
       "model" = model,
       "data_lik" = data_lik,
       "cor.type" = cor.type,
       "coords.obsd.loc" = loc.o,
       "coords.pred.loc" = loc.p
       ))
}





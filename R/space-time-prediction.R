#' Spatio-Temporal Prediction Based on the Baseline Model
#'
#' This function performs spatio-temporal prediction based on the baseline model for different types of predictions: spatial interpolation at new locations, forecasting at observed sites, or forecasting at new spatial locations.
#'
#' @param nt Integer, total number of temporal points at which we fit the model.
#' @param ns Integer, number of spatial locations where the data is observed and the model is fit.
#' @param beta0.samples A scalar, simulated samples of the overall intercept at a given iteration.
#' @param tau2.samples A scalar, simulated samples of the variance of the nugget.
#' @param pred_type Character, types of prediction. Options are `"spatInt"` for spatial prediction at new spatial locations, `"forecast"` for future forecasting at observed sites, and `"spatInt_forecast"` for future forecasting at new spatial locations.
#' @param spatInt.ind Integer vector, indices of the new locations for spatial interpolation or future forecasting.
#' @param forcast.ind Integer vector, time indices for spatial predictions (e.g., `1` means 1-step ahead prediction; `1:7` means prediction for the next 7 time points).
#' @param data_lik Character, likelihood type. Options include `"Poisson"`, `"NegB"`, and `"lNormal"`.
#' @param X.intpl A matrix of covariates for interpolation.
#' @param X.frcast A matrix of covariates for forecasting.
#' @param X.st A matrix of covariates for spatio-temporal predictions.
#' @param beta.samples A vector of simulated covariate coefficients.
#' @param r.samples A scalar, simulated shape parameter for the negative binomial distribution.
#' @param k.samples A scalar, simulated precision parameter for the log-normal distribution.
#'
#' @return A list containing predicted `y.pred.samples` and `lambda.pred.samples`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' nt <- 100
#' ns <- 50
#' Xs <- matrix(runif((ns + 10) * 5), ncol = 5)
#' Xt <- matrix(runif((nt + 5) * 4), ncol = 4)
#' beta0.samples <- 0.5
#' tau2.samples <- 0.1
#' beta.samples <- runif(5)
#' pred_type <- "forecast"
#' y_pred <- pred.reg(nt, ns, "Poisson", 1:10, 1:5, pred_type, Xs, Xt, Xs, beta0.samples, tau2.samples, beta.samples, NULL, NULL)
#' }
pred.reg<- function(nt,
                    ns,
                    data_lik,
                    spatInt.ind,
                    forcast.ind,
                    pred_type = c("spatInt", "forecast", "spatInt_forecast"),
                    X.intpl,
                    X.frcast,
                    X.st,
                    beta0.samples,
                    tau2.samples,
                    beta.samples,
                    r.samples,
                    k.samples
){
  if(pred_type=="spatInt"){
    ns.pred<- length(spatInt.ind)
    comp.cov<- matrix(c(X.intpl %*% beta.samples), nrow = nt, ncol=ns.pred)

    lambda.samples<- beta0.samples + comp.cov +
      matrix(rnorm(nt * ns.pred, sd=sqrt(tau2.samples)), nrow = nt, ncol = ns.pred)

    if(data_lik=="Poisson"){
    y.pred<- matrix(rpois(n=nt * ns.pred, lambda=exp(lambda.samples)), nrow = nt, ncol = ns.pred)
    } else if(data_lik=="NegB"){
    y.pred<- matrix(rnbinom(n=nt * ns.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt, ncol = ns.pred)
    } else if(data_lik=="lNormal"){
    y.pred<- matrix(exp(rnorm(n=nt * ns.pred, mean = lambda.samples, sd = sqrt(k.samples)))-1, nrow = nt, ncol = ns.pred)
    }

  } else if(pred_type=="forecast"){
    nt.pred<- length(forcast.ind)
    comp.cov<- matrix(c(X.frcast %*% beta.samples), nrow = nt.pred, ncol = ns)

    lambda.samples<- beta0.samples + comp.cov +
      matrix(rnorm(nt.pred * ns, sd=sqrt(tau2.samples)), nrow = nt.pred, ncol = ns)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=ns * nt.pred, lambda=exp(lambda.samples)), nrow = nt.pred, ncol = ns)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=ns * nt.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt.pred, ncol = ns)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=ns * nt.pred, mean = lambda.samples, sd = sqrt(k.samples)))-1, nrow = nt.pred, ncol = ns)
    }


  } else if(pred_type=="spatInt_forecast") {
    nt.pred<- length(forcast.ind)
    ns.pred<- length(spatInt.ind)

    comp.cov<- matrix(c(X.st %*% beta.samples), nrow = nt.pred, ncol = ns.pred)
    lambda.samples<- beta0.samples + comp.cov +
      matrix(rnorm(nt.pred * ns.pred, sd=sqrt(tau2.samples)), nrow = nt.pred, ncol = ns.pred)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=ns.pred * nt.pred, lambda=exp(lambda.samples)), nrow = nt.pred, ncol = ns.pred)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=ns.pred * nt.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt.pred, ncol = ns.pred)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=ns.pred * nt.pred, mean = lambda.samples, sd = sqrt(k.samples)))-1, nrow = nt.pred, ncol = ns.pred)
    }
  }

  return(list(y.pred.samples =  y.pred, lambda.pred.samples=lambda.samples))
}

#' Spatio-Temporal Prediction Based on the Fitted DGLM Full Model with Dense Correlation Structures
#'
#' This function performs spatio-temporal prediction based on the fitted DGLM full model with dense correlation structures. It supports three types of predictions: spatial interpolation at new locations, forecasting at observed sites, and forecasting at new spatial locations.
#'
#' @param nt Integer, total number of temporal points at which the model is fitted.
#' @param ns Integer, number of spatial locations where the data is observed and the model is fitted.
#' @param tau2.samples A scalar, posterior samples of the variance of the nugget term.
#' @param beta.samples A vector, posterior samples of covariate coefficients corresponding to design matrix `X`.
#' @param pred_type Character, type of prediction. Options are `"spatInt"` for spatial prediction at new spatial locations, `"forecast"` for future forecasting at observed sites, and `"spatInt_forecast"` for future forecasting at new spatial locations.
#' @param spatInt.ind Integer vector, indices of new locations for spatial interpolation or future forecasting.
#' @param forcast.ind Integer vector, time indices for spatial predictions (e.g., `1` means 1-step ahead prediction; `1:7` means prediction for the next 7 time points).
#' @param G A matrix of dimension `p x p`, state-evolving matrix.
#' @param kappa.samples A scalar, posterior samples of the range parameter.
#' @param sigma2.samples A scalar, posterior samples of the variance of the spatial covariance function.
#' @param theta.samples A matrix of dimension `nt x p`, posterior samples of state vectors.
#' @param Wt.samples A vector of length `p`, posterior samples of the variance of state vectors.
#' @param mus.samples A matrix of dimension `nt x ns`, posterior samples of spatially varying intercept parameters.
#' @param data_lik Character, likelihood type. Options include `"Poisson"`, `"NegB"`, and `"lNormal"`.
#' @param X.intpl A matrix of covariates for interpolation.
#' @param X.frcast A matrix of covariates for forecasting.
#' @param X.st A matrix of covariates for spatio-temporal predictions.
#' @param r.samples A scalar, posterior samples of the shape parameter for the negative binomial distribution.
#' @param k.samples A scalar, posterior samples of the precision parameter for the log-normal distribution.
#' @param Ft.o A matrix of temporal covariates at observed time points.
#' @param Ft.frcast A matrix of temporal covariates for future forecasting.
#' @param dist.mat.all A distance matrix of the spatial locations.
#' @param cor.type Character, specifying the correlation type for the dense model. Options include `"Matern0.5"`, `"Matern1"`, `"Matern1.5"`, `"Matern2.5"`, and `"MaternInf"`.
#'
#' @return A list containing predicted `y.pred` and `lambda.samples`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' nt <- 100
#' ns <- 50
#' beta.samples <- runif(5)
#' kappa.samples <- 0.5
#' sigma2.samples <- 0.2
#' tau2.samples <- 0.1
#' theta.samples <- matrix(runif(nt * 3), ncol = 3)
#' pred.dense(nt, ns, "Poisson", 1:10, 1:5, "forecast", X.intpl, X.frcast, X.st, Ft, Ft, G, dist.mat, kappa.samples, sigma2.samples, tau2.samples, theta.samples, Wt.samples, beta.samples, mus.samples)
#' }
pred.dense <- function(nt,
                       ns,
                       data_lik,
                       spatInt.ind,
                       forcast.ind,
                       cor.type,
                       pred_type = c("spatInt", "forecast", "spatInt_forecast"),
                       X.intpl,
                       X.frcast,
                       X.st,
                       Ft.o,
                       Ft.frcast,
                       G,
                       dist.mat.all,
                       kappa.samples,
                       sigma2.samples,
                       tau2.samples,
                       theta.samples,
                       Wt.samples,
                       beta.samples,
                       mus.samples,
                       r.samples,
                       k.samples) {


  if(pred_type=="spatInt"){
    ns.pred <- length(spatInt.ind)
    Sigma.kappa <- sigma2.samples * cor.fun(cor.type = cor.type, dist.mat = dist.mat.all, kappa = kappa.samples) # sigma2.samples * exp(-dist.mat.all / kappa.samples)
    Sigma.pred <- Sigma.kappa[spatInt.ind, spatInt.ind]
    Sigma.obs <- Sigma.kappa[-spatInt.ind, -spatInt.ind]
    inv_Sigma.obs <- solve(Sigma.obs)
    Sigma.pred.obs <- Sigma.kappa[spatInt.ind, -spatInt.ind]
    cov_sim <- Sigma.pred - Sigma.pred.obs %*% inv_Sigma.obs %*% t(Sigma.pred.obs)
    mus.mean <- t(Sigma.pred.obs %*% inv_Sigma.obs %*% t(mus.samples))

    mus.pred <- mus.mean #c(mvtnorm::rmvnorm(n = 1, mean = mus.mean, sigma =  as.matrix(Matrix::forceSymmetric(cov_sim))))
    Ft.mean <- matrix(rowSums(Ft.o[1:nt,] * theta.samples), nrow = nt, ncol = ns.pred)
    comp.cov<- matrix(c(X.intpl %*% beta.samples), nrow = nt, ncol=ns.pred)

    lambda.samples<- mus.pred + comp.cov + Ft.mean +
      matrix(rnorm(nt * ns.pred, sd=sqrt(tau2.samples)), nrow = nt, ncol = ns.pred)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=nt * ns.pred, lambda=exp(lambda.samples)), nrow = nt, ncol = ns.pred)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=nt * ns.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt, ncol = ns.pred)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=nt * ns.pred, mean = lambda.samples, sd = sqrt(k.samples)))-1, nrow = nt, ncol = ns.pred)
    }

  } else if(pred_type=="forecast"){
    mu.nt<- mus.samples[nt,]
    theta0<- theta.samples[nt,] ### the last time of the sequences would be the initial time for the forecasting
    p<- length(theta0)
    nt.pred<- length(forcast.ind)
    theta.pred<- matrix(NA, nrow=nt.pred, ncol=p)
    theta.pred[1,]<- c(G %*% theta0) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    for (t in 2:nt.pred) {
      theta.pred[t,]<- c(G %*% theta.pred[t-1,]) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    }

    ns.pred <- length(spatInt.ind)
    Sigma.kappa <- sigma2.samples * cor.fun(cor.type = cor.type, dist.mat = dist.mat.all, kappa = kappa.samples) # sigma2.samples * exp(-dist.mat.all / kappa.samples)
    Sigma.obs <- Sigma.kappa[-spatInt.ind, -spatInt.ind]
    mus.pred<-  matrix(NA, nrow=nt.pred, ncol=ns)
    mus.pred[1,]<- mu.nt + rnorm(n=ns, sd = sqrt(sigma2.samples))
    for (t in 2:nt.pred) {
      mus.pred[t,]<-  mus.pred[t-1,] + rnorm(n=ns, sd = sqrt(sigma2.samples))
    }
    Ft.mean<- matrix(rowSums(Ft.frcast * theta.pred), nrow = nt.pred, ncol = ns)
    comp.cov<- matrix(c(X.frcast %*% beta.samples), nrow = nt.pred, ncol = ns)

    lambda.samples<-  mus.pred + comp.cov + Ft.mean +
      matrix(rnorm(nt.pred * ns, sd=sqrt(tau2.samples)), nrow = nt.pred, ncol = ns)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=ns * nt.pred, lambda=exp(lambda.samples)), nrow = nt.pred, ncol = ns)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=ns * nt.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt.pred, ncol = ns)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=ns * nt.pred, mean = lambda.samples, sd = sqrt(k.samples))) -1, nrow = nt.pred, ncol = ns)
    }

  } else if(pred_type=="spatInt_forecast"){
    ns.pred <- length(spatInt.ind)
    nt.pred<- length(forcast.ind)
    Sigma.kappa <- sigma2.samples * cor.fun(cor.type = cor.type, dist.mat = dist.mat.all, kappa = kappa.samples) # sigma2.samples * exp(-dist.mat.all / kappa.samples)
    Sigma.pred <- Sigma.kappa[spatInt.ind, spatInt.ind]
    Sigma.obs <- Sigma.kappa[-spatInt.ind, -spatInt.ind]
    inv_Sigma.obs <- solve(Sigma.obs)
    Sigma.pred.obs <- Sigma.kappa[spatInt.ind, -spatInt.ind]
    cov_sim <- Sigma.pred - Sigma.pred.obs %*% inv_Sigma.obs %*% t(Sigma.pred.obs)
    mus.mean.nt <-  c(Sigma.pred.obs %*% inv_Sigma.obs %*% (mus.samples[nt, ]))  # c(Sigma.pred.obs %*% inv_Sigma.obs %*% mus.samples)
    mus.pred.nt <- mus.mean.nt # c(mvtnorm::rmvnorm(n = 1, mean = mus.mean.nt, sigma = as.matrix(Matrix::forceSymmetric(cov_sim))))

    mus.pred<-  matrix(NA, nrow=nt.pred, ncol=ns.pred)
    mus.pred[1,]<- mus.pred.nt + mvtnorm::rmvnorm(n=1, sigma = as.matrix(Matrix::forceSymmetric(cov_sim)))
    for (t in 2:nt.pred) {
      mus.pred[t,]<-  mus.pred[t-1,] + mvtnorm::rmvnorm(n=1, sigma = as.matrix(Matrix::forceSymmetric(cov_sim)))
    }

    comp.cov<- matrix(c(X.st %*% beta.samples), nrow = nt.pred, ncol = ns.pred)

    theta0<- theta.samples[nt,]
    p<- length(theta0)
    theta.pred<- matrix(NA, nrow=nt.pred, ncol=p)
    theta.pred[1,]<- c(G %*% theta0) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    for (t in 2:nt.pred) {
      theta.pred[t,]<- c(G %*% theta.pred[t-1,]) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    }

    Ft.mean<- matrix(rowSums(Ft.frcast * theta.pred), nrow = nt.pred, ncol = ns.pred)

    lambda.samples<- mus.pred + Ft.mean + comp.cov +
      matrix(rnorm(nt.pred * ns.pred, sd=sqrt(tau2.samples)), nrow = nt.pred, ncol = ns.pred)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=ns.pred * nt.pred, lambda=exp(lambda.samples)), nrow = nt.pred, ncol = ns.pred)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=ns.pred * nt.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt.pred, ncol = ns.pred)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=ns.pred * nt.pred, mean = lambda.samples, sd = sqrt(k.samples))) -1, nrow = nt.pred, ncol = ns.pred)
    }

    }
  return(list(y.pred = y.pred, lambda.samples=lambda.samples))
}


#' Spatio-Temporal Prediction Based on the Fitted Sparse SPDE DGLM Model
#'
#' This function performs spatio-temporal predictions using the fitted sparse SPDE DGLM model. It supports three types of predictions: spatial interpolation at new locations, forecasting at observed sites, and forecasting at new spatial locations.
#'
#' @param nt Integer, total number of temporal points at which the model is fitted.
#' @param ns Integer, number of spatial locations where the data is observed and the model is fitted.
#' @param tau2.samples Numeric, posterior samples of the variance of the nugget term.
#' @param beta.samples Numeric vector, posterior samples of the covariate coefficients corresponding to design matrix `X`.
#' @param pred_type Character, type of prediction. Options are `"spatInt"` for spatial prediction at new spatial locations, `"forecast"` for future forecasting at observed sites, and `"spatInt_forecast"` for future forecasting at new spatial locations.
#' @param spatInt.ind Integer vector, indices of new locations for spatial interpolation or future forecasting.
#' @param forcast.ind Integer vector, time indices for spatial predictions (e.g., `1` means 1-step ahead prediction; `1:7` means prediction for the next 7 time points).
#' @param G Matrix of dimension `p x p`, state-evolving matrix.
#' @param kappa.samples Numeric, posterior samples of the range parameter.
#' @param sigma2.samples Numeric, posterior samples of the variance of the spatial covariance function.
#' @param theta.samples Matrix of dimension `nt x p`, posterior samples of state vectors.
#' @param Wt.samples Numeric vector of length `p`, posterior samples of the variance of state vectors.
#' @param Rs.samples Numeric vector of length `nb`, posterior samples of SPDE basis weights.
#' @param nb Integer, number of mesh nodes.
#' @param X.intpl Numeric matrix, covariate values for spatial interpolation.
#' @param X.frcast Numeric matrix, covariate values for future time points.
#' @param X.st Numeric matrix, covariate values for both space and time.
#' @param Ft.o Numeric matrix, temporal covariates at observed times.
#' @param Ft.frcast Numeric matrix, temporal covariates for future forecasting.
#' @param A.proj.o Numeric matrix of dimension `ns x nb`, projection matrix for observed locations.
#' @param A.proj.p Numeric matrix of dimension `length(spatInt.ind) x nb`, projection matrix for new spatial locations.
#' @param c.mat Numeric matrix, INLA mesh components.
#' @param g1.mat Numeric matrix, INLA mesh components.
#' @param g2.mat Numeric matrix, INLA mesh components.
#' @param r.samples Numeric, posterior samples of the shape parameter for the negative binomial distribution.
#' @param k.samples Numeric, posterior samples of the precision parameter for the log-normal distribution.
#'
#' @return A list containing predicted `y.pred` and `lambda.samples`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' pred.sparse(nt = 100, ns = 50, nb = 30, spatInt.ind = 1:10, forcast.ind = 1:5,
#'             pred_type = "forecast", X.intpl = X.intpl, X.frcast = X.frcast,
#'             X.st = X.st, Ft.o = Ft.o, Ft.frcast = Ft.frcast, G = G, A.proj.o = A.proj.o,
#'             A.proj.p = A.proj.p, c.mat = c.mat, g1.mat = g1.mat, g2.mat = g2.mat,
#'             kappa.samples = kappa.samples, sigma2.samples = sigma2.samples,
#'             tau2.samples = tau2.samples, theta.samples = theta.samples,
#'             Wt.samples = Wt.samples, beta.samples = beta.samples, Rs.samples = Rs.samples)
#' }
pred.sparse <- function(nt,
                        ns,
                        nb,
                        data_lik,
                        spatInt.ind,
                        forcast.ind,
                        pred_type = c("spatInt", "forecast", "spatInt_forecast"),
                        X.intpl,
                        X.frcast,
                        X.st,
                        Ft.o,
                        Ft.frcast,
                        G,
                        A.proj.o,
                        A.proj.p,
                        c.mat,
                        g1.mat,
                        g2.mat,
                        kappa.samples,
                        sigma2.samples,
                        tau2.samples,
                        theta.samples,
                        Wt.samples,
                        beta.samples,
                        Rs.samples,
                        r.samples,
                        k.samples
                        ) {

 # browser()
  if(pred_type=="spatInt"){
    ns.pred <- length(spatInt.ind)

    mus.pred <- t(A.proj.p %*% t(Rs.samples))
    Ft.mean <- matrix(rowSums(Ft.o * theta.samples), nrow = nt, ncol = ns.pred)
    comp.cov<- matrix(c(X.intpl %*% beta.samples), nrow = nt, ncol=ns.pred)

    lambda.samples<- mus.pred + Ft.mean + comp.cov +
      matrix(rnorm(nt * ns.pred, sd=sqrt(tau2.samples)), nrow = nt, ncol = ns.pred)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=nt * ns.pred, lambda=exp(lambda.samples)), nrow = nt, ncol = ns.pred)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=nt * ns.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt, ncol = ns.pred)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=nt * ns.pred, mean = lambda.samples, sd = sqrt(k.samples))) -1, nrow = nt, ncol = ns.pred)
    }

  } else if(pred_type=="forecast"){
    theta0<- theta.samples[nt,] ### the last time of the sequences would be the initial time for the forecasting
    p<- length(theta0)
    nt.pred<- length(forcast.ind)

    ### temporal components of harmonics and other components
    theta.pred<- matrix(NA, nrow=nt.pred, ncol=p)
    theta.pred[1,]<- c(G %*% theta0) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    for (t in 2:nt.pred) {
      theta.pred[t,]<- c(G %*% theta.pred[t-1,]) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    }
    Ft.mean<- matrix(rowSums(Ft.frcast * theta.pred), nrow = nt.pred, ncol = ns)

    ### temporal prediction of space-time varying components
    Q.rho.det <- covmat.inv.update.inla(rho = kappa.samples,  mar.std = sqrt(sigma2.samples), c.mat = c.mat,
                                        g1.mat = g1.mat, g2.mat = g2.mat, alpha = 2, logDet = FALSE)
    Q.Z0 <- as.matrix(Q.rho.det$cormat.inv)

    R.st<- matrix(NA, nrow=nt.pred, ncol=nb)
    mus<- matrix(NA, nrow=nt.pred, ncol=ns)

    R.st[1,] = Rs.samples[nt,]+ rnorm(n=nb, sd=sqrt(sigma2.samples))
    mus[1,]<- as.numeric(A.proj.o %*% R.st[1,])
    for (t in 2:nt.pred) {
      R.st[t,]<-  R.st[t-1,] + rnorm(n=nb, sd=sqrt(sigma2.samples))
      mus[t,]<- as.numeric(A.proj.o %*% R.st[t,])
    }
    comp.cov<- matrix(c(X.frcast %*% beta.samples), nrow = nt.pred, ncol = ns)

    lambda.samples<- mus + comp.cov + Ft.mean +
      matrix(rnorm(nt.pred * ns, sd=sqrt(tau2.samples)), nrow = nt.pred, ncol = ns)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=ns * nt.pred, lambda=exp(lambda.samples)), nrow = nt.pred, ncol = ns)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=ns * nt.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt.pred, ncol = ns)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=ns * nt.pred, mean = lambda.samples, sd = sqrt(k.samples))) -1, nrow = nt.pred, ncol = ns)
    }

  } else if(pred_type=="spatInt_forecast"){

    ### temporal prediction of space-time varying components
    ns.pred <- length(spatInt.ind)
    nt.pred<- length(forcast.ind)

    Q.rho.det <- covmat.inv.update.inla(rho = kappa.samples,  mar.std = sqrt(sigma2.samples), c.mat = c.mat,
                                        g1.mat = g1.mat, g2.mat = g2.mat, alpha = 2, logDet = FALSE)
    Q.Z0 <- as.matrix(Q.rho.det$cormat.inv)

    R.st<- matrix(NA, nrow=nt.pred, ncol=nb)
    mus.pred<- matrix(NA, nrow=nt.pred, ncol=ns.pred)

    R.st[1,] = Rs.samples[nt,] + c(rMVNormP(n = 1, mu = rep(0, nb), Omega = Q.Z0))
    mus.pred[1,]<- as.numeric(A.proj.p %*% R.st[1,])
    for (t in 2:nt.pred) {
      R.st[t,]<-  R.st[t-1,] + c(rMVNormP(n = 1, mu = rep(0, nb), Omega = Q.Z0))
      mus.pred[t,]<- as.numeric(A.proj.p %*% R.st[t,])
    }

    theta0<- theta.samples[nt,]
    p<- length(theta0)

    theta.pred<- matrix(NA, nrow=nt.pred, ncol=p)
    theta.pred[1,]<- c(G %*% theta0) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    for (t in 2:nt.pred) {
      theta.pred[t,]<- c(G %*% theta.pred[t-1,]) + c(mvtnorm::rmvnorm(n=1, sigma = diag(Wt.samples, p)))
    }

    Ft.mean<- matrix(rowSums(Ft.frcast * theta.pred), nrow = nt.pred, ncol = ns.pred)
    comp.cov<- matrix(c(X.st %*% beta.samples), nrow = nt.pred, ncol = ns.pred)

    lambda.samples<-  mus.pred + comp.cov + Ft.mean +
      matrix(rnorm(nt.pred * ns.pred, sd=sqrt(tau2.samples)), nrow = nt.pred, ncol = ns.pred)

    if(data_lik=="Poisson"){
      y.pred<- matrix(rpois(n=ns.pred * nt.pred, lambda=exp(lambda.samples)), nrow = nt.pred, ncol = ns.pred)
    } else if(data_lik=="NegB"){
      y.pred<- matrix(rnbinom(n=ns.pred * nt.pred, size  = r.samples, mu=exp(lambda.samples)), nrow = nt.pred, ncol = ns.pred)
    } else if(data_lik=="lNormal"){
      y.pred<- matrix(exp(rnorm(n=ns.pred * nt.pred, mean = lambda.samples, sd = sqrt(k.samples))) -1, nrow = nt.pred, ncol = ns.pred)
    }
  }

  return(list(y.pred = y.pred, lambda.samples=lambda.samples))
}

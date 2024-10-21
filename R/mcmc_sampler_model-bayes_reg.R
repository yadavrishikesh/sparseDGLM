#' MCMC Sampler for Bayesian Regression Model With Only Fixed Effects
#'
#' This function performs MCMC sampling for a Bayesian regression model with dynamic temporal covariates. The sampler supports various likelihood functions and handles missing data, spatial interpolation, and forecasting.
#'
#' @param data_lik Character, specifying the data likelihood, e.g., `"Poisson"`, `"NegB"`, etc.
#' @param nt Integer, the number of time points in the dataset.
#' @param ns Integer, the number of spatial locations.
#' @param p Integer, the number of dynamic temporal covariates.
#' @param q Integer, the number of fixed covariates.
#' @param Y.o Matrix of observed data (excluding spatial interpolation indices).
#' @param data.log.mean Numeric, the mean of the log-transformed observed data.
#' @param data.mean Numeric, the mean of the observed data.
#' @param Y.intpl Matrix of data for spatial interpolation.
#' @param Y.frcast Matrix of forecasted data for future time points.
#' @param Y.st Matrix of spatio-temporal data.
#' @param ind_NA_Y Logical matrix indicating missing values in `Y.o`.
#' @param X.o Matrix of fixed covariates for observed data.
#' @param quad.X.o Quadratic form of the design matrix `X.o`.
#' @param X.intpl Matrix of fixed covariates for spatial interpolation.
#' @param X.frcast Matrix of fixed covariates for forecasting.
#' @param X.st Matrix of fixed covariates for spatio-temporal data.
#' @param X.p Matrix of fixed covariates for prediction.
#' @param spatInt.ind Vector of spatial interpolation indices.
#' @param forcast.ind Vector of forecast indices.
#' @param samples.store Integer, the number of MCMC samples to store.
#' @param N.MCMC Integer, the total number of MCMC iterations.
#' @param burn_in1 Integer, the number of iterations for the first burn-in period.
#' @param burn_in2 Integer, the number of iterations for the second burn-in period.
#' @param adapt_seq Vector, indicating the sequence of iterations for adaptation.
#' @param thin Integer, the thinning interval for MCMC samples.
#' @param adapt Integer, the adaptation interval.
#' @param tun_r Numeric, tuning parameter for `r` in the negative binomial likelihood.
#' @param tun_lambda Numeric, tuning parameter for `lambda`.
#' @param print.result Logical, if `TRUE`, prints the MCMC results.
#' @param traceplot Logical, if `TRUE`, generates trace plots of MCMC samples.
#' @param true.values List of true parameter values (optional).
#' @param simulation Logical, if `TRUE`, the function runs in simulation mode.
#' @param init.seed Integer, seed for random number generation to ensure reproducibility.
#'
#' @return A list containing MCMC samples, trace plots, and convergence diagnostics.
#'
#' @export
#'
#' @examples
#' data_lik <- "Poisson"
#' nt <- 10
#' ns <- 5
#' p <- 3
#' q <- 2
#' Y.o <- matrix(rpois(50, lambda = 3), nrow = 10)
#' data.log.mean <- log(mean(Y.o))
#' data.mean <- mean(Y.o)
#' Y.intpl <- matrix(rpois(20, lambda = 3), nrow = 4)
#' Y.frcast <- matrix(rpois(10, lambda = 3), nrow = 2)
#' Y.st <- matrix(rpois(20, lambda = 3), nrow = 4)
#' X.o <- matrix(rnorm(100), nrow = 10)
#' quad.X.o <- t(X.o) %*% X.o
#' X.intpl <- matrix(rnorm(40), nrow = 4)
#' X.frcast <- matrix(rnorm(20), nrow = 2)
#' X.st <- matrix(rnorm(40), nrow = 4)
#' spatInt.ind <- 1:2
#' forcast.ind <- 1:2
#' samples.store <- 100
#' N.MCMC <- 1000
#' burn_in1 <- 200
#' burn_in2 <- 500
#' adapt_seq <- seq(100, 1000, by = 100)
#' thin <- 10
#' adapt <- 100
#' tun_r <- 0.1
#' tun_lambda <- 0.1
#' print.result <- TRUE
#' traceplot <- TRUE
#' true.values <- NULL
#' simulation <- FALSE
#' init.seed <- 123
#' result <- MCMC.sampler_model.reg(data_lik, nt, ns, p, q, Y.o, data.log.mean, data.mean,
#'                                  Y.intpl, Y.frcast, Y.st, NULL, X.o, quad.X.o, X.intpl,
#'                                  X.frcast, X.st, X.o, spatInt.ind, forcast.ind, samples.store,
#'                                  N.MCMC, burn_in1, burn_in2, adapt_seq, thin, adapt, tun_r,
#'                                  tun_lambda, print.result, traceplot, true.values, simulation, init.seed)
MCMC.sampler_model.reg<-function(model,
                                data_lik,
                                nt,
                                ns,
                                p,
                                q,
                                Y.o,
                                data.log.mean,
                                data.mean,
                                Y.intpl,
                                Y.frcast,
                                Y.st,
                                ind_NA_Y,
                                X.o,
                                quad.X.o,
                                X.intpl,
                                X.frcast,
                                X.st,
                                X.p,
                                spatInt.ind,
                                forcast.ind,
                                samples.store,
                                N.MCMC,
                                burn_in1,
                                burn_in2,
                                adapt_seq,
                                thin,
                                adapt,
                                tun_r,
                                tun_lambda,
                                print.result,
                                traceplot,
                                true.values,
                                simulation,
                                init.seed)
{
  tun_lambda<- matrix(rep(tun_lambda, nt * ns), nrow = nt, ncol = ns)
  if(data_lik=="NegB"){tun_r= tun_r} else {tun_r=NULL}
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt), ncol= 2 + length(tun_r), byrow=TRUE)
  sigma.matrix[1,]<-c(tun_lambda[1], tun_lambda[2], tun_r)

  init.names<- init_fun(p=p, q=q, ns = ns, nt = nt, nb = NULL, Y=Y.o, model=model,
                        data_lik = data_lik, delta=delta, seed=init.seed)

  cur.samples.r<- init.names$init$r
  cur.samples.k<- init.names$init$k
  cur.samples.beta0<- init.names$init$beta0
  cur.samples.tau2<- init.names$init$tau2
  cur.samples.beta<- init.names$init$beta
  cur.samples.lambda<- init.names$init$lambda_init

  samples <- matrix(nrow=floor(N.MCMC/thin), ncol= length(cur.samples.r) + length(cur.samples.k) +
                      2 + q + 4, byrow=TRUE)
  samples[1,]<- c(cur.samples.r, cur.samples.k, cur.samples.beta0,
                  cur.samples.tau2, cur.samples.beta,
                  cur.samples.lambda[c(1, 3, nt+1, nt+3)]) ## saving only few samples

  #### saving the samples of the few parameters for the prediction later
  no.samples<- samples.store
  samples.save<-  floor(seq(from=burn_in1+burn_in2+1, to=N.MCMC, length.out=no.samples))
  samples.save[no.samples]<- samples.save[no.samples]-1 ### last samples
  samples.miss.imput.Y.save<- array(NA, dim = c(no.samples, sum(ind_NA_Y)))
  samples.miss.imput.lambda.save<- array(NA, dim = c(no.samples, sum(ind_NA_Y)))
  samples.spat.pred.save<- array(NA, dim = c(no.samples, nt, length(spatInt.ind)))
  samples.forecast.save<- array(NA, dim = c(no.samples, length(forcast.ind), ns))
  samples.spatInt_forecast.save<- array(NA, dim = c(no.samples, length(forcast.ind), length(spatInt.ind)))

  ### assessing within sample predictions
  sum.post.Y.save<- rep(0, sum(!ind_NA_Y))
  sum.post.Y2.save<- rep(0, sum(!ind_NA_Y))
  sum.post.lambda.save<- matrix(0, nrow = nt, ncol=ns)
  sum.post.lambda2.save<- matrix(0, nrow = nt, ncol = ns)

  ls<- 1

  comp.cov<-  matrix(c(X.o %*% cur.samples.beta), nrow = nt, ncol = ns)

  j<-1
  l<-1
  m<-1
  k<-1
  for (i in 1:(N.MCMC-1)){
    if(((i%%(adapt))-1==0) & (i< (burn_in1+burn_in2+2))){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
      rate.lambda<- matrix(0, nrow = nt, ncol = ns)
      rate.r<- 0

    }

    #browser()

    ######## imputations #########
      if(data_lik=="Poisson"){
        y.imputed<- matrix(rpois(n  = ns * nt, lambda = exp(cur.samples.lambda)), nrow = nt, ncol=ns)
      } else if(data_lik=="NegB"){
        y.imputed<- matrix(rnbinom(n  = ns * nt, size = cur.samples.r, mu = exp(cur.samples.lambda)), nrow = nt, ncol=ns)
      } else if(data_lik =="lNormal"){
        y.imputed<-  matrix(exp(rnorm(n = ns * nt, mean = cur.samples.lambda, sd = sqrt(cur.samples.k))) - 1, nrow = nt, ncol=ns)
    }
    Y.o[ind_NA_Y]<- y.imputed[ind_NA_Y]



    ### update r
    if(data_lik == "NegB"){
      r.all<- update_r(Y = Y.o,
                       lambda= cur.samples.lambda,
                       r = cur.samples.r,
                       tun.r = tun_r)
      rate.r<- ifelse(r.all$ind, rate.r + 1, rate.r)
      tun_r<- adpative_function(index_MCMC_iter=i, sigma2_adapt=tun_r, target_accept=0.40,
                                    rate_adapt=rate.r, burn_in1=burn_in1, burn_in2=burn_in2, adapt_seq=adapt_seq,
                                    adapt=adapt, adpat_param=1, lower.acc=0.35, upper.acc=0.50)
      cur.samples.r<- r.all$r
    } else{
      cur.samples.r <- NULL
    }

    ### update k
    if(data_lik=="lNormal"){
    cur.samples.k<- update_var(ns=ns,
                               nt=nt,
                               param.hyperprior=c(0.25, 0.25),
                               qud_sum = sum((log(1+Y.o) - cur.samples.lambda)^2)
    )
    } else{
      cur.samples.k<- NULL
    }

    #cur.samples.k<- hyper$k

    #### update tau2 #############
    cur.samples.tau2<- update_var(ns=ns,
                                  nt=nt,
                                  param.hyperprior=c(0.25, 0.25),
                                  qud_sum = sum((cur.samples.lambda - cur.samples.beta0 - comp.cov)^2)
    )
    #cur.samples.tau2<- hyper$tau2

    ## update overall intercepts
    cur.samples.beta0<- update_intercept_reg(nt = nt,
                                        ns = ns,
                                        lambda.mean.diff =  cur.samples.lambda - comp.cov,
                                        tau2 = cur.samples.tau2,
                                        prior.var = 3
                                        )



     ### update betas
     cur.samples.beta<- update_space_time_beta(lambda = cur.samples.lambda,
                                             m = cur.samples.beta0,
                                             X = X.o,
                                             quad.X = quad.X.o,
                                             tau2 = cur.samples.tau2
                                             )

    comp.cov<-  matrix(c(X.o%*% cur.samples.beta), nrow = nt, ncol = ns)

    ########## update lambda
    cur.samples.lambda.all<- update_lambda_MMALA(ns = ns,
                                                 nt=nt,
                                                 Y = Y.o,
                                                 cur.lambda = cur.samples.lambda,
                                                 mean.lambda = cur.samples.beta0 + comp.cov,
                                                 tau2 = cur.samples.tau2,
                                                 r=  cur.samples.r,
                                                 k = cur.samples.k,
                                                 data_lik = data_lik,
                                                 tun_lambda = tun_lambda)
    rate.lambda<- ifelse(cur.samples.lambda.all$ind.acc, rate.lambda+1, rate.lambda)
    tun_lambda<- adpative_function(index_MCMC_iter=i, sigma2_adapt=tun_lambda, target_accept=0.60,
                                   rate_adapt=rate.lambda, burn_in1=burn_in1, burn_in2=burn_in2,
                                   adapt=adapt, adpat_param=1, adapt_seq=adapt_seq, lower.acc=0.50, upper.acc=0.70)

    cur.samples.lambda<- cur.samples.lambda.all$cur_lambda

    #cur.samples.lambda<- sim_M0$lambda[(1:nt), -spatInt.ind]

    #browser()

    ## storing within samples prediction performances
    if(i > burn_in1 + burn_in2){
      sum.post.Y.save<- sum.post.Y.save + (y.imputed[!ind_NA_Y])/data.mean
      sum.post.Y2.save<- sum.post.Y2.save + ((y.imputed[!ind_NA_Y])^2)/(data.mean^2)
      sum.post.lambda.save<- sum.post.lambda.save + cur.samples.lambda/data.log.mean
      sum.post.lambda2.save<- sum.post.lambda2.save + ((cur.samples.lambda)^2)/ (data.log.mean^2)

    }

    if(i == samples.save[ls]){
      samples.miss.imput.Y.save[ls,]<- y.imputed[ind_NA_Y]
      samples.spat.pred.save[ls,,]<- pred.reg(nt = nt, ns = ns,
                                              data_lik =  data_lik,
                                              spatInt.ind = spatInt.ind,
                                              forcast.ind = forcast.ind,
                                              pred_type = "spatInt",
                                              X.intpl = X.intpl,
                                              X.frcast = X.frcast,
                                              X.st = X.st,
                                              beta0.samples = cur.samples.beta0,
                                              tau2.samples = cur.samples.tau2,
                                              beta.samples = cur.samples.beta,
                                              r.samples = cur.samples.r,
                                              k.samples = cur.samples.k)$y.pred

      samples.forecast.save[ls,,]<- pred.reg(nt = nt, ns = ns,
                                             data_lik =  data_lik,
                                             spatInt.ind = spatInt.ind,
                                             forcast.ind = forcast.ind,
                                             pred_type = "forecast",
                                             X.intpl = X.intpl,
                                             X.frcast = X.frcast,
                                             X.st = X.st,
                                             beta0.samples = cur.samples.beta0,
                                             tau2.samples = cur.samples.tau2,
                                             beta.samples = cur.samples.beta,
                                             r.samples = cur.samples.r,
                                             k.samples = cur.samples.k)$y.pred

      samples.spatInt_forecast.save[ls,,]<- pred.reg(nt = nt, ns = ns,
                                                     data_lik =  data_lik,
                                                     spatInt.ind = spatInt.ind,
                                                     forcast.ind = forcast.ind,
                                                     pred_type = "spatInt_forecast",
                                                     X.intpl = X.intpl,
                                                     X.frcast = X.frcast,
                                                     X.st = X.st,
                                                     beta0.samples = cur.samples.beta0,
                                                     tau2.samples = cur.samples.tau2,
                                                     beta.samples = cur.samples.beta,
                                                     r.samples = cur.samples.r,
                                                     k.samples = cur.samples.k)$y.pred
      ls<- ls + 1
    }


    ### Saving the samples after thinning the samples at every thiniterations .
    if((i%%thin)-1==0){
      samples[j,]<- c( cur.samples.r, cur.samples.k, cur.samples.beta0,
                       cur.samples.tau2, cur.samples.beta,
                      cur.samples.lambda[c(1,10,nt+10, nt+15)])


      j=j+1
    }
    if((i%%thin)==0 | i==1){  ## this is the thinning steps and adapt is the number of iterations after which i update the variance of MALA and random walk algorithms
      par(mfrow=c(5,4),oma=c(0,0,2,0),mar=c(4,5,1,1))
      if(print.result==TRUE){
        if(i< (burn_in1+burn_in2+2)){
          cat(paste0("Iteration: ",i, "\n",
                     " rMSE = ",round(sqrt(mean((Y.o -  y.imputed)^2, na.rm=TRUE)), digits = 5), " ",
                     "  |  MAE = ",round((mean(abs(Y.o -  y.imputed), na.rm=TRUE)), digits = 5), "\n",
                     " Accep rate lambda[1] = ", round(rate.lambda[1]/((i%%(adapt))), digits = 3),
                      " | tuning lambda[1]=",round(tun_lambda[1], digits = 5), "\n",
                     " Accep rate lambda[2] = ", round(rate.lambda[2]/((i%%(adapt))), digits = 3),
                     " | tuning lambda[2]=",round(tun_lambda[2], digits = 5), "\n",
                     " Accep rate r = ", if(data_lik=="NegB"){round(rate.r/((i%%(adapt))), digits = 3)} else{NULL},
                     " | tuning r=", if(data_lik=="NegB"){round(tun_r, digits = 5)} else{NULL}, "\n",
                     "-------------------------------------------------------------",  "\n",
                     sep=""))
        } else{
          cat(paste0("Iteration: ",i, "\n",
                     "  rMSE = ",round(sqrt(mean((Y.o -  y.imputed)^2, na.rm=TRUE)), digits = 5), " ",
                     " | MAE = ",round((mean(abs(Y.o -  y.imputed), na.rm=TRUE)), digits = 5), "\n",
                     " Accep rate lambda[1]=", round(rate.lambda[1]/(i-(burn_in1+burn_in2)), digits = 3),
                     " | tuning lambda[1] =", round(tun_lambda[1], digits = 5), "\n",
                     "  Accep rate lambda[2]=", round(rate.lambda[2]/(i-(burn_in1+burn_in2)), digits = 3),
                     " | tuning lambda[2] =", round(tun_lambda[2], digits = 5), "\n",
                     " Accep rate r = ", if(data_lik=="NegB"){round(rate.r/(i-(burn_in1+burn_in2)), digits = 3)} else{NULL},
                     " | tuning r=", if(data_lik=="NegB") {round(tun_r, digits = 5)} else{NULL}, "\n",
                     "-------------------------------------------------------------",  "\n",
                     sep=""))
        }
      }
      if(traceplot){
        model.param.name<- init.names$param.name
        for (lls  in 1: length(model.param.name)) {
          if(simulation==TRUE){
            plot(thin*c(0:(l-1))+1,samples[1:l,lls],type = "l",xlab="MCMC iteration",ylab=model.param.name[lls])
            abline(h=true.values[lls], col=2)
          }  else {
            plot(thin*c(0:(l-1))+1,samples[1:l,lls],type = "l",xlab="MCMC iteration",ylab=model.param.name[lls])

          }
        }
      }

      if((i%%adapt)-1==0){
        k<-k+1
      }

      l<-l+1
    }

      if((i%%adapt)-1==0){
      sigma.matrix[m,]<- c(tun_lambda[1], tun_lambda[2], tun_r)
      m=m+1
    }

  }

  return(
    list("samples.trace.plot" = list("samples.r" =  if(data_lik=="NegB"){samples[,1]} else{NULL},
                                          "samples.k" =  if(data_lik=="lNormal"){samples[,1]} else{NULL},
                                          "samples.beta0" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,2]} else{samples[,1]},
                                          "samples.tau2" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,3]} else{samples[,2]},
                                          "samples.beta" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[, 4:(3+q)]} else{samples[, 3:(2+q)]},
                                          "samples.few.lambdas" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[, (4+q):(7+q)]} else{samples[, (3+q):(6+q)]},
                                          "samples.kappa" = NULL,
                                          "samples.sigma2" = NULL,
                                          "sample.Wt" = NULL,
                                          "samples.few.mu.st" = NULL
  ),
  "within.sample.pred"=list("pred.cont.mean" = data.mean * (sum.post.Y.save/(N.MCMC - (burn_in1 + burn_in2))),
                            "pred.cont.sd"=  sqrt(data.mean^2 * (sum.post.Y2.save/(N.MCMC - (burn_in1 + burn_in2))) - data.mean^2 * (sum.post.Y.save/(N.MCMC - (burn_in1 + burn_in2)))^2),
                            "pred.lambda.mean" = data.log.mean * sum.post.lambda.save/(N.MCMC - (burn_in1 + burn_in2)),
                            "pred.lambda.sd" = data.log.mean^2 * (sum.post.lambda2.save/(N.MCMC - (burn_in1 + burn_in2))) - data.log.mean^2 * (sum.post.lambda.save/(N.MCMC - (burn_in1 + burn_in2)))^2
  ),
  "miss.value.imputation.result" = list("samples.miss.imput" = samples.miss.imput.Y.save
  ),
  "spat.temp.pred" = list("samples.spat.intpol" = samples.spat.pred.save,
                          "samples.future.forecast" = samples.forecast.save,
                          "samples.spat.future.forecast" = samples.spatInt_forecast.save
  ),
  "tun.param.info" = sigma.matrix)
  )

}



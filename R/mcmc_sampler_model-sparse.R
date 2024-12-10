#' MCMC Sampler for Sparse Spatio-Temporal Dynamic Generalized Linear Model
#'
#' This function performs MCMC sampling for a sparse spatio-temporal dynamic generalized linear model (DGLM) using a Gaussian random field approximation through SPDE. It handles spatial interpolation and forecasting and allows for various data likelihoods.
#'
#' @param data_lik Character, specifying the data likelihood (e.g., `"Poisson"`, `"NegB"`, `"lNormal"`, etc.).
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
#' @param Ft.o Matrix of dynamic temporal covariates for observed data.
#' @param F.mat.o 3D array of dynamic temporal covariates for observed data.
#' @param Ft.frcast Matrix of dynamic temporal covariates for forecasting.
#' @param m0.R Numeric, initial value for the SPDE spatial random effect.
#' @param C0.R Numeric, prior variance for the SPDE spatial random effect.
#' @param m0.theta Vector, initial value for dynamic temporal covariates.
#' @param C0.theta Matrix, prior covariance matrix for dynamic temporal covariates.
#' @param G Transition matrix for the state-space model.
#' @param delta Numeric, the maximum distance between the observed locations.
#' @param A.proj.o Projection matrix for the observed spatial locations.
#' @param A.proj.p Projection matrix for the predicted spatial locations.
#' @param c.mat.o Matrix representing the 0th-order basis functions for SPDE.
#' @param g1.mat.o Matrix representing the 1st-order basis functions for SPDE.
#' @param g2.mat.o Matrix representing the 2nd-order basis functions for SPDE.
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
#' @param tun_kappa Numeric, tuning parameter for `kappa` (spatial range parameter).
#' @param tun_lambda Numeric, tuning parameter for `lambda`.
#' @param print.result Logical, if `TRUE`, prints the MCMC results.
#' @param traceplot Logical, if `TRUE`, generates trace plots of MCMC samples.
#' @param true.values List of true parameter values (optional).
#' @param init_beta initial values for beta set by pre-run GLM counterpart
#' @param init_lambda initial values for linear predictors set by pre-run GLM counterpart
#' @param simulation Logical, if `TRUE`, the function runs in simulation mode.
#' @param init.seed Integer, seed for random number generation to ensure reproducibility.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A list containing MCMC samples, trace plots, convergence diagnostics, and predictions.
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
#' Ft.o <- matrix(rnorm(30), nrow = 10)
#' F.mat.o <- array(rnorm(150), dim = c(10, 5, 3))
#' Ft.frcast <- matrix(rnorm(6), nrow = 2)
#' m0.R <- 0
#' C0.R <- 1
#' m0.theta <- rep(0, 3)
#' C0.theta <- diag(5, 3)
#' G <- diag(1, 3)
#' delta <- 0.1
#' A.proj.o <- matrix(rnorm(30), nrow = 10)
#' A.proj.p <- matrix(rnorm(20), nrow = 5)
#' c.mat.o <- matrix(rnorm(25), nrow = 5)
#' g1.mat.o <- matrix(rnorm(25), nrow = 5)
#' g2.mat.o <- matrix(rnorm(25), nrow = 5)
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
#' tun_kappa <- 0.1
#' tun_lambda <- 0.1
#' print.result <- TRUE
#' traceplot <- TRUE
#' true.values <- NULL
#' simulation <- FALSE
#' init.seed <- 123
#' result <- MCMC.sampler_model.sparse(data_lik, nt, ns, p, q, Y.o, data.log.mean, data.mean,
#'                                     Y.intpl, Y.frcast, Y.st, NULL, X.o, quad.X.o, X.intpl,
#'                                     X.frcast, X.st, X.o, Ft.o, F.mat.o, Ft.frcast, m0.R, C0.R,
#'                                     m0.theta, C0.theta, G, delta, A.proj.o, A.proj.p,
#'                                     c.mat.o, g1.mat.o, g2.mat.o, spatInt.ind, forcast.ind,
#'                                     samples.store, N.MCMC, burn_in1, burn_in2, adapt_seq,
#'                                     thin, adapt, tun_r, tun_kappa, tun_lambda, print.result,
#'                                     traceplot, true.values, simulation, init.seed)
MCMC.sampler_model.sparse <- function(model,
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
                                     Ft.o,
                                     F.mat.o,
                                     Ft.frcast,
                                     m0.R,
                                     C0.R,
                                     m0.theta,
                                     C0.theta,
                                     G,
                                     delta,
                                     A.proj.o,
                                     A.proj.p,
                                     c.mat.o,
                                     g1.mat.o,
                                     g2.mat.o,
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
                                     tun_kappa,
                                     tun_lambda,
                                     print.result,
                                     traceplot,
                                     true.values,
                                     init_beta,
                                     init_lambda,
                                     simulation,
                                     init.seed,
                                     ...)
{
  #browser()
  nb<- ncol(A.proj.o)
  A.proj.SPDE<- as.matrix(A.proj.o)
  tA.proj<- as.matrix(t(A.proj.SPDE))
  t.A.A<- t(A.proj.SPDE) %*% A.proj.SPDE

  tun_lambda<- matrix(rep(tun_lambda, nt * ns), nrow = nt, ncol = ns)
  if(data_lik=="NegB"){tun_r= tun_r} else {tun_r=NULL}
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt),
                       ncol= 2  + length(tun_r) + 1,
                       byrow=TRUE) ## storing the adaptive scale parameters ## storing the adaptive scale parameters
  sigma.matrix[1,]<-c(tun_lambda[1], tun_lambda[2], tun_kappa, tun_r)

  init.names<- init_fun(p=p, q=q, ns = ns, nt = nt, nb=nb, Y=Y.o, model=model,
                        data_lik = data_lik, delta=delta, seed=init.seed)

  cur.samples.r<- init.names$init$r
  cur.samples.k<- init.names$init$k
  cur.samples.kappa<- init.names$init$kappa
  cur.samples.sigma2<- init.names$init$sigma2
  cur.samples.tau2<- init.names$init$tau2
  cur.samples.Wt<- init.names$init$Wt
  cur.samples.R<- init.names$init$R.st
  cur.samples.beta<- init_beta # init.names$init$beta
  cur.samples.theta<- init.names$init$theta
  cur.samples.lambda<- init_lambda # init.names$init$lambda_init

  samples.theta.save<- cur.samples.theta
  ### saving the Chain for all the hyper parameters and some latent parameters
  samples <- matrix(nrow=floor(N.MCMC/thin), ncol=  length(cur.samples.r) + length(cur.samples.k) +
                      3 + p + 3 + q + 4 + 4, byrow=TRUE) # storing the samples
  samples[1,]<- c(cur.samples.r, cur.samples.k, cur.samples.kappa,
                  cur.samples.sigma2, cur.samples.tau2, cur.samples.Wt,
                  cur.samples.R[1:3], cur.samples.beta,
                  cur.samples.theta[1:4], cur.samples.lambda[c(1, 3, nt+10, nt+15)]) ## saving only few samples

  #### saving the samples of the few parameters for the prediction later
  no.samples<- samples.store
  samples.save<-  floor(seq(from=burn_in1+burn_in2+1, to=N.MCMC, length.out=no.samples))
  samples.save[no.samples]<- samples.save[no.samples]-1 ### last samples
  samples.miss.imput.Y.save<- array(NA, dim = c(no.samples, sum(ind_NA_Y)))
  samples.spat.pred.save<- array(NA, dim = c(no.samples, nt, length(spatInt.ind)))
  samples.forecast.save<- array(NA, dim = c(no.samples, length(forcast.ind), ns))
  samples.spatInt_forecast.save<- array(NA, dim = c(no.samples, length(forcast.ind), length(spatInt.ind)))


  ### assessing within sample predictions
  sum.post.Y.save<- rep(0, sum(!ind_NA_Y))
  sum.post.Y2.save<- rep(0, sum(!ind_NA_Y))
  sum.post.lambda.save<- matrix(0, nrow = nt, ncol = ns)
  sum.post.lambda2.save<- matrix(0, nrow = nt, ncol = ns)

  sum.post.R.save<- matrix(0, nrow = nt, ncol=nb)
  sum.post.R2.save<- matrix(0, nrow = nt, ncol=nb)
  sum.post.mu.st.save<- matrix(0, nrow = nt, ncol=ns)
  sum.post.mu.st2.save<- matrix(0, nrow = nt, ncol=ns)

  ### run-time for kappas and updating mus
  run_time_kappa<- 0
  run_time_Rs<- 0

  Ft.mean<- matrix(rowSums(Ft.o * cur.samples.theta), nrow = nt, ncol=ns)
  comp.cov<- matrix(c(X.o %*% cur.samples.beta), nrow = nt, ncol = ns)
  cor.mat.inv.and.log.det<- cormat.inv.update.inla(rho = cur.samples.kappa, c.mat=c.mat.o,
                                                   g1.mat=g1.mat.o, g2.mat=g2.mat.o, alpha = 2)
  cor.mat.inv= cor.mat.inv.and.log.det$cormat.inv
  log.det.cor.mat = cor.mat.inv.and.log.det$cormat.logdet

  #browser()

  j<-1
  l<-1
  m<-1
  k<-1
  ls<- 1
  ls.w<- 1
  for (i in 1:(N.MCMC-1)){
    if(((i%%(adapt))-1==0) & (i< (burn_in1+burn_in2+2))){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
      rate.kappa<- 0
      rate.sigma2<- 0
      rate.lambda<- matrix(0, nrow = nt, ncol = ns)
      rate.r<- 0
    }

    # browser()
    #
    # print(i)
    ######## imputations #########
      if(data_lik=="Poisson"){
        y.imputed<- matrix(rpois(n  = ns * nt, lambda = exp(cur.samples.lambda)), nrow = nt, ncol=ns)
      } else if(data_lik=="NegB"){
        y.imputed<- matrix(rnbinom(n  = ns * nt, size = cur.samples.r, mu = exp(cur.samples.lambda)), nrow = nt, ncol=ns)
      } else if(data_lik =="lNormal"){
        y.imputed<-  matrix(exp(rnorm(n = ns * nt, mean = cur.samples.lambda, sd = sqrt(cur.samples.k))) -1, nrow = nt, ncol=ns)
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


    ###### update kappa ##########
    start_time_kappa<- Sys.time()
    cur.samples.kappa.all<- update_kappa.sparse(nt=nt,
                                                cur.kappa = cur.samples.kappa,
                                                sigma2 = cur.samples.sigma2,
                                                c.mat = c.mat.o,
                                                g1.mat = g1.mat.o,
                                                g2.mat = g2.mat.o,
                                                R = cur.samples.R,
                                                R0=m0.R,
                                                cor.mat.inv = cor.mat.inv,
                                                log.det.cor.mat= log.det.cor.mat,
                                                delta = delta ,
                                                tun.kappa =  tun_kappa
                                                )

    cor.mat.inv= cur.samples.kappa.all$cor.mat.inv
    log.det.cor.mat = cur.samples.kappa.all$log.det.cor.mat
    rate.kappa<- ifelse(cur.samples.kappa.all$ind, rate.kappa+1 , rate.kappa)
    tun_kappa<- adpative_function(index_MCMC_iter=i, sigma2_adapt=tun_kappa, target_accept=0.40,
                                  rate_adapt=rate.kappa, burn_in1=burn_in1, burn_in2=burn_in2, adapt_seq=adapt_seq,
                                  adapt=adapt, adpat_param=1, lower.acc=0.35, upper.acc=0.50)

    cur.samples.kappa<- cur.samples.kappa.all$cur.kappa
    end_time_kappa<- Sys.time()
    run_time_kappa<- run_time_kappa + (end_time_kappa - start_time_kappa)


    #### update sigma2 #############
    mu_diff<- cur.samples.R[2:nt,] - cur.samples.R[1:(nt-1),]
    qud_sum<- t(cur.samples.R[1,] - m0.R) %*% cor.mat.inv %*% (cur.samples.R[1,] - m0.R) +  sum((mu_diff %*% cor.mat.inv) * mu_diff)
    cur.samples.sigma2<- update_var(ns=nb,
                                    nt=nt,
                                    param.hyperprior=c(0.25, 0.25),
                                    qud_sum = as.numeric(qud_sum)
    )

    #cur.samples.sigma2<- hyper$sigma2
    #### update tau2 #############
    cur.samples.tau2<- update_var(ns=ns,
                                  nt=nt,
                                  param.hyperprior=c(0.25, 0.25),
                                  qud_sum = sum((cur.samples.lambda - cur.samples.R %*% tA.proj -
                                                   comp.cov - Ft.mean)^2)
    )


    #### update spatial-temporal intercepts #############
    start_time_Rs<- Sys.time()
    cur.samples.R<- update_intercept.sparse(nt = nt,
                                                 ns = ns,
                                                 nb = nb,
                                                 mu_current = cur.samples.R,
                                                 R0=m0.R,
                                                 lambdas_mean = cur.samples.lambda - Ft.mean - comp.cov,
                                                 Q = cor.mat.inv,
                                                 t.A.A = t.A.A,
                                                 t.A = A.proj.SPDE,
                                                 sigma2 =cur.samples.sigma2,
                                                 tau2 = cur.samples.tau2)

    end_time_Rs<- Sys.time()
    run_time_Rs<- run_time_Rs + (end_time_Rs - start_time_Rs)

   # cur.samples.R<- sim_M0$R.st[-(forcast.ind),]
    ### update betas
    cur.samples.beta<- update_space_time_beta(lambda = cur.samples.lambda,
                                              m = cur.samples.R %*% tA.proj + Ft.mean,
                                              X = X.o,
                                              quad.X = quad.X.o,
                                              tau2 = cur.samples.tau2
    )
    comp.cov<-  matrix(c(X.o%*% cur.samples.beta), nrow = nt, ncol = ns)


    #### update Wt ############
    cur.samples.theta0<- m0.theta
    cur.samples.Wt<- update.Wt(nt = nt,
                               p =p,
                               theta0 =cur.samples.theta0,
                               theta = cur.samples.theta,
                               G = G,
                               param.hyperprior = c(2,0.1))

    ########## update theta using FFBS #########

    cur.samples.theta<- ffbs_thetas(p=p,
                                        mus = cur.samples.lambda - cur.samples.R %*% tA.proj - comp.cov,
                                        G=G,
                                        Sigma= diag(cur.samples.Wt, p),
                                        H= F.mat.o,
                                        tau2 = cur.samples.tau2,
                                        ns=ns,
                                        nt=nt)

    Ft.mean<-  matrix(rowSums(Ft.o * cur.samples.theta), nrow = nt, ncol=ns)


    ########## update lambda
    cur.samples.lambda.all<- update_lambda_MMALA(ns = ns,
                                                 nt=nt,
                                                 Y = Y.o,
                                                 cur.lambda = cur.samples.lambda,
                                                 mean.lambda = cur.samples.R %*% tA.proj +  Ft.mean + comp.cov,
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

    #cur.samples.lambda<- sim_M0$lambda
    #browser()

    ## storing within samples prediction performances
    if(i > burn_in1 + burn_in2){
      sum.post.Y.save<- sum.post.Y.save + (y.imputed[!ind_NA_Y])/data.mean
      sum.post.Y2.save<- sum.post.Y2.save + ((y.imputed[!ind_NA_Y])^2)/(data.mean^2)
      sum.post.lambda.save<- sum.post.lambda.save + cur.samples.lambda/data.log.mean
      sum.post.lambda2.save<- sum.post.lambda2.save + ((cur.samples.lambda)^2)/ (data.log.mean^2)

      sum.post.R.save<- sum.post.R.save + cur.samples.R
      sum.post.R2.save<- sum.post.R2.save +  cur.samples.R^2

      sum.post.mu.st.save<- sum.post.mu.st.save + cur.samples.R %*% tA.proj
      sum.post.mu.st2.save<- sum.post.mu.st2.save + (cur.samples.R %*% tA.proj)^2

    }
    #### just in case posteriors means is inflating the results, save the samples to asses the within sample assessments and calculate medians


    if(i == samples.save[ls]){
      samples.miss.imput.Y.save[ls,]<- y.imputed[ind_NA_Y]
      samples.spat.pred.save[ls,,]<- pred.sparse(nt = nt, ns = ns, nb=nb,
                                                 data_lik = data_lik,
                                                            spatInt.ind = spatInt.ind,
                                                            forcast.ind = forcast.ind, pred_type = "spatInt",
                                                            X.intpl = X.intpl,
                                                            X.frcast = X.frcast,
                                                            X.st = X.st,
                                                            Ft.o = Ft.o,
                                                            Ft.frcast = Ft.frcast,
                                                            G = G,
                                                            A.proj.o = A.proj.o,
                                                            A.proj.p = A.proj.p,
                                                            c.mat = c.mat.o,
                                                            g1.mat = g1.mat.o,
                                                            g2.mat = g2.mat.o,
                                                            kappa.samples = cur.samples.kappa,
                                                            sigma2.samples = cur.samples.sigma2,
                                                            tau2.samples = cur.samples.tau2,
                                                            theta.samples = cur.samples.theta,
                                                            Wt.samples = cur.samples.Wt,
                                                            beta.samples = cur.samples.beta,
                                                            Rs.samples = cur.samples.R, r.samples = cur.samples.r,
                                                            k.samples = cur.samples.k)$y.pred

      samples.forecast.save[ls,,]<- pred.sparse(nt = nt, ns = ns, nb=nb,
                                                data_lik = data_lik,
                                                spatInt.ind = spatInt.ind,
                                                forcast.ind = forcast.ind, pred_type = "forecast",
                                                X.intpl = X.intpl,
                                                X.frcast = X.frcast,
                                                X.st = X.st,
                                                Ft.o = Ft.o,
                                                Ft.frcast = Ft.frcast,
                                                G = G,
                                                A.proj.o = A.proj.o,
                                                A.proj.p = A.proj.p,
                                                c.mat = c.mat.o,
                                                g1.mat = g1.mat.o,
                                                g2.mat = g2.mat.o,
                                                kappa.samples = cur.samples.kappa,
                                                sigma2.samples = cur.samples.sigma2,
                                                tau2.samples = cur.samples.tau2,
                                                theta.samples = cur.samples.theta,
                                                Wt.samples = cur.samples.Wt,
                                                beta.samples = cur.samples.beta,
                                                Rs.samples = cur.samples.R, r.samples = cur.samples.r,
                                                k.samples = cur.samples.k)$y.pred



      samples.spatInt_forecast.save[ls,,]<- pred.sparse(nt = nt, ns = ns, nb=nb,
                                                        data_lik = data_lik,
                                                        spatInt.ind = spatInt.ind,
                                                        forcast.ind = forcast.ind, pred_type = "spatInt_forecast",
                                                        X.intpl = X.intpl,
                                                        X.frcast = X.frcast,
                                                        X.st = X.st,
                                                        Ft.o = Ft.o,
                                                        Ft.frcast = Ft.frcast,
                                                        G = G,
                                                        A.proj.o = A.proj.o,
                                                        A.proj.p = A.proj.p,
                                                        c.mat = c.mat.o,
                                                        g1.mat = g1.mat.o,
                                                        g2.mat = g2.mat.o,
                                                        kappa.samples = cur.samples.kappa,
                                                        sigma2.samples = cur.samples.sigma2,
                                                        tau2.samples = cur.samples.tau2,
                                                        theta.samples = cur.samples.theta,
                                                        Wt.samples = cur.samples.Wt,
                                                        beta.samples = cur.samples.beta,
                                                        Rs.samples = cur.samples.R, r.samples = cur.samples.r,
                                                        k.samples = cur.samples.k)$y.pred
      ls<- ls + 1
    }

    #print(i)
    ### Saving the samples after thinning the samples at every thiniterations .
    if((i%%thin)-1==0){
      samples[j,]<- c(cur.samples.r, cur.samples.k,
                      cur.samples.kappa, cur.samples.sigma2, cur.samples.tau2, cur.samples.Wt,
                      cur.samples.R[1:3], cur.samples.beta,
                      cur.samples.theta[1:4], cur.samples.lambda[c(1,10,nt+10, nt+15)])

        samples.theta.save<- abind::abind(samples.theta.save, cur.samples.theta, along = 3)

      j=j+1
    }
    if((i%%thin)==0 | i==1){  ## this is the thinning steps and adapt is the number of iterations after which i update the variance of MALA and random walk algorithms
      par(mfrow=c(6,6),oma=c(0,0,2,0),mar=c(4,5,1,1))
      if(print.result==TRUE){
        if(i< (burn_in1+burn_in2+2)){
          cat(paste0("Iteration: ",i, "\n",
                     " rMSE = ",round(sqrt(mean((Y.o -  y.imputed)^2, na.rm=TRUE)), digits = 5), " ",
                     " | MAE = ",round((mean(abs(Y.o -  y.imputed), na.rm=TRUE)), digits = 5), "\n",
                     " Accep rate kappa = ", round(rate.kappa/((i%%(adapt))), digits = 3),
                     " | tuning kappa = ",round(tun_kappa, digits = 5), "\n",
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
                     " Accep rate kappa =", round(rate.kappa/(i-(burn_in1+burn_in2)), digits = 3),
                     " | tuning kappa =", round(tun_kappa, digits = 5), "\n",
                     "  Accep rate lambda[1]=", round(rate.lambda[1]/(i-(burn_in1+burn_in2)), digits = 3),
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
            plot(thin*c(0:(l-1))+1,samples[1:l,lls],type = "l",xlab="MCMC iteration",ylab=model.param.name[lls]) # Plot for alpha.tilde
            abline(h=true.values[lls], col=2)
          }  else {
            plot(thin*c(0:(l-1))+1,samples[1:l,lls],type = "l",xlab="MCMC iteration",ylab=model.param.name[lls]) # Plot for alpha.tilde

          }
        }
      }

      l<-l+1
    }


    ### storing the adaptive tuning parameters
    if((i%%adapt)-1==0){ # to save allexp the scale parameter of the MALA
      sigma.matrix[m,]<- c(tun_lambda[1], tun_lambda[2], tun_r, tun_kappa)
      m=m+1
    }

  }
  return(list("samples.trace.plot" = list("samples.r" =  if(data_lik=="NegB"){samples[,1]} else{NULL},
                                           "samples.k" =  if(data_lik=="lNormal"){samples[,1]} else{NULL},
                                           "samples.kappa" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,2]} else{samples[,1]},
                                           "samples.sigma2" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,3]} else{samples[,2]},
                                           "samples.tau2" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,4]} else{samples[,3]},
                                           "sample.Wt" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,5:(4+p)]} else{samples[, 4:(3+p)]},
                                           "samples.few.R.st" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,(5+p):(7+p)]} else{ samples[, (4+p):(6+p)]},
                                           "samples.beta" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,(8+p):(7+p+q)]} else{samples[,(7+p):(6+p +q)]},
                                           "samples.few.lambdas" = if(data_lik=="NegB" | data_lik=="lNormal"){samples[,(4+8+p+q):(4+11+p+q)]} else{ samples[,(4+7+p+q):(4+10+p+q)]},
                                           "samples.beta0" = NULL
  ),
  "samples.theta" = samples.theta.save,
  "within.sample.pred"=list("pred.cont.mean" = data.mean * (sum.post.Y.save/(N.MCMC - (burn_in1 + burn_in2))),
                            "pred.cont.sd"=  sqrt(data.mean^2 * (sum.post.Y2.save/(N.MCMC - (burn_in1 + burn_in2))) - data.mean^2 * (sum.post.Y.save/(N.MCMC - (burn_in1 + burn_in2)))^2),
                            "pred.lambda.mean" = data.log.mean * sum.post.lambda.save/(N.MCMC - (burn_in1 + burn_in2)),
                            "pred.lambda.sd" = data.log.mean^2 * (sum.post.lambda2.save/(N.MCMC - (burn_in1 + burn_in2))) - data.log.mean^2 * (sum.post.lambda.save/(N.MCMC - (burn_in1 + burn_in2)))^2,
                            "pred.R.mean" =  sum.post.R.save / (N.MCMC - (burn_in1 + burn_in2)),
                            "pred.R.sd" = sqrt(sum.post.R2.save / (N.MCMC - (burn_in1 + burn_in2)) - (sum.post.R.save / (N.MCMC - (burn_in1 + burn_in2)))^2),
                            "pred.mu.st.mean" =  sum.post.mu.st.save / (N.MCMC - (burn_in1 + burn_in2)),
                            "pred.mu.st.sd" = sqrt(sum.post.mu.st2.save / (N.MCMC - (burn_in1 + burn_in2)) - (sum.post.mu.st.save / (N.MCMC - (burn_in1 + burn_in2)))^2)

  ),
  "miss.value.imputation.result" = list("samples.miss.imput" = samples.miss.imput.Y.save
  ),
  "spat.temp.pred" = list("samples.spat.intpol" = samples.spat.pred.save,
                          "samples.future.forecast" = samples.forecast.save,
                          "samples.spat.future.forecast" = samples.spatInt_forecast.save
  ),
  "tun.param.info" = sigma.matrix,
  "runtime" = list("spat_range" =  as.numeric(run_time_kappa)/(N.MCMC),
                   "spat_intercept" = as.numeric(run_time_Rs)/(N.MCMC))
  ))


}



#' Forward Filtering Backward Sampling (FFBS) for State Vectors (Thetas)
#'
#' This function implements the Forward Filtering Backward Sampling (FFBS) algorithm to draw samples of state vectors (thetas) from their posterior distribution.
#'
#' @param p Integer, the dimension of the state vector.
#' @param mus A matrix of dimension (nt x ns), the mean of the observations.
#' @param m0 Initial mean vector of the state vector at time t = 0.
#' @param m0.bar Posterior mean vector of the state vector at time t = 0.
#' @param G A matrix of dimension (p x p), the state-evolution matrix.
#' @param Sigma A covariance matrix of the state vector.
#' @param H Observation matrix, mapping state vectors to observations.
#' @param tau2 Observation noise variance.
#' @param ns Integer, number of spatial locations.
#' @param nt Integer, number of time points.
#'
#' @return A matrix of sampled state vectors (thetas) of dimension (nt x p).
#' @export
#'
#' @examples
#' \dontrun{
#' p <- 2
#' ns <- 5
#' nt <- 100
#' mus <- matrix(runif(nt * ns), nrow = nt, ncol = ns)
#' m0 <- rep(0, p)
#' m0.bar <- rep(0, p)
#' G <- diag(p)
#' Sigma <- diag(p)
#' H <- matrix(runif(ns * p), nrow = ns, ncol = p)
#' tau2 <- 0.1
#' theta_samples <- ffbs_thetas(p, mus, m0, m0.bar, G, Sigma, H, tau2, ns, nt)
#' }
#' @examples
ffbs_thetas<- function (p,
                        mus,
                        m0,
                        m0.bar,
                        G,
                        Sigma,
                        H,
                        tau2,
                        ns,
                        nt)
{
##browser()
  G <- as.matrix(G)
  mtt1 <- array(0, c(nt, p))
  mtt <- array(0, c(nt + 1, p))
  Rtt1 <- array(0, c(nt, p, p))
  Rtt <- array(0, c(nt + 1, p, p))
  simAlpha <- array(0, c(nt + 1, p))
  Innt <- array(0, c(nt, ns))
  PrecInnt <- array(0, c(nt, ns, ns))
  Rtt[1, , ] <- as.matrix(Sigma)
  for (t in 1:nt) {
    mtt1[t, ] <- G %*% mtt[t, ]
    Rtt1[t, , ] <- Sigma + G %*% Rtt[t, , ] %*% t(G)
    Rtt1[t, , ] <- (Rtt1[t, , ] + t(Rtt1[t, , ]))/2
    TM1 <- matrix(c(H[t,1,] %*% Rtt1[t, , ]), nrow = ns, ncol=p, byrow = TRUE)   ## outer(rep(1,N), (c(H[t,1,] %*% Rtt1[t, , ])));  H[t,,] %*% Rtt1[t, , ]
    TM1.1 <-  c(matrix(H[t,1,], ncol=p) %*% Rtt1[t, , ] %*% t(matrix(H[t,1,], ncol=p)))
    Mti<- inverse_SMW(alpha =TM1.1, tau = tau2, n=ns) # solve(diag(tau2, N) + TM1 %*% t(H[t,,]))

    PrecInnt[t, , ] <- Mti
    TM2 <- Rtt1[t, , ] %*% outer(c(H[t,1,]), c(colSums(Mti))) # Rtt1[t, , ] %*% t(H[t,,]) %*% Mti
    Innt[t, ] <- mus[t,] - sum(H[t,1,] * mtt1[t, ])  ## mus[t,] - H[t,,] %*% mtt1[t, ]
    mtt[t + 1, ] <- mtt1[t, ] + TM2[,1] * sum(Innt[t, ])  ## mtt1[t, ] + TM2 %*% (Innt[t, ])
    Rtt[t + 1, , ] <- Rtt1[t, , ] - ns * outer(TM2[,1], TM1[1,]) ## Rtt1[t, , ] - TM2 %*% TM1
    Rtt[t + 1, , ] <- (Rtt[t + 1, , ] + t(Rtt[t + 1, , ]))/2
  }
  #browser()
    Rtt[nt + 1, , ] = (Rtt[nt + 1, , ] + t(Rtt[nt + 1, , ]))/2

    #St.pd<- Matrix::nearPD(Rtt[T + 1, , ])
    simAlpha[nt + 1, ] <- mvtnorm::rmvnorm(1, mean = mtt[nt + 1, ],
                                          sigma = Rtt[nt + 1, , ], method = "chol")
    for (t in nt:1) {
      tm <- Rtt[t, , ] %*% t(G) %*% solve(Rtt1[t, , ])
      Rt <- Rtt[t, , ] - tm %*% G %*% Rtt[t, , ]
      Rt <- (Rt + t(Rt))/2
      mt <- mtt[t, ] + tm %*% (simAlpha[t + 1, ] - mtt1[t, ])

      #Rt.pd<- Matrix::nearPD(Rt)
      simAlpha[t, ] <- mvtnorm::rmvnorm(1, mean = mt, sigma = Rt, #as.matrix(Rt.pd$mat),
                                        method = "chol")
    }
  retsimAlpha = simAlpha[-1, ]
  return(retsimAlpha)
}


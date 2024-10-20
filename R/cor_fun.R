#' Matern Correlation Function for Different Choices of Smoothness Parameters
#'
#' This function computes the Matern correlation for various choices of smoothness parameters based on the input correlation type.
#'
#' @param dist.mat A numeric matrix representing the distance matrix of dimensions n x n, where n is the total number of spatial locations.
#' @param kappa A numeric value representing the range parameter.
#' @param cor.type A character string specifying the type of Matern correlation. There are 5 possible choices: "Matern0.5", "Matern1", "Matern1.5", "Matern2.5", and "MaternInf".
#'
#' @return A numeric matrix of the same dimensions as `dist.mat`, containing the computed correlations.
#'
#' @details
#' - "Matern0.5" corresponds to the exponential correlation.
#' - "Matern1.5" and "Matern2.5" correspond to specific Matern correlations with different smoothness.
#' - "MaternInf" is also known as the squared exponential correlation, which is obtained as the smoothness parameter tends to infinity.
#' - "Matern1" includes the Bessel function for correlation computation, and a small value (1e-10) is added to avoid numerical instability.
#'
#' @export
#'
#' @examples
#' dist.mat <- matrix(runif(100), nrow = 10, ncol = 10)
#' kappa <- 1
#' cor.type <- "Matern1.5"
#' cor.fun(cor.type, dist.mat, kappa)
cor.fun<- function(cor.type, dist.mat, kappa){

  dist.mat<- dist.mat/kappa
  if(cor.type=="Matern0.5"){
    corr<-  exp(-dist.mat)
  } else if(cor.type=="Matern1.5"){
    corr<- (1 + dist.mat) * exp(-dist.mat)
  } else if(cor.type=="Matern2.5"){
    corr<- (1 + dist.mat + dist.mat^2/3) * exp(-dist.mat)
  }  else if(cor.type=="MaternInf"){ ## also called squared exponential and obtained as the smoothness parameter tends to infinity
    corr<- exp(-dist.mat^2)
  } else if(cor.type=="Matern1"){
    dist.mat[dist.mat==0]<- 1e-10  ## to avoid any numerical instability in calculation of bessel function
    corr<- dist.mat  * besselK(dist.mat, 1)
  }
  return(corr)
}


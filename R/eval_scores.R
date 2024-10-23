#' Compute Summary (Median) Statistics for Prediction Accuracy
#'
#' This function calculates various summary statistics to evaluate the accuracy of predictions 
#' using different scoring rules. It computes statistics like mean squared prediction error, 
#' mean absolute error, continuous ranked probability score (CRPS), and more.
#'
#' @param pred_val A vector of predicted values.
#' @param true_val A vector of true values, same length as \code{pred_val}.
#' @param std_pred A vector of standard deviations of the predicted values, same length as \code{pred_val}.
#' @param sumry_type A string specifying the type of summary statistic to compute. 
#'        Options include:
#'        \itemize{
#'          \item \code{"rmspe"}: Root mean squared prediction error.
#'          \item \code{"mae"}: Mean absolute error.
#'          \item \code{"IS"}: Interval score.
#'          \item \code{"mpiw"}: Mean prediction interval width.
#'          \item \code{"picp"}: Prediction interval coverage probability.
#'          \item \code{"normlized_PI"}: Normalized prediction residuals.
#'          \item \code{"crps"}: Continuous ranked probability score.
#'          \item \code{"norm_pr"}: Normalized prediction residuals.
#'        }
#' @param alpha A numeric value for the confidence level, default is 0.05 for a 95\% coverage probability.
#'
#' @return The calculated summary statistic based on the specified \code{sumry_type}.
#'
#' @details
#' The summary statistics are used to evaluate the performance of the predictive distribution 
#' based on the true observations. All the scores are negatively oriented, meaning lower values 
#' indicate better performance. 
#'
#' The different scoring methods include:
#' \itemize{
#'   \item \code{mspe or MSE}: Mean squared error, \code{S(F, y) = (y - m)^2}.
#'   \item \code{mae}: Mean absolute error, \code{S(F, y) = |y - m|}.
#'   \item \code{normlized_PI or IGN}: Ignorance score, \code{S(F, y) = (y - m)^2 / s^2 / 2 + log(s)}.
#'   \item \code{crps or CRPS}: Continuous ranked probability score.
#'   \item \code{IS}: Interval score.
#'   \item \code{mpiw}: Mean prediction interval width.
#'   \item \code{picp or Coverage}: Prediction interval coverage, \code{I(y in [m - q*s, m + q*s])}, 
#'         where \code{q} is the quantile corresponding to \code{1 - alpha/2}.
#' }
#'
#' For more details on scoring rules, see Gneiting & Raftery, JASA, 2007.
#'
#' @note
#' Be cautious of using improper scores. For example, \code{(y - m)^2 + s^2} is not a proper score 
#' and should not be used.
#'
#' @examples
#' # Example usage of eval_scores_median
#' pred_val <- c(1.2, 2.3, 3.1)
#' true_val <- c(1.0, 2.0, 3.0)
#' std_pred <- c(0.1, 0.2, 0.3)
#' eval_scores_median(pred_val, true_val, std_pred, sumry_type = "rmspe")
#'
#' @export
eval_scores_median<- function(pred_val, 
                              true_val, 
                              std_pred,
                              sumry_type="rmspe",
                              alpha=0.05
){
  #browser()
  if(sumry_type=="rmspe"){ ## Compute the sqerr-score: or mean squared prediction scores
    z <- as.numeric((true_val - pred_val))
    scores <- sqrt(median(z^2, na.rm=TRUE))
  } else if(sumry_type=="mae"){   ## Compute the absolute error score:
    z <- as.numeric((true_val - pred_val))
    scores <- median(abs(z), na.rm=TRUE)
  } else if(sumry_type=="IS"){  ## Compute the median prediction Interval score:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- median(2 * hw + (2/alpha) * (((pred_val - hw) - true_val) * (true_val < pred_val - hw) +
                                             (true_val - (pred_val + hw)) * (true_val > pred_val + hw)), na.rm=TRUE)
  } else if(sumry_type=="mpiw"){  ## Compute the median prediction Interval score:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- median(2 * hw, na.rm=TRUE)
  } else if(sumry_type=="picp"){ ## Compute the median prediction coverage:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- median((pred_val - hw <= true_val) & (true_val <= pred_val + hw), na.rm=TRUE)
  } else if (sumry_type=="normlized_PI") { ## Compute the Compute normalised prediction residuals
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <- median(z^2 / 2 + log(std_pred), na.rm=TRUE)
  } else if (sumry_type=="crps") {  ## Compute the crps-score:
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <-  median(std_pred * (z *(2 * pnorm(z, 0, 1) - 1) +
                                    2 * dnorm(z, 0, 1) - 1/sqrt(pi)), na.rm=TRUE)
  } else { ## Compute normalised prediction residuals
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <- median(pnorm(z, 0, 1), na.rm=TRUE)
  }
  return(scores)
}



#' Compute Summary(Mean) Statistics for Prediction Accuracy Using Mean
#'
#' This function computes various summary statistics to evaluate the accuracy of predictions 
#' using different scoring rules, based on the mean of the errors. It calculates statistics like 
#' mean squared prediction error, mean absolute error, continuous ranked probability score (CRPS), 
#' and more.
#'
#' @param pred_val A vector of predicted values.
#' @param true_val A vector of true values, same length as \code{pred_val}.
#' @param std_pred A vector of standard deviations of the predicted values, same length as \code{pred_val}.
#' @param sumry_type A string specifying the type of summary statistic to compute. 
#'        Options include:
#'        \itemize{
#'          \item \code{"rmspe"}: Root mean squared prediction error.
#'          \item \code{"mae"}: Mean absolute error.
#'          \item \code{"IS"}: Interval score.
#'          \item \code{"mpiw"}: Mean prediction interval width.
#'          \item \code{"picp"}: Prediction interval coverage probability.
#'          \item \code{"normlized_PI"}: Normalized prediction residuals.
#'          \item \code{"crps"}: Continuous ranked probability score.
#'          \item \code{"norm_pr"}: Normalized prediction residuals.
#'        }
#' @param alpha A numeric value for the confidence level, default is 0.05 for a 95\% coverage probability.
#'
#' @return The calculated summary statistic based on the specified \code{sumry_type}.
#'
#' @details
#' The summary statistics are used to evaluate the performance of the predictive distribution 
#' based on the true observations. All the scores are negatively oriented, meaning lower values 
#' indicate better performance. 
#'
#' The different scoring methods include:
#' \itemize{
#'   \item \code{mspe or MSE}: Mean squared error, \code{S(F, y) = (y - m)^2}.
#'   \item \code{mae}: Mean absolute error, \code{S(F, y) = |y - m|}.
#'   \item \code{normlized_PI or IGN}: Ignorance score, \code{S(F, y) = (y - m)^2 / s^2 / 2 + log(s)}.
#'   \item \code{crps or CRPS}: Continuous ranked probability score.
#'   \item \code{IS}: Interval score.
#'   \item \code{mpiw}: Mean prediction interval width.
#'   \item \code{picp or Coverage}: Prediction interval coverage, \code{I(y in [m - q*s, m + q*s])}, 
#'         where \code{q} is the quantile corresponding to \code{1 - alpha/2}.
#' }
#'
#' For more details on scoring rules, see Gneiting & Raftery, JASA, 2007.
#'
#' @note
#' Be cautious of using improper scores. For example, \code{(y - m)^2 + s^2} is not a proper score 
#' and should not be used.
#'
#' @examples
#' # Example usage of eval_scores_mean
#' pred_val <- c(1.2, 2.3, 3.1)
#' true_val <- c(1.0, 2.0, 3.0)
#' std_pred <- c(0.1, 0.2, 0.3)
#' eval_scores_mean(pred_val, true_val, std_pred, sumry_type = "rmspe")
#'
#' @export
eval_scores_mean<- function(pred_val, 
                            true_val, 
                            std_pred,
                            sumry_type="rmspe",
                            alpha=0.05
){
  #browser()
  if(sumry_type=="rmspe"){ ## Compute the sqerr-score: or mean squared prediction scores
    z <- as.numeric((true_val - pred_val))
    scores <- sqrt(mean(z^2, na.rm=TRUE))
  } else if(sumry_type=="mae"){   ## Compute the absolute error score:
    z <- as.numeric((true_val - pred_val))
    scores <- mean(abs(z), na.rm=TRUE)
  } else if(sumry_type=="IS"){  ## Compute the mean prediction Interval score:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- mean(2 * hw + (2/alpha) * (((pred_val - hw) - true_val) * (true_val < pred_val - hw) +
                                           (true_val - (pred_val + hw)) * (true_val > pred_val + hw)), na.rm=TRUE)
  } else if(sumry_type=="mpiw"){  ## Compute the mean prediction Interval score:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- mean(2 * hw, na.rm=TRUE)
  } else if(sumry_type=="picp"){ ## Compute the mean prediction coverage:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- mean((pred_val - hw <= true_val) & (true_val <= pred_val + hw), na.rm=TRUE)
  } else if (sumry_type=="normlized_PI") { ## Compute the Compute normalised prediction residuals
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <- mean(z^2 / 2 + log(std_pred), na.rm=TRUE)
  } else if (sumry_type=="crps") {  ## Compute the crps-score:
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <-  mean(std_pred * (z *(2 * pnorm(z, 0, 1) - 1) +
                                  2 * dnorm(z, 0, 1) - 1/sqrt(pi)), na.rm=TRUE)
  } else { ## Compute normalised prediction residuals
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <- mean(pnorm(z, 0, 1), na.rm=TRUE)
  }
  return(scores)
}





#' Estimate the Average Annual Daily Bike (AADB)
#'
#' estimate for the Average Annual Daily Bike (AADB).
#'
#' @param data A numeric vector containing the data for which the AADB is to be estimated.
#' @return A numeric value representing the estimated AADB.
#' @export
#' @examples
#' # Example usage:
#' data <- c(1, 2, 3, 4, 5)
#' est_AADB_simple(data)
est_AADB_simple <- function(data) {
  aadb <- mean(data)  
  return(aadb)
}

#' Generate Harmonic Matrices for a Given Period and Number of Harmonics
#'
#' This function generates two matrices, `Ft.all` and `G.all`, used in harmonic modeling.
#' `Ft.all` contains alternating columns of 1's and 0's for each harmonic, and `G.all`
#' is a block diagonal matrix based on cosine and sine transformations for each harmonic
#' over the specified period.
#'
#' @param nt Integer. The number of time points for each harmonic component.
#' @param num_harmonics Integer. The number of harmonics to generate.
#' @param period Numeric. The period used for harmonic transformations in `G.all`.
#'
#' @return A list containing:
#'   \item{Ft.all}{A matrix with `nt` rows and `2 * num_harmonics` columns. Each pair of columns
#'                 contains alternating 1's and 0's for each harmonic component.}
#'   \item{G.all}{A block diagonal matrix with harmonic transformations based on the specified period.}
#'
#' @examples
#' # Example with 10 time points, 3 harmonics, and a period of 7
#' harmonics <- generate_harmonic_matrices(nt = 10, num_harmonics = 3, period = 7)
#'
#' # Access Ft.all and G.all
#' Ft.all <- harmonics$Ft.all
#' G.all <- harmonics$G.all
#'
#' # Print Ft.all and G.all
#' print(Ft.all)
#' print(G.all)
#'
#' @export
generate_harmonic_matrices <- function(nt, num_harmonics, period) {
  # Create Ft.all matrix with alternating 1 and 0 for each harmonic
  Ft.all <- do.call(cbind, lapply(1:num_harmonics, function(h) {
    cbind(rep(1, nt), rep(0, nt))
  }))
  colnames(Ft.all) <- unlist(lapply(1:num_harmonics, function(h) {
    c(paste0("harm", h, ".1"), paste0("harm", h, ".2"))
  }))
  G.all <- as.matrix(Matrix::bdiag(lapply(1:num_harmonics, function(h) {
    angle <- 2 * pi * h / period
    matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2, byrow = TRUE)
  })))
  
  list(Ft.all = Ft.all, G.all = G.all)
}


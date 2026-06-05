#' Create a matrix of any size with a random sample of infected starting cells
#'
#' @param row Number of rows for matrix
#' @param col Number of columns for matrix
#' @param start_inf Probability that a cell will start as infected
#'
#' @return Infection Matrix
#'
#' @examples
#' create_matrix(5,5)
#' create_matrix(15,15,0.25)
#' @export
create_matrix <- function(row, col, start_inf = 0.125) {
  matrix(
    sample(c(0, 1), row * col, replace = TRUE, prob = c(1 - start_inf, start_inf)),
    nrow = row,
    ncol = col
  )
}

#' Create a matrix of any size with only the four corners starting as infected cells
#'
#' @param row Number of rows for matrix
#' @param col Number of columns for matrix
#'
#' @return Infection Matrix
#'
#' @examples
#' create_crnr_matrix(10,10)
#' @export
create_crnr_matrix <- function(row, col) {
  inf_matrix <- matrix(0, row, col)
  inf_matrix[1, 1] <- 1
  inf_matrix[1, col] <- 1
  inf_matrix[row, 1] <- 1
  inf_matrix[row, col] <- 1
  inf_matrix
}

#' Create a matrix of any size with only a single cell in the center starts as infected
#'
#' @param row Number of rows for matrix
#' @param col Number of columns for matrix
#'
#' @return Infection Matrix
#'
#' @examples
#' create_cntr_matrix(10,10)
#' @export
create_cntr_matrix <- function(row, col) {
  inf_matrix <- matrix(0, nrow = row, ncol = col)
  inf_matrix[ceiling(row / 2), ceiling(col / 2)] <- 1
  inf_matrix
}

#' Pad the input matrix so we don't have to worry about the edges
#'
#' @param input Input matrix
#'
#' @return input_pad
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' pad(x)
#' @noRd
pad <- function(input) {
  nr <- nrow(input)
  nc <- ncol(input)
  input_pad <- matrix(SUSCEPTIBLE, nrow = nr + 2, ncol = nc + 2)
  input_pad[2:(nr + 1), 2:(nc + 1)] <- input
  input_pad
}

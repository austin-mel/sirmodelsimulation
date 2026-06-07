#' Create a matrix of any size with random infected starting cells
#'
#' @param row Number of rows for matrix
#' @param col Number of columns for matrix
#' @param start_infected Probability that each cell starts infected when between
#'   0 and 1, or the exact number of randomly selected infected cells when 1 or
#'   greater.
#' @param seed Optional random seed for reproducible sampling
#'
#' @return Infection Matrix
#'
#' @examples
#' create_random_matrix(5, 5)
#' create_random_matrix(15, 15, 0.25)
#' create_random_matrix(15, 15, 10)
#' @export
create_random_matrix <- function(row, col, start_infected = 0.125, seed = NULL) {
  total_cells <- row * col

  if (!is.numeric(start_infected) || length(start_infected) != 1 ||
      is.na(start_infected) || start_infected < 0 || start_infected > total_cells) {
    stop("start_infected must be one number between 0 and row * col.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (start_infected == 0) {
    return(matrix(0, nrow = row, ncol = col))
  }

  if (start_infected < 1) {
    return(matrix(
      sample(c(0, 1), total_cells, replace = TRUE, prob = c(1 - start_infected, start_infected)),
      nrow = row,
      ncol = col
    ))
  }

  if (start_infected != floor(start_infected)) {
    stop("start_infected must be a whole number when it is 1 or greater.", call. = FALSE)
  }

  inf_matrix <- matrix(0, nrow = row, ncol = col)
  infected_cells <- sample(seq_len(total_cells), start_infected)
  inf_matrix[infected_cells] <- 1
  inf_matrix
}

#' Create a matrix of any size with random infected starting cells
#'
#' `create_matrix()` is kept for backward compatibility. Prefer
#' [create_random_matrix()] for new code.
#'
#' @param row Number of rows for matrix
#' @param col Number of columns for matrix
#' @param start_inf Probability that each cell starts infected when between 0
#'   and 1, or the exact number of randomly selected infected cells when 1 or
#'   greater.
#' @param seed Optional random seed for reproducible sampling
#'
#' @return Infection Matrix
#'
#' @examples
#' create_matrix(5, 5)
#' create_matrix(15, 15, 0.25)
#' @export
create_matrix <- function(row, col, start_inf = 0.125, seed = NULL) {
  create_random_matrix(row, col, start_infected = start_inf, seed = seed)
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
#' x <- create_random_matrix(5, 5)
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

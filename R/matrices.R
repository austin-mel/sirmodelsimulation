#' Create a matrix of any size with random infected starting cells
#'
#' @param row Number of rows for matrix
#' @param col Number of columns for matrix
#' @param start_infected Exact number of randomly selected cells that start infected
#' @param seed Optional random seed for reproducible sampling
#'
#' @return Infection Matrix
#'
#' @examples
#' create_random_matrix(5, 5)
#' create_random_matrix(15, 15, 10)
#' @export
create_random_matrix <- function(row, col, start_infected = 1, seed = NULL) {
  total_cells <- row * col

  if (!is.numeric(start_infected) || length(start_infected) != 1 ||
      is.na(start_infected) || start_infected < 0 || start_infected > total_cells) {
    stop("start_infected must be one number between 0 and row * col.", call. = FALSE)
  }

  if (start_infected != floor(start_infected)) {
    stop("start_infected must be a whole number.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (start_infected == 0) {
    return(matrix(0, nrow = row, ncol = col))
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
#' @param start_inf Exact number of randomly selected cells that start infected
#' @param seed Optional random seed for reproducible sampling
#'
#' @return Infection Matrix
#'
#' @examples
#' create_matrix(5, 5)
#' create_matrix(15, 15, 10)
#' @export
create_matrix <- function(row, col, start_inf = 1, seed = NULL) {
  create_random_matrix(row, col, start_infected = start_inf, seed = seed)
}

#' Create a matrix of any size with only the four corners starting as infected cells
#'
#' @param row Number of rows for matrix
#' @param col Number of columns for matrix
#' @param start_infected Accepted for consistency with [create_random_matrix()] and ignored
#' @param seed Accepted for consistency with [create_random_matrix()] and ignored
#'
#' @return Infection Matrix
#'
#' @examples
#' create_corner_matrix(10, 10)
#' @export
create_corner_matrix <- function(row, col, start_infected = 1, seed = NULL) {
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
#' @param start_infected Accepted for consistency with [create_random_matrix()] and ignored
#' @param seed Accepted for consistency with [create_random_matrix()] and ignored
#'
#' @return Infection Matrix
#'
#' @examples
#' create_center_matrix(10, 10)
#' @export
create_center_matrix <- function(row, col, start_infected = 1, seed = NULL) {
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

#' Create a random starting matrix
#'
#' Create a matrix for simulation input where `0` means susceptible and `1`
#' means infected. The requested number of cells are randomly selected to start
#' infected.
#'
#' @param row Number of rows for the matrix.
#' @param col Number of columns for the matrix.
#' @param start_infected Exact number of randomly selected cells that start
#'   infected.
#' @param seed Optional random seed for reproducible sampling.
#'
#' @return A numeric matrix with susceptible cells encoded as `0` and infected
#'   cells encoded as `1`.
#'
#' @examples
#' create_random_matrix(5, 5, start_infected = 1, seed = 123)
#' create_random_matrix(15, 15, start_infected = 10, seed = 123)
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

#' Create a random starting matrix with the legacy argument name
#'
#' `create_matrix()` is kept for backward compatibility. Prefer
#' [create_random_matrix()] for new code.
#'
#' @param row Number of rows for the matrix.
#' @param col Number of columns for the matrix.
#' @param start_inf Exact number of randomly selected cells that start infected.
#' @param seed Optional random seed for reproducible sampling.
#'
#' @return A numeric matrix with susceptible cells encoded as `0` and infected
#'   cells encoded as `1`.
#'
#' @examples
#' create_matrix(5, 5, start_inf = 1, seed = 123)
#' create_matrix(15, 15, start_inf = 10, seed = 123)
#' @export
create_matrix <- function(row, col, start_inf = 1, seed = NULL) {
  create_random_matrix(row, col, start_infected = start_inf, seed = seed)
}

#' Create a starting matrix with the four corners infected
#'
#' @param row Number of rows for the matrix.
#' @param col Number of columns for the matrix.
#' @param start_infected Accepted for consistency with
#'   [create_random_matrix()] and ignored.
#' @param seed Accepted for consistency with [create_random_matrix()] and
#'   ignored.
#'
#' @return A numeric matrix with the four corner cells encoded as infected
#'   (`1`) and all other cells encoded as susceptible (`0`).
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

#' Create a starting matrix with the center cell infected
#'
#' @param row Number of rows for the matrix.
#' @param col Number of columns for the matrix.
#' @param start_infected Accepted for consistency with
#'   [create_random_matrix()] and ignored.
#' @param seed Accepted for consistency with [create_random_matrix()] and
#'   ignored.
#'
#' @return A numeric matrix with the center cell encoded as infected (`1`) and
#'   all other cells encoded as susceptible (`0`).
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

#' Count the number cells that are infected near a specific cell
#'
#' @param x Input matrix
#' @param i Row coordinate
#' @param j Column coordinate
#'
#' @return Count
#'
#' @examples
#' x <- create_random_matrix(5, 5)
#' count_infected_neighbors(x, 2, 2)
#' @noRd
count_infected_neighbors <- function(x, i, j) {
  if (x[i, j] != SUSCEPTIBLE) {
    return(0)
  }

  neighbors <- x[(i - 1):(i + 1), (j - 1):(j + 1)]
  sum(neighbors == INFECTED) - 1L * (x[i, j] == INFECTED)
}

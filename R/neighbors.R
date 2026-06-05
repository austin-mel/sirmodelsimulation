#' Count the number cells that are infected near a specific cell
#'
#' @param x Input matrix
#' @param i Row coordinate
#' @param j Column coordinate
#'
#' @return Count
#'
#' @examples
#' x <- create_matrix(5,5)
#' count_infected_neighbors(x, 2, 2)
#' @noRd
count_infected_neighbors <- function(x, i, j) {
  count <- 0
  target <- x[i, j]

  if (target == 0) {
    # top
    if (x[i, j - 1] == 1) {
      count <- count + 1
    }
    # bot
    if (x[i, j + 1] == 1) {
      count <- count + 1
    }
    # right
    if (x[i + 1, j] == 1) {
      count <- count + 1
    }
    # left
    if (x[i - 1, j] == 1) {
      count <- count + 1
    }
    # tl
    if (x[i - 1, j - 1] == 1) {
      count <- count + 1
    }
    # tr
    if (x[i + 1, j - 1] == 1) {
      count <- count + 1
    }
    # bl
    if (x[i - 1, j + 1] == 1) {
      count <- count + 1
    }
    # br
    if (x[i + 1, j + 1] == 1) {
      count <- count + 1
    }
  }
  count
}

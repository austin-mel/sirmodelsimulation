#' Plot the matrix
#'
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics title
#' @importFrom graphics legend
#' @param input Input matrix
#'
#' @return Matrix Image
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' dis_plot(x)
#' @noRd
dis_plot <- function(input) {
  nr <- nrow(input)
  nc <- ncol(input)
  par(yaxt = "n", xaxt = "n", ann = FALSE)

  input <- apply(input, 2, rev)

  title_string <- paste("SIR Model Simulation with", nr, "x", nc, "Input Matrix")
  image(1:nc, 1:nr, t(input), zlim = c(0, 3), col = c("white", "red", "grey", "black"))
  title(title_string)
  legend("bottom", legend = c("Susceptible", "Infected", "Immune", "Desceased"), fill = c("white", "red", "grey", "black"), horiz = TRUE, inset = c(0, -.15), xpd = TRUE)
}

#' Run a given simulation multiple times to produce a heatmap to identify which cells are more likely to get infected
#'
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics title
#' @param prob_infect Probability of infection spreading
#' @param input_matrix Input matrix
#' @param runs Number of simulation runs
#' @param auto_immunity Enable auto immunity after infection
#' @param imm_prob Probability of getting immunity
#' @param allow_death Enable death
#' @param fat_prob Probability of fatality
#'
#' @return Heatmap image
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' multiple_run_heatmap(0.25,x,10)
#' multiple_run_heatmap(0.15,x,10,TRUE)
#' multiple_run_heatmap(0.35,x,10,FALSE,0.25,FALSE)
#' multiple_run_heatmap(0.25,x,10,FALSE,0.4,TRUE,0.08)
#' @export
multiple_run_heatmap <- function(prob_infect = 0.25, input_matrix, runs = 10, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15) {
  simulations <- replicate(runs, simulate_sir(prob_infect, input_matrix, allow_death = allow_death, auto_immunity = auto_immunity, imm_prob = imm_prob, fat_prob = fat_prob), simplify = FALSE)

  sx <- 0
  # Total up the final matrix values from every simulation
  for (s in simulations) {
    sx <- sx + s[["final_matrix"]]
  }

  sx <- apply(sx, 2, rev)

  # Create a heat map of most likely cells infected (higher numbers are darker colors)
  title_string <- paste("Heat map for infection simulation with", runs, "total runs")
  par(yaxt = "n", xaxt = "n", ann = FALSE)
  image(1:ncol(sx), 1:nrow(sx), t(sx))
  title(title_string)
}

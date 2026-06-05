validate_sir_inputs <- function(input_matrix, prob_infect, model = "SIR", imm_prob = 0.5) {
  if (!is.matrix(input_matrix)) {
    stop("input_matrix must be a matrix.", call. = FALSE)
  }

  if (!is.numeric(prob_infect) || length(prob_infect) != 1 || is.na(prob_infect) ||
      prob_infect < 0 || prob_infect > 1) {
    stop("prob_infect must be a number between 0 and 1.", call. = FALSE)
  }

  valid_states <- c(SUSCEPTIBLE, INFECTED, RECOVERED)
  if (any(!input_matrix %in% valid_states)) {
    stop("input_matrix contains invalid SIR state values.", call. = FALSE)
  }

  if (!model %in% c("SIR", "SIS", "SIRS")) {
    stop('model must be "SIR", "SIS", or "SIRS".', call. = FALSE)
  }

  if (!is.numeric(imm_prob) || length(imm_prob) != 1 || is.na(imm_prob) ||
      imm_prob < 0 || imm_prob > 1) {
    stop("imm_prob must be a number between 0 and 1.", call. = FALSE)
  }
}

set_optional_seed <- function(seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
}

sir_counts <- function(x, step) {
  data.frame(
    step = step,
    susceptible = sum(x == SUSCEPTIBLE),
    infected = sum(x == INFECTED),
    recovered = sum(x == RECOVERED)
  )
}

plot_sir_matrix <- function(input, main = "SIR simulation matrix") {
  state_colors <- c("white", "red", "grey", "black")
  state_labels <- c("Susceptible", "Infected", "Recovered", "Deceased")
  display_matrix <- apply(input, 2, rev)

  graphics::par(yaxt = "n", xaxt = "n", ann = FALSE)
  graphics::image(
    seq_len(ncol(input)),
    seq_len(nrow(input)),
    t(display_matrix),
    zlim = c(SUSCEPTIBLE, DECEASED),
    col = state_colors
  )
  graphics::title(main)
  graphics::legend(
    "bottom",
    legend = state_labels,
    fill = state_colors,
    horiz = TRUE,
    inset = c(0, -0.15),
    xpd = TRUE
  )
}

#' Single time step in simulation
#'
#' @param prob_infect Probability of infection spreading
#' @param input Input matrix
#' @param model Simulation model. Use "SIR", "SIS", or "SIRS".
#' @param imm_prob Probability that an infected cell becomes recovered in SIRS.
#'
#' @return new_inf_matrix
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' infect_step(0.25,x)
#' @noRd
infect_step <- function(prob_infect = 0.25, input, model = "SIR", imm_prob = 0.5) {
  validate_sir_inputs(input, prob_infect, model, imm_prob)

  nr <- nrow(input)
  nc <- ncol(input)
  input_pad <- pad(input)
  new_x <- input

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (input[i, j] == SUSCEPTIBLE) {
        infected_neighbors <- count_infected_neighbors(input_pad, i + 1, j + 1)

        if (infected_neighbors > 0) {
          infection_probability <- 1 - (1 - prob_infect)^infected_neighbors
          if (runif(1) < infection_probability) {
            new_x[i, j] <- INFECTED
          }
        }
      } else if (input[i, j] == INFECTED) {
        if (model == "SIS") {
          new_x[i, j] <- SUSCEPTIBLE
        } else if (model == "SIRS") {
          if (runif(1) < imm_prob) {
            new_x[i, j] <- RECOVERED
          } else {
            new_x[i, j] <- SUSCEPTIBLE
          }
        } else {
          new_x[i, j] <- RECOVERED
        }
      }
    }
  }

  new_x
}

#' Run simulation of infection spread until there are no longer any infected cells
#'
#' @param prob_infect Probability of infection spreading
#' @param input_matrix Input matrix
#' @param model Simulation model. Use "SIR", "SIS", or "SIRS".
#' @param imm_prob Probability that an infected cell becomes recovered in SIRS.
#' @param plot Plot each simulation step
#' @param seed Optional random seed for reproducible simulation
#'
#' @return Results List
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' simulate_sir(0.4,x)
#' @export
simulate_sir <- function(prob_infect = 0.25, input_matrix, model = "SIR", imm_prob = 0.5, plot = FALSE, seed = NULL) {
  validate_sir_inputs(input_matrix, prob_infect, model, imm_prob)
  set_optional_seed(seed)

  x <- input_matrix
  step_count <- 0
  history <- sir_counts(x, step_count)

  if (isTRUE(plot)) {
    plot_sir_matrix(x, main = paste("SIR simulation step", step_count))
  }

  while (any(x == INFECTED)) {
    x <- infect_step(prob_infect, x, model = model, imm_prob = imm_prob)
    step_count <- step_count + 1
    history <- rbind(history, sir_counts(x, step_count))

    if (isTRUE(plot)) {
      plot_sir_matrix(x, main = paste("SIR simulation step", step_count))
    }
  }

  num_recovered <- sum(x == RECOVERED)
  num_total <- ncol(x) * nrow(x)
  inf_prop <- num_recovered / num_total

  list(
    steps = step_count,
    prob_infect = prob_infect,
    inf_prop = inf_prop,
    final_matrix = x,
    history = history
  )
}

#' Run a given simulation multiple times and get the average steps and average proportion of cells infected.
#'
#' @param prob_infect Probability of infection spreading
#' @param input_matrix Input matrix
#' @param runs Number of simulation runs
#' @param model Simulation model. Use "SIR", "SIS", or "SIRS".
#' @param imm_prob Probability that an infected cell becomes recovered in SIRS.
#' @param seed Optional random seed for reproducible simulations
#'
#' @return Results Dataframe
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' simulate_many_runs(0.25,x,10)
#' @export
simulate_many_runs <- function(prob_infect = 0.25, input_matrix, runs = 10, model = "SIR", imm_prob = 0.5, seed = NULL) {
  validate_sir_inputs(input_matrix, prob_infect, model, imm_prob)
  set_optional_seed(seed)

  simulations <- replicate(
    runs,
    simulate_sir(prob_infect = prob_infect, input_matrix = input_matrix, model = model, imm_prob = imm_prob),
    simplify = FALSE
  )

  data.frame(
    total_runs = runs,
    avg_steps = mean(vapply(simulations, `[[`, numeric(1), "steps")),
    avg_prop_infected = mean(vapply(simulations, `[[`, numeric(1), "inf_prop"))
  )
}

#' Run a given simulation multiple times using a sequence of increasing infection probabilities.
#'
#' @param input_matrix Input matrix
#' @param step Value to step by in sequence of infection probabilities
#' @param model Simulation model. Use "SIR", "SIS", or "SIRS".
#' @param imm_prob Probability that an infected cell becomes recovered in SIRS.
#' @param seed Optional random seed for reproducible simulations
#'
#' @return Results Dataframe
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' simulate_inf_seq(x,0.1)
#' @export
simulate_inf_seq <- function(input_matrix, step = 0.1, model = "SIR", imm_prob = 0.5, seed = NULL) {
  validate_sir_inputs(input_matrix, 0, model, imm_prob)
  set_optional_seed(seed)

  prop_infect <- seq(from = 0.1, to = 0.999, by = step)
  results <- lapply(prop_infect, simulate_sir, input_matrix = input_matrix, model = model, imm_prob = imm_prob)

  data.frame(
    prob_infection = vapply(results, `[[`, numeric(1), "prob_infect"),
    total_steps = vapply(results, `[[`, numeric(1), "steps"),
    infected_prop = vapply(results, `[[`, numeric(1), "inf_prop")
  )
}

#' Run a given simulation multiple times to summarize which cells were previously infected
#'
#' @param prob_infect Probability of infection spreading
#' @param input_matrix Input matrix
#' @param runs Number of simulation runs
#' @param model Simulation model. Use "SIR", "SIS", or "SIRS".
#' @param imm_prob Probability that an infected cell becomes recovered in SIRS.
#' @param plot Plot the heatmap
#' @param seed Optional random seed for reproducible simulations
#'
#' @return Matrix of recovered counts by cell
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' multiple_run_heatmap(0.25,x,10)
#' @export
multiple_run_heatmap <- function(prob_infect = 0.25, input_matrix, runs = 10, model = "SIR", imm_prob = 0.5, plot = FALSE, seed = NULL) {
  validate_sir_inputs(input_matrix, prob_infect, model, imm_prob)
  set_optional_seed(seed)

  simulations <- replicate(
    runs,
    simulate_sir(prob_infect = prob_infect, input_matrix = input_matrix, model = model, imm_prob = imm_prob),
    simplify = FALSE
  )

  heatmap_matrix <- matrix(0, nrow = nrow(input_matrix), ncol = ncol(input_matrix))
  for (simulation in simulations) {
    heatmap_matrix <- heatmap_matrix + (simulation[["final_matrix"]] == RECOVERED)
  }

  if (isTRUE(plot)) {
    display_matrix <- apply(heatmap_matrix, 2, rev)
    graphics::par(yaxt = "n", xaxt = "n", ann = FALSE)
    graphics::image(
      seq_len(ncol(heatmap_matrix)),
      seq_len(nrow(heatmap_matrix)),
      t(display_matrix)
    )
    graphics::title(paste("Recovered cell count across", runs, "SIR runs"))
  }

  heatmap_matrix
}

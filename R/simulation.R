#' Single time step in simulation
#'
#' @param prob_infect Probability of infection spreading
#' @param input Input matrix
#' @param auto_immunity Enable auto immunity after infection
#' @param imm_prob Probability of getting immunity
#' @param allow_death Enable death
#' @param fat_prob Probability of fatality
#'
#' @return new_inf_matrix
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' infect_step(0.25,x)
#' infect_step(0.25,x,FALSE,0.25)
#' infect_step(0.25,x,FALSE,0.3,FALSE)
#' infect_step(0.25,x,FALSE,0.25,TRUE,0.15)
#' infect_step(0.25,x,TRUE,1.0,TRUE,0.09)
#' @noRd
infect_step <- function(prob_infect = 0.25, input, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15) {
  nr <- nrow(input)
  nc <- ncol(input)
  input_pad <- pad(input)
  new_x <- input_pad

  midrows <- 2:(nr + 1)
  midcols <- 2:(nc + 1)
  for (i in midrows) {
    for (j in midcols) {
      if (input_pad[i, j] == 0) {
        count <- count_infected_neighbors(input_pad, i, j)
        if (0 < count) {
          prob0 <- (1 - prob_infect)^count
          infect <- sample(c(TRUE, FALSE), prob = c(1 - prob0, prob0))
          if (infect[1]) {
            new_x[i, j] <- 1
          }
        }
      }
    }
  }

  # If auto immunity is disabled & death is disabled
  # Allow infected cells to become Susceptible or Immune after infection
  if (identical(FALSE, auto_immunity) && identical(FALSE, allow_death)) {
    new_x[input_pad == 1] <- 3
    for (y in midrows) {
      for (z in midcols) {
        if (new_x[y, z] == 3) {
          result <- sample(0:1, 1, prob = c(1 - imm_prob, imm_prob))
          if (result == 1) {
            result <- 2
          }
          new_x[y, z] <- result
        }
      }
    }
  }
  # If auto immunity is disabled & death is enabled
  # Allow infected cells to become Susceptible, Immune, or Deceased after infection
  else if (identical(FALSE, auto_immunity) && identical(TRUE, allow_death)) {
    new_x[input_pad == 1] <- 4

    for (y in midrows) {
      for (z in midcols) {
        if (new_x[y, z] == 4) {
          result <- sample(3:4, 1, prob = c(fat_prob, 1 - fat_prob))
          if (result == 4) {
            result <- sample(1:2, 1, prob = c(imm_prob, 1 - imm_prob))
            if (result == 1) {
              result <- 0
            }
          }
          new_x[y, z] <- result
        }
      }
    }
  }
  # If auto immunity is enabled & death is enabled
  # Allow infected cells to become Immune or Deceased after infection
  else if (identical(TRUE, auto_immunity) && identical(TRUE, allow_death)) {
    new_x[input_pad == 1] <- 4

    for (y in midrows) {
      for (z in midcols) {
        if (new_x[y, z] == 4) {
          result <- sample(3:4, 1, prob = c(fat_prob, 1 - fat_prob))
          if (result == 4) {
            result <- 2
          }
          new_x[y, z] <- result
        }
      }
    }
  }
  # Else auto immunity is enabled & death is disabled
  else {
    new_x[input_pad == 1] <- 2
  }

  completed_step <- new_x[midrows, midcols]

  # Visualize matrix
  dis_plot(completed_step)

  # Return matrix
  return(completed_step)
}

#' Run simulation of infection spread until there are no longer any infected cells
#'
#' @param prob_infect Probability of infection spreading
#' @param input_matrix Input matrix
#' @param auto_immunity Enable auto immunity after infection
#' @param imm_prob Probability of getting immunity
#' @param allow_death Enable death
#' @param fat_prob Probability of fatality
#'
#' @return Results List
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' simulate_sir(0.4,x)
#' simulate_sir(0.15,x,TRUE)
#' simulate_sir(0.35,x,FALSE,0.5,FALSE)
#' simulate_sir(0.55,x,FALSE,0.35,TRUE,0.1)
#' @export
simulate_sir <- function(prob_infect = 0.25, input_matrix, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15) {
  x <- input_matrix

  step_count <- 0

  while (any(x == 1)) {
    # If allow_death == F require only immunity probability argument
    if (identical(FALSE, allow_death)) {
      x <- infect_step(prob_infect, x, auto_immunity, imm_prob)
    }
    # If allow_death == T require immunity probability and fatality probability argument
    else {
      x <- infect_step(prob_infect, x, auto_immunity, imm_prob, allow_death, fat_prob)
    }

    # Increment step count for each run
    step_count <- step_count + 1

    Sys.sleep(0.1)
  }

  num_inf <- sum(x == 1) + sum(x == 2)
  num_healthy <- sum(x == 0)
  num_total <- ncol(x) * nrow(x)

  inf_prop <- num_inf / num_total

  # Print list of results
  if (identical(FALSE, auto_immunity) && identical(FALSE, allow_death)) {
    list(steps = step_count, prob_infect = prob_infect, inf_prop = inf_prop, imm_prob = imm_prob, final_matrix = x)
  }
  else if (identical(FALSE, auto_immunity) && identical(TRUE, allow_death)) {
    list(steps = step_count, prob_infect = prob_infect, inf_prop = inf_prop, imm_prob = imm_prob, fat_prob = fat_prob, final_matrix = x)
  }
  else if (identical(TRUE, auto_immunity) && identical(TRUE, allow_death)) {
    list(steps = step_count, prob_infect = prob_infect, inf_prop = inf_prop, fat_prob = fat_prob, final_matrix = x)
  }
  else {
    list(steps = step_count, prob_infect = prob_infect, inf_prop = inf_prop, final_matrix = x)
  }
}

#' Run a given simulation multiple times using a sequence of increasing infection probabilities.
#'
#' @param input_matrix Input matrix
#' @param step Value to step by in sequence of infection probabilities
#' @param auto_immunity Enable auto immunity after infection
#' @param imm_prob Probability of getting immunity
#' @param allow_death Enable death
#' @param fat_prob Probability of fatality
#'
#' @return Results Dataframe
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' simulate_inf_seq(x,0.1,TRUE)
#' simulate_inf_seq(x,0.2,FALSE,0.35,FALSE)
#' simulate_inf_seq(x,0.1,FALSE,0.35,TRUE,0.12)
#' @export
simulate_inf_seq <- function(input_matrix, step = 0.1, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15) {
  prop_infect <- seq(from = 0.1, to = 0.999, by = step)

  # Use each probability in defined sequence in a given simulation run
  # Save results via list
  results <- lapply(prop_infect, simulate_sir, input_matrix = input_matrix, auto_immunity = auto_immunity, imm_prob = imm_prob, allow_death = allow_death, fat_prob = fat_prob)

  # Save results in data frame
  if (identical(FALSE, auto_immunity) && identical(FALSE, allow_death)) {
    results2 <- data.frame(
      prob_infection = sapply(results, `[[`, "prob_infect"),
      total_steps = sapply(results, `[[`, "steps"),
      infected_prop = sapply(results, `[[`, "inf_prop"),
      immunity_prob = sapply(results, `[[`, "imm_prob")
    )
  }
  else if (identical(FALSE, auto_immunity) && identical(TRUE, allow_death)) {
    results2 <- data.frame(
      prob_infection = sapply(results, `[[`, "prob_infect"),
      total_steps = sapply(results, `[[`, "steps"),
      infected_prop = sapply(results, `[[`, "inf_prop"),
      immunity_prob = sapply(results, `[[`, "imm_prob"),
      fatality_prob = sapply(results, `[[`, "fat_prob")
    )
  }
  else if (identical(TRUE, auto_immunity) && identical(TRUE, allow_death)) {
    results2 <- data.frame(
      prob_infection = sapply(results, `[[`, "prob_infect"),
      total_steps = sapply(results, `[[`, "steps"),
      infected_prop = sapply(results, `[[`, "inf_prop"),
      fatality_prob = sapply(results, `[[`, "fat_prob")
    )
  }
  else {
    results2 <- data.frame(
      prob_infection = sapply(results, `[[`, "prob_infect"),
      total_steps = sapply(results, `[[`, "steps"),
      infected_prop = sapply(results, `[[`, "inf_prop")
    )
  }

  print(results2)
}

#' Run a given simulation multiple times and get the average steps and average proportion of cells infected.
#'
#' @param prob_infect Probability of infection spreading
#' @param input_matrix Input matrix
#' @param runs Number of simulation runs
#' @param auto_immunity Enable auto immunity after infection
#' @param imm_prob Probability of getting immunity
#' @param allow_death Enable death
#' @param fat_prob Probability of fatality
#'
#' @return Results Dataframe
#'
#' @examples
#' x <- create_matrix(5,5)
#'
#' simulate_many_runs(0.25,x,10)
#' simulate_many_runs(0.15,x,10,TRUE)
#' simulate_many_runs(0.35,x,10,FALSE,0.25,FALSE)
#' simulate_many_runs(0.25,x,10,FALSE,0.4,TRUE,0.08)
#' @export
simulate_many_runs <- function(prob_infect = 0.25, input_matrix, runs = 10, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15) {
  simulations <- replicate(runs, simulate_sir(prob_infect, input_matrix, allow_death = allow_death, auto_immunity = auto_immunity, imm_prob = imm_prob, fat_prob = fat_prob), simplify = FALSE)

  sx <- 0
  sp <- 0

  # Total up the steps and proportion of infected cells from every simulation
  for (s in simulations) {
    sx <- sx + s[["steps"]]
    sp <- sp + s[["inf_prop"]]
  }

  steps_avg <- sx / runs
  prop_avg <- sp / runs

  # Save results in data frame
  results <- data.frame(
    total_runs = runs,
    avg_steps = steps_avg,
    avg_prop_infected = prop_avg
  )

  print(results)
}

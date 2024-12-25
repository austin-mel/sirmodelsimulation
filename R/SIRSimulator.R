# SIR Model Simulation - Austin Melendez

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
create_matrix <- function(row,col,start_inf=0.125){
  set.seed(94128)
  inf_matrix <- matrix(sample(c(0,1), row*col, replace = TRUE, prob=c(1-start_inf,start_inf)), nrow=row, ncol=col)
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
create_crnr_matrix <- function(row,col){
  inf_matrix <- matrix(0, row, col)
  inf_matrix[1,1] <- 1
  inf_matrix[1,col] <- 1
  inf_matrix[row,1] <- 1
  inf_matrix[row,col] <- 1
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
create_cntr_matrix <- function(row,col){
  inf_matrix <- matrix(0, nrow=row, ncol=col)
  inf_matrix[(row/2),(col/2)] <- 1
  inf_matrix
}

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
dis_plot <- function(input){
  nr <- nrow(input)
  nc <- ncol(input)
  par(yaxt="n", xaxt="n",ann=FALSE)

  input <- apply(input, 2, rev)

  title_string <- paste("SIR Model Simulation with",nr,"x",nc,"Input Matrix" )
  image(1:nc,1:nr,t(input), zlim = c(0,3), col = c("white","red","grey","black"))
  title(title_string)
  legend("bottom", legend=c("Susceptible","Infected","Immune","Desceased"), fill=c("white","red","grey","black"), horiz=TRUE, inset=c(0,-.15), xpd=TRUE)
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
pad <- function(input){
  nr <- nrow(input)
  nc <- ncol(input)
  input_pad <- matrix(0, nrow = nr+2, ncol = nc+2)
  input_pad[2:(nr + 1), 2:(nc + 1)] <- input
  input_pad
}

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
count_infected_neighbors <- function(x, i, j){
  count <- 0
  target <- x[i,j]

  if(target == 0){
    #top
    if(x[i,j-1] == 1){
      count <- count + 1
    }
    #bot
    if(x[i,j+1] == 1){
      count <- count + 1
    }
    #right
    if(x[i+1,j] == 1){
      count <- count + 1
    }
    #left
    if(x[i-1,j] == 1){
      count <- count + 1
    }
    #tl
    if(x[i-1,j-1] == 1){
      count <- count + 1
    }
    #tr
    if(x[i+1,j-1] == 1){
      count <- count + 1
    }
    #bl
    if(x[i-1,j+1] == 1){
      count <- count + 1
    }
    #br
    if(x[i+1,j+1] == 1){
      count <- count + 1
    }
  }
  count
}

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
infect_step <- function(prob_infect = 0.25, input, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15){
  nr <- nrow(input)
  nc <- ncol(input)
  input_pad <- pad(input)
  new_x <- input_pad

  midrows <- 2:(nr + 1)
  midcols <- 2:(nc + 1)
  for(i in midrows){
    for(j in midcols){
      if(input_pad[i,j] == 0){ # Cell is susceptible
        count <- count_infected_neighbors(input_pad, i, j)
        if(0 < count){ # Has infected neighbors
          prob0 <- (1 - prob_infect)^count
          infect <- sample(c(TRUE,FALSE), prob = c(1 - prob0, prob0))
          if(infect[1]){ # Cell gets infected
            new_x[i,j] <- 1
          }
        }
      }
    }
  }

  # If auto immunity is disabled & death is disabled
  # Allow infected cells to become Susceptible or Immune after infection
  if(identical(FALSE,auto_immunity) && identical(FALSE, allow_death)){
    new_x[input_pad == 1] <- 3
    for(y in midrows){
      for(z in midcols){
        if(new_x[y,z] == 3){
          result <- sample(0:1, 1, prob = c(1 - imm_prob, imm_prob))
          if(result == 1){
            result <- 2
          }
          new_x[y,z] <- result
        }
      }
    }
  }
  # If auto immunity is disabled & death is enabled
  # Allow infected cells to become Susceptible, Immune, or Deceased after infection
  else if(identical(FALSE,auto_immunity) && identical(TRUE,allow_death)){
    new_x[input_pad == 1] <- 4

    for(y in midrows){
      for(z in midcols){
        if(new_x[y,z] == 4){
          result <- sample(3:4, 1, prob = c(fat_prob, 1 - fat_prob))
          if(result == 4){
            result <- sample(1:2, 1, prob = c(imm_prob, 1 - imm_prob))
            if(result == 1){
              result <- 0
            }
          }
          new_x[y,z] <- result
        }
      }
    }
  }
  # If auto immunity is enabled & death is enabled
  # Allow infected cells to become Immune or Deceased after infection
  else if(identical(TRUE,auto_immunity) && identical(TRUE,allow_death)){
    new_x[input_pad == 1] <- 4

    for(y in midrows){
      for(z in midcols){
        if(new_x[y,z] == 4){
          result <- sample(3:4, 1, prob = c(fat_prob, 1 - fat_prob))
          if(result == 4){
            result <- 2
          }
          new_x[y,z] <- result
        }
      }
    }
  }
  # Else auto immunity is enabled & death is disabled
  else{
    new_x[input_pad == 1] <- 2
  }

  completed_step <- new_x[midrows,midcols]

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
simulate_sir <- function(prob_infect = 0.25, input_matrix, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15){
  x <- input_matrix

  step_count <- 0

  while(any(x == 1)){
    # If allow_death == F require only immunity probability argument
    if(identical(FALSE,allow_death)){
      x <- infect_step(prob_infect, x, auto_immunity, imm_prob)
    }
    # If allow_death == T require immunity probability and fatality probability argument
    else{
      x <- infect_step(prob_infect, x, auto_immunity, imm_prob, allow_death, fat_prob)
    }

    #Increment step count for each run
    step_count <- step_count + 1

    Sys.sleep(0.1)
  }

  num_inf <- sum(x==1) + sum(x==2)
  num_healthy <- sum(x==0)
  num_total <- ncol(x) * nrow(x)

  inf_prop <- num_inf / num_total

  #Print list of results
  if(identical(FALSE,auto_immunity) && identical(FALSE,allow_death)){
    list(steps = step_count, prob_infect = prob_infect, inf_prop = inf_prop, imm_prob = imm_prob, final_matrix = x)
  }
  else if(identical(FALSE,auto_immunity) && identical(TRUE,allow_death)){
    list(steps = step_count, prob_infect = prob_infect, inf_prop = inf_prop, imm_prob = imm_prob, fat_prob = fat_prob, final_matrix = x)
  }
  else if(identical(TRUE,auto_immunity) && identical(TRUE,allow_death)){
    list(steps = step_count, prob_infect = prob_infect, inf_prop = inf_prop, fat_prob = fat_prob, final_matrix = x)
  }
  else{
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
simulate_inf_seq <- function(input_matrix, step = 0.1, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15){
  prop_infect <- seq(from = 0.1, to = 0.999, by = step)

  # Use each probability in defined sequence in a given simulation run
  # Save results via list
  results <- lapply(prop_infect, simulate_sir, input_matrix = input_matrix, auto_immunity = auto_immunity, imm_prob = imm_prob, allow_death = allow_death, fat_prob = fat_prob)

  # Save results in data frame
  if(identical(FALSE,auto_immunity) && identical(FALSE,allow_death)){
    results2 <- data.frame(
      prob_infection = sapply(results, `[[`, "prob_infect"),
      total_steps = sapply(results, `[[`, "steps"),
      infected_prop = sapply(results, `[[`, "inf_prop"),
      immunity_prob = sapply(results, `[[`, "imm_prob")
    )
  }
  else if(identical(FALSE,auto_immunity) && identical(TRUE,allow_death)){
    results2 <- data.frame(
      prob_infection = sapply(results, `[[`, "prob_infect"),
      total_steps = sapply(results, `[[`, "steps"),
      infected_prop = sapply(results, `[[`, "inf_prop"),
      immunity_prob = sapply(results, `[[`, "imm_prob"),
      fatality_prob = sapply(results, `[[`, "fat_prob")
    )
  }
  else if(identical(TRUE,auto_immunity) && identical(TRUE,allow_death)){
    results2 <- data.frame(
      prob_infection = sapply(results, `[[`, "prob_infect"),
      total_steps = sapply(results, `[[`, "steps"),
      infected_prop = sapply(results, `[[`, "inf_prop"),
      fatality_prob = sapply(results, `[[`, "fat_prob")
    )
  }
  else{
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
simulate_many_runs <- function(prob_infect = 0.25, input_matrix, runs = 10, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15){
  simulations <- replicate(runs, simulate_sir(prob_infect, input_matrix, allow_death = allow_death, auto_immunity = auto_immunity, imm_prob = imm_prob, fat_prob = fat_prob), simplify = FALSE)

  sx <- 0
  sp <- 0

  # Total up the steps and proportion of infected cells from every simulation
  for(s in simulations){
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
multiple_run_heatmap <- function(prob_infect = 0.25, input_matrix, runs = 10, auto_immunity = FALSE, imm_prob = 0.5, allow_death = FALSE, fat_prob = 0.15){
  simulations <- replicate(runs, simulate_sir(prob_infect, input_matrix, allow_death = allow_death, auto_immunity = auto_immunity, imm_prob = imm_prob, fat_prob = fat_prob), simplify = FALSE)

  sx <- 0
  # Total up the final matrix values from every simulation
  for(s in simulations){
    sx <- sx + s[["final_matrix"]]
  }

  sx <- apply(sx, 2, rev)


  # Create a heat map of most likely cells infected (higher numbers are darker colors)
  title_string <- paste("Heat map for infection simulation with",runs,"total runs")
  par(yaxt="n", xaxt="n",ann=FALSE)
  image(1:ncol(sx),1:nrow(sx),t(sx))
  title(title_string)
}

#' Create corner starting matrix
#'
#' @param row Number of rows in matrix
#' @param col Number of columns in matrix
#'
#' @return SIRmatrix
#' @export
#'
#' @examples
#' create_crnr_matrix(20, 30)
create_crnr_matrix <- function(row=20,col=30){
  inf_matrix <- matrix(0, nrow=row, ncol=col)
  inf_matrix[1,1] <- 1
  inf_matrix[1,col] <- 1
  inf_matrix[row,1] <- 1
  inf_matrix[row,col] <- 1
  print(inf_matrix)
}

#' Create infection matrix
#'
#' @param row Number of rows in matrix
#' @param col Number of columns in matrix
#' @param start_prop Probability of cells initially infected
#'
#' @return SIR Matrix
#' @export
#'
#' @examples
#' init_matrix(15, 15, 0.125)
init_matrix <- function(row,col,start_prop=0.125){
  inf_matrix <- matrix(sample(c(0,1), row*col, replace = TRUE, prob=c(1-start_prop,start_prop)), nrow=row, ncol=col)
}

#' Single step in simulation
#'
#' @param input Input SIR Matrix
#' @param inf_prob Probability of spreading infection
#'
#' @return Updated SIR Matrix after infection spread during single step
#' @export
#'
#' @examples
#' infection_step(, 0.5)
infection_step <- function(input = init_matrix(5,5), inf_prob = 0.5){
  x <- input

  cc_RL <- init_matrix(nrow(x),ncol(x)-1,inf_prob)
  cc_TB <- init_matrix(nrow(x)-1,ncol(x),inf_prob)
  cc_DIAG <- init_matrix(nrow(x)-1,ncol(x)-1,inf_prob)

  max_col <- ncol(x)
  max_row <- nrow(x)

  xnew <- x

  #ABOVE
  new_inf_above <- ((x[-1, ] == 1) & cc_TB)
  xnew[-max_row, ][new_inf_above] <- 1

  xnew <- x

  #BELOW
  new_inf_below <- ((x[-max_row, ] == 1) & cc_TB)
  xnew[-1, ][new_inf_below] <- 1

  xnew <- x

  #RIGHT
  new_inf_right <- ((x[ ,-max_col] == 1) & cc_RL)
  xnew[, -1][new_inf_right] <- 1

  xnew <- x

  #LEFT
  new_inf_left <- ((x[ ,-1] == 1) & cc_RL)
  xnew[, -max_col][new_inf_left] <- 1

  xnew <- x

  #TOP RIGHT
  new_inf_top_right <- ((x[-1,-max_col] == 1) & cc_DIAG)
  xnew[-max_row, -1][new_inf_top_right] <- 1

  xnew <- x

  #TOP LEFT
  new_inf_top_left <- ((x[-1,-1] == 1) & cc_DIAG)
  xnew[-1, -max_col][new_inf_top_left] <- 1

  xnew <- x

  #BOTTOM LEFT
  new_inf_bot_left <- ((x[-max_row,-1] == 1) & cc_DIAG)
  xnew[-1,-max_col][new_inf_bot_left] <- 1

  xnew <- x

  #BOTTOM RIGHT
  new_inf_bot_right <- ((x[-max_row,-max_col] == 1) & cc_DIAG)
  xnew[-1,-1][new_inf_bot_right] <- 1

  #Set previously infected cells to 2
  xnew[x == 1] <- 2

  #print(xnew)
  return(xnew)
}

#' Run SIR simulation
#'
#' @param inf_prob Probability of infection spreading
#' @param row Number of rows in matrix
#' @param col Number of columns in matrix
#' @param start_inf Probability of cells initially infected
#'
#' @return List of total steps taken, starting infection probability, infection spread probability, and final SIR matrix result.
#' @export
#'
#' @examples
#' simulate_sir(0.5, 20, 30, 0.125)
simulate_sir <- function(inf_prob = 0.5, row, col, start_inf = 0.125){
  x <- matrix(sample(c(0,1), row*col, replace = TRUE, prob=c(1-start_inf,start_inf)), nrow=row, ncol=col)

  step_count <- 0

  while(any(x == 1)){
    x <- infection_step(x, inf_prob)
    #Increment step count for each run
    step_count <- step_count + 1
  }

  num_inf <- sum(x==1) + sum(x==2)
  num_healthy <- sum(x==0)
  num_total <- ncol(x) * nrow(x)

  prop <- num_inf / num_total

  #Print list of results
  list(steps = step_count, start_inf = start_inf, inf_prob = inf_prob, prop = prop, final_matrix = x)
}

#' Run simulation with different infection probabilities
#'
#' @param prob Probability of spreading infection
#' @param row Number of rows in matrix
#' @param col Number of columns in matrix
#' @param step Value to increment probability range by
#'
#' @return Dataframe with statistics about all simulations run
#' @export
#'
#' @examples
#' simulation(0.5, 1000, 20, 30)
simulation <- function(prob, step = 0.1, row, col){
  prop_infect <- seq(from = 0.1, to = 0.9, by = step)

  results <- lapply(prop_infect, simulate_sir, row=row, col=row)

  results2 <- data.frame(
    prob = sapply(results, `[[`, "inf_prob"),
    steps = sapply(results, `[[`, "steps"),
    prop = sapply(results, `[[`, "prop")
  )

  print(results2)
}


#' Run multiple simulations
#'
#' @importFrom graphics par
#' @importFrom graphics image
#' @param prob Probability of spreading infection to other cells
#' @param runs Number of simulation runs
#' @param row Number of rows in matrix
#' @param col Number of columns in matrix
#'
#' @return List with each simulation statistics and final infection matrix
#' @export
#'
#' @examples
#' multiple_sims(0.5, 1000, 15, 15)
multiple_sims <- function(prob, runs, row, col){
  simulations <- replicate(runs, simulate_sir(prob,row,col), simplify = FALSE)

  sx <- 0
  for(s in simulations){
    sx <- sx + s[["final_matrix"]]
  }

  sx <- apply(sx, 2, rev)

  par(yaxt="n", xaxt="n",ann=FALSE)
  image(1:ncol(sx),1:nrow(sx),t(sx))
}

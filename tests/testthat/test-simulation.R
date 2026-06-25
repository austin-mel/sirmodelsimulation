test_that("infect_step infects susceptible cells and recovers infected cells in SIR", {
  input <- matrix(
    c(
      0, 0, 0,
      0, 1, 0,
      0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )

  result <- SIRSsim:::infect_step(prob_infect = 1, input = input, model = "SIR")

  expect_equal(result[2, 2], 2)
  expect_true(all(result[cbind(c(1, 1, 1, 2, 2, 3, 3, 3), c(1, 2, 3, 1, 3, 1, 2, 3))] == 1))
})

test_that("simulate_sir returns summary fields and per-step history", {
  input <- matrix(c(1, 0, 0, 0), nrow = 2)

  result <- simulate_sir(prob_infect = 0, input_matrix = input)

  expect_named(result, c("steps", "prob_infect", "inf_prop", "final_matrix", "history"))
  expect_equal(result$steps, 1)
  expect_equal(result$prob_infect, 0)
  expect_equal(result$inf_prop, 0.25)
  expect_equal(result$final_matrix, matrix(c(2, 0, 0, 0), nrow = 2))
  expect_equal(result$history$infected, c(1, 0))
})

test_that("simulate_sir omits full_log by default", {
  input <- matrix(c(1, 0, 0, 0), nrow = 2)

  result <- simulate_sir(prob_infect = 0, input_matrix = input)

  expect_false("full_log" %in% names(result))
})

test_that("simulate_sir includes per-cell per-step full_log when requested", {
  input <- matrix(c(1, 0, 0, 0), nrow = 2)

  result <- simulate_sir(
    prob_infect = 0,
    input_matrix = input,
    model = "SIR",
    full_log = TRUE
  )

  expect_true("full_log" %in% names(result))
  expect_s3_class(result$full_log, "data.frame")
  expect_equal(nrow(result$full_log), length(input) * (result$steps + 1))
  expect_equal(
    names(result$full_log),
    c(
      "step",
      "cell_id",
      "row",
      "col",
      "state",
      "model",
      "prob_infect",
      "auto_immunity",
      "imm_prob",
      "allow_death",
      "fat_prob",
      "was_susceptible",
      "was_infected",
      "was_immune",
      "was_deceased"
    )
  )
  expect_equal(sort(unique(result$full_log$step)), 0:result$steps)
  expect_equal(result$full_log$state[result$full_log$step == 0 & result$full_log$cell_id == 1], 1)
  expect_equal(result$full_log$state[result$full_log$step == 1 & result$full_log$cell_id == 1], 2)
})

test_that("full_log contains enough information to derive survival endpoints", {
  input <- matrix(c(1, 0, 0, 0), nrow = 2)

  result <- simulate_sir(prob_infect = 0, input_matrix = input, full_log = TRUE)
  cell_one <- result$full_log[result$full_log$cell_id == 1, ]

  first_infected_step <- min(cell_one$step[cell_one$was_infected])
  recovered_step <- min(cell_one$step[cell_one$was_immune])
  death_step <- if (any(cell_one$was_deceased)) min(cell_one$step[cell_one$was_deceased]) else NA
  last_observed_step <- max(cell_one$step)
  final_state <- cell_one$state[cell_one$step == last_observed_step]

  expect_equal(first_infected_step, 0)
  expect_equal(recovered_step, 1)
  expect_true(is.na(death_step))
  expect_equal(last_observed_step, result$steps)
  expect_equal(final_state, 2)
  expect_true(any(cell_one$was_infected))
  expect_true(any(cell_one$was_immune))
  expect_false(any(cell_one$was_deceased))
})

test_that("simulate_many_runs returns averaged SIR results", {
  input <- matrix(c(1, 0, 0, 0), nrow = 2)

  result <- simulate_many_runs(prob_infect = 0, input_matrix = input, runs = 3)

  expect_s3_class(result, "data.frame")
  expect_equal(result$total_runs, 3)
  expect_equal(result$avg_steps, 1)
  expect_equal(result$avg_prop_infected, 0.25)
})

test_that("simulate_inf_seq returns one row per probability", {
  input <- matrix(1, nrow = 1, ncol = 1)

  result <- simulate_inf_seq(input_matrix = input, step = 0.45)

  expect_s3_class(result, "data.frame")
  expect_equal(result$prob_infection, c(0.1, 0.55))
  expect_equal(result$total_steps, c(1, 1))
  expect_equal(result$infected_prop, c(1, 1))
})

test_that("multiple_run_heatmap returns ever-infected counts by cell", {
  input <- matrix(c(1, 0, 0, 0), nrow = 2)

  result <- multiple_run_heatmap(prob_infect = 0, input_matrix = input, runs = 3)

  expect_equal(result, matrix(c(3, 0, 0, 0), nrow = 2))
})

test_that("multiple_run_heatmap counts infected cells that do not end recovered", {
  input <- matrix(1, nrow = 1, ncol = 1)

  result <- multiple_run_heatmap(
    prob_infect = 0,
    input_matrix = input,
    runs = 3,
    model = "SIS"
  )

  expect_equal(result, matrix(3, nrow = 1, ncol = 1))
})

test_that("model SIR recovers infected cells", {
  input <- matrix(1, nrow = 1, ncol = 1)

  result <- SIRSsim:::infect_step(prob_infect = 0, input = input, model = "SIR")

  expect_equal(result[1, 1], 2)
})

test_that("model SIS returns infected cells to susceptible", {
  input <- matrix(1, nrow = 1, ncol = 1)

  result <- SIRSsim:::infect_step(prob_infect = 0, input = input, model = "SIS")

  expect_equal(result[1, 1], 0)
})

test_that("model SIRS uses imm_prob to choose recovered or susceptible", {
  input <- matrix(1, nrow = 1, ncol = 1)

  recovered <- SIRSsim:::infect_step(prob_infect = 0, input = input, model = "SIRS", imm_prob = 1)
  susceptible <- SIRSsim:::infect_step(prob_infect = 0, input = input, model = "SIRS", imm_prob = 0)

  expect_equal(recovered[1, 1], 2)
  expect_equal(susceptible[1, 1], 0)
})

test_that("allow_death and fat_prob control mortality edge cases", {
  input <- matrix(1, nrow = 1, ncol = 1)

  deceased <- SIRSsim:::infect_step(
    prob_infect = 0,
    input = input,
    model = "SIR",
    allow_death = TRUE,
    fat_prob = 1
  )
  recovered <- SIRSsim:::infect_step(
    prob_infect = 0,
    input = input,
    model = "SIR",
    allow_death = TRUE,
    fat_prob = 0
  )

  expect_equal(deceased[1, 1], 3)
  expect_equal(recovered[1, 1], 2)
})

test_that("invalid probability inputs are rejected", {
  input <- matrix(1, nrow = 1, ncol = 1)

  expect_error(simulate_sir(prob_infect = -0.1, input_matrix = input), "prob_infect")
  expect_error(simulate_sir(prob_infect = 1.1, input_matrix = input), "prob_infect")
  expect_error(simulate_sir(prob_infect = 0.5, input_matrix = input, imm_prob = -0.1), "imm_prob")
  expect_error(simulate_sir(prob_infect = 0.5, input_matrix = input, fat_prob = 1.1), "fat_prob")
})

test_that("invalid matrix inputs are rejected", {
  expect_error(simulate_sir(prob_infect = 0.5, input_matrix = c(1, 0)), "input_matrix")
  expect_error(simulate_sir(prob_infect = 0.5, input_matrix = matrix(4, nrow = 1)), "invalid SIR state")
})

test_that("invalid model names are rejected", {
  input <- matrix(1, nrow = 1, ncol = 1)

  expect_error(simulate_sir(prob_infect = 0.5, input_matrix = input, model = "SEIR"), "model")
})

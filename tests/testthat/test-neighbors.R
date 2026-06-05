test_that("count_infected_neighbors counts all eight neighboring cells", {
  input <- matrix(
    c(
      1, 1, 1,
      1, 0, 1,
      1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  padded <- SIRSimulation:::pad(input)

  expect_equal(SIRSimulation:::count_infected_neighbors(padded, 3, 3), 8)
})

test_that("count_infected_neighbors returns zero for non-susceptible targets", {
  input <- matrix(
    c(
      1, 1, 1,
      1, 1, 1,
      1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  padded <- SIRSimulation:::pad(input)

  expect_equal(SIRSimulation:::count_infected_neighbors(padded, 3, 3), 0)
})

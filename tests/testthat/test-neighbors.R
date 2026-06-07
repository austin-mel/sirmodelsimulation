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
  padded <- SIRSsim:::pad(input)

  expect_equal(SIRSsim:::count_infected_neighbors(padded, 3, 3), 8)
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
  padded <- SIRSsim:::pad(input)

  expect_equal(SIRSsim:::count_infected_neighbors(padded, 3, 3), 0)
})

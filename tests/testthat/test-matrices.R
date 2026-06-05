test_that("create_matrix returns a reproducible susceptible/infected matrix", {
  first <- create_matrix(row = 3, col = 4, start_inf = 0.5, seed = 123)
  second <- create_matrix(row = 3, col = 4, start_inf = 0.5, seed = 123)

  expect_equal(dim(first), c(3, 4))
  expect_true(all(first %in% c(0, 1)))
  expect_equal(first, second)
})

test_that("create_crnr_matrix infects only the four corners", {
  result <- create_crnr_matrix(row = 4, col = 5)
  expected <- matrix(0, nrow = 4, ncol = 5)
  expected[1, 1] <- 1
  expected[1, 5] <- 1
  expected[4, 1] <- 1
  expected[4, 5] <- 1

  expect_equal(result, expected)
})

test_that("create_cntr_matrix infects the center cell", {
  result <- create_cntr_matrix(row = 5, col = 5)
  expected <- matrix(0, nrow = 5, ncol = 5)
  expected[3, 3] <- 1

  expect_equal(result, expected)
})

test_that("pad surrounds a matrix with susceptible cells", {
  input <- matrix(c(1, 0, 0, 1), nrow = 2)
  result <- SIRSimulation:::pad(input)
  expected <- matrix(0, nrow = 4, ncol = 4)
  expected[2:3, 2:3] <- input

  expect_equal(result, expected)
})

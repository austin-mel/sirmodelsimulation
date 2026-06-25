test_that("create_random_matrix returns exact-count matrices reproducibly", {
  first <- create_random_matrix(row = 3, col = 4, start_infected = 5, seed = 123)
  second <- create_random_matrix(row = 3, col = 4, start_infected = 5, seed = 123)

  expect_equal(dim(first), c(3, 4))
  expect_true(all(first %in% c(0, 1)))
  expect_equal(sum(first == 1), 5)
  expect_equal(first, second)
})

test_that("create_random_matrix supports exact infected cell counts", {
  result <- create_random_matrix(row = 3, col = 4, start_infected = 5, seed = 123)

  expect_equal(dim(result), c(3, 4))
  expect_equal(sum(result == 1), 5)
  expect_equal(sum(result == 0), 7)
})

test_that("create_random_matrix supports zero infected cells", {
  result <- create_random_matrix(row = 3, col = 4, start_infected = 0, seed = 123)

  expect_equal(result, matrix(0, nrow = 3, ncol = 4))
})

test_that("create_random_matrix rejects invalid start_infected values", {
  expect_error(create_random_matrix(3, 4, start_infected = NA), "start_infected")
  expect_error(create_random_matrix(3, 4, start_infected = -1), "start_infected")
  expect_error(create_random_matrix(3, 4, start_infected = 13), "start_infected")
  expect_error(create_random_matrix(3, 4, start_infected = 1.5), "whole number")
})

test_that("create_matrix remains a reproducible exact-count wrapper", {
  first <- create_matrix(row = 3, col = 4, start_inf = 5, seed = 123)
  second <- create_matrix(row = 3, col = 4, start_inf = 5, seed = 123)

  expect_equal(dim(first), c(3, 4))
  expect_true(all(first %in% c(0, 1)))
  expect_equal(sum(first == 1), 5)
  expect_equal(first, second)
})

test_that("create_corner_matrix infects only the four corners", {
  result <- create_corner_matrix(row = 4, col = 5)
  expected <- matrix(0, nrow = 4, ncol = 5)
  expected[1, 1] <- 1
  expected[1, 5] <- 1
  expected[4, 1] <- 1
  expected[4, 5] <- 1

  expect_equal(result, expected)
})

test_that("create_corner_matrix accepts shared constructor arguments without changing pattern", {
  result <- create_corner_matrix(row = 4, col = 5, start_infected = 2, seed = 123)
  expected <- matrix(0, nrow = 4, ncol = 5)
  expected[1, 1] <- 1
  expected[1, 5] <- 1
  expected[4, 1] <- 1
  expected[4, 5] <- 1

  expect_equal(result, expected)
})

test_that("create_center_matrix infects the center cell", {
  result <- create_center_matrix(row = 5, col = 5)
  expected <- matrix(0, nrow = 5, ncol = 5)
  expected[3, 3] <- 1

  expect_equal(result, expected)
})

test_that("create_center_matrix accepts shared constructor arguments without changing pattern", {
  result <- create_center_matrix(row = 5, col = 5, start_infected = 4, seed = 123)
  expected <- matrix(0, nrow = 5, ncol = 5)
  expected[3, 3] <- 1

  expect_equal(result, expected)
})

test_that("pad surrounds a matrix with susceptible cells", {
  input <- matrix(c(1, 0, 0, 1), nrow = 2)
  result <- SIRSsim:::pad(input)
  expected <- matrix(0, nrow = 4, ncol = 4)
  expected[2:3, 2:3] <- input

  expect_equal(result, expected)
})

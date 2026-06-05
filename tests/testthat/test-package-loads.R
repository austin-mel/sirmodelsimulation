test_that("SIRSimulation exports the existing public functions", {
  expected_exports <- c(
    "create_cntr_matrix",
    "create_crnr_matrix",
    "create_matrix",
    "multiple_run_heatmap",
    "simulate_inf_seq",
    "simulate_many_runs",
    "simulate_sir"
  )

  expect_true(requireNamespace("SIRSimulation", quietly = TRUE))
  expect_setequal(getNamespaceExports("SIRSimulation"), expected_exports)
})

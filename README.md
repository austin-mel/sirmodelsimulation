# SIRSimulation

SIRSimulation is a small R package for running the existing SIR model
simulation functions.

Load the package with:

```r
library(SIRSimulation)
```

Public functions:

- `create_matrix()`
- `create_crnr_matrix()`
- `create_cntr_matrix()`
- `simulate_sir()`
- `simulate_inf_seq()`
- `simulate_many_runs()`
- `multiple_run_heatmap()`

## SIR simulation workflow

Create an initial population matrix with `0` for susceptible cells and `1` for
infected cells:

```r
initial <- create_matrix(row = 10, col = 10, start_inf = 0.1, seed = 94128)
```

Run one basic SIR simulation:

```r
result <- simulate_sir(prob_infect = 0.25, input_matrix = initial, seed = 94128)
result$steps
result$inf_prop
result$final_matrix
```

Inspect the per-step susceptible, infected, and recovered counts:

```r
result$history
```

Run repeated simulations and average the results:

```r
simulate_many_runs(prob_infect = 0.25, input_matrix = initial, runs = 10, seed = 94128)
```

Sweep across infection probabilities:

```r
simulate_inf_seq(input_matrix = initial, step = 0.1, seed = 94128)
```

Summarize how often each cell ended recovered across repeated SIR runs:

```r
heatmap_counts <- multiple_run_heatmap(
  prob_infect = 0.25,
  input_matrix = initial,
  runs = 10,
  plot = TRUE,
  seed = 94128
)
```

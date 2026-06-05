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

## State Encoding

Simulation matrices use numeric state values:

- `0`: susceptible
- `1`: infected
- `2`: recovered
- `3`: deceased

Initial matrices usually start with only susceptible and infected cells:

```r
initial <- create_matrix(row = 10, col = 10, start_inf = 0.1, seed = 94128)
```

## Models

Use `model` to choose the transition behavior for infected cells:

- `"SIR"`: infected cells become recovered.
- `"SIS"`: infected cells become susceptible again.
- `"SIRS"`: infected cells become recovered with probability `imm_prob`; otherwise they become susceptible again.

Optional mortality is available for all models with `allow_death = TRUE`.
When mortality is enabled, infected cells become deceased with probability
`fat_prob`; otherwise they follow the active model transition.

Run one SIR simulation:

```r
result <- simulate_sir(
  prob_infect = 0.25,
  input_matrix = initial,
  model = "SIR",
  seed = 94128
)
result$steps
result$inf_prop
result$final_matrix
```

Inspect the per-step susceptible, infected, and recovered counts:

```r
result$history
```

Run an SIRS simulation with a custom immunity probability:

```r
simulate_sir(
  prob_infect = 0,
  input_matrix = matrix(1, nrow = 1, ncol = 1),
  model = "SIRS",
  imm_prob = 0.7,
  seed = 94128
)
```

Enable mortality:

```r
simulate_sir(
  prob_infect = 0,
  input_matrix = matrix(1, nrow = 1, ncol = 1),
  model = "SIR",
  allow_death = TRUE,
  fat_prob = 1,
  seed = 94128
)
```

Run repeated simulations and average the results:

```r
simulate_many_runs(
  prob_infect = 0.25,
  input_matrix = initial,
  runs = 10,
  model = "SIR",
  seed = 94128
)
```

Sweep across infection probabilities:

```r
simulate_inf_seq(
  input_matrix = initial,
  step = 0.1,
  model = "SIR",
  seed = 94128
)
```

Summarize how often each cell ended recovered across repeated SIR runs:

```r
heatmap_counts <- multiple_run_heatmap(
  prob_infect = 0.25,
  input_matrix = initial,
  runs = 10,
  model = "SIR",
  plot = TRUE,
  seed = 94128
)
```

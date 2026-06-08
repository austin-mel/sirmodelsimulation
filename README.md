# SIRSsim

SIRSsim is a small R package for running SIR, SIS, and SIRS model simulation
functions.

Load the package with:

```r
library(SIRSsim)
```

Public functions:

- `create_random_matrix()`
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

Initial matrices usually start with only susceptible and infected cells.
Use a fractional `start_infected` value to treat it as an infection
probability:

```r
initial <- create_random_matrix(
  row = 10,
  col = 10,
  start_infected = 0.1,
  seed = 94128
)
```

Use a whole number to infect exactly that many randomly selected cells:

```r
initial_count <- create_random_matrix(
  row = 10,
  col = 10,
  start_infected = 8,
  seed = 94128
)
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

## Full Logs for Survival Analysis

Set `full_log = TRUE` to return one row per cell per simulation step. The log
includes step `0`, before any transitions occur.

```r
logged_result <- simulate_sir(
  prob_infect = 0.25,
  input_matrix = initial,
  model = "SIR",
  seed = 94128,
  full_log = TRUE
)

head(logged_result$full_log)
```

The `full_log` table includes cell identity, position, state, model settings,
and state indicators:

```r
names(logged_result$full_log)
```

Downstream survival-analysis code can derive per-cell endpoints such as first
infection, recovery, death, last observed step, event indicators, time to event,
and final state:

```r
cell_log <- logged_result$full_log

survival_ready <- do.call(
  rbind,
  lapply(split(cell_log, cell_log$cell_id), function(cell) {
    first_infected_step <- if (any(cell$was_infected)) min(cell$step[cell$was_infected]) else NA
    recovered_step <- if (any(cell$was_immune)) min(cell$step[cell$was_immune]) else NA
    death_step <- if (any(cell$was_deceased)) min(cell$step[cell$was_deceased]) else NA
    last_observed_step <- max(cell$step)

    data.frame(
      cell_id = cell$cell_id[1],
      first_infected_step = first_infected_step,
      recovered_step = recovered_step,
      death_step = death_step,
      last_observed_step = last_observed_step,
      event_infection = !is.na(first_infected_step),
      event_recovery = !is.na(recovered_step),
      event_death = !is.na(death_step),
      time_to_infection = first_infected_step,
      time_to_recovery = recovered_step,
      time_to_death = death_step,
      final_state = cell$state[cell$step == last_observed_step]
    )
  })
)
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

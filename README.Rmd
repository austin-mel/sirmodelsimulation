# SIRSsim

SIRSsim is a small R package for simulating how an infection can move through a
grid of cells. Each cell can be read as a person, location, or unit in a
population. At each step, infected cells may spread infection to nearby
susceptible cells, then move to their next state according to the selected
model.

These simulations are useful for teaching, exploration, and code experiments.
They are not calibrated public-health forecasts.

## Installation

Install the current package from GitHub:

```r
install.packages("remotes")
remotes::install_github("austin-mel/sirmodelsimulation")
```

The GitHub repository name is `sirmodelsimulation`, while the package loaded in
R is `SIRSsim`, as listed in `DESCRIPTION`.

```r
library(SIRSsim)
```

If you are working from a local clone, you can also install the package from
the project folder:

```r
remotes::install_local(".")
```

## What The Model Names Mean

SIRSsim supports three related disease-spread models:

- `SIR`: susceptible cells can become infected; infected cells then recover and
  stay recovered.
- `SIS`: susceptible cells can become infected; infected cells then become
  susceptible again, so they can be infected again later.
- `SIRS`: susceptible cells can become infected; infected cells recover with
  probability `imm_prob` or otherwise become susceptible again. In this package,
  recovered cells stay recovered once they reach that state.

Optional mortality can be layered onto any model with `allow_death = TRUE`.
When mortality is enabled, infected cells can become deceased with probability
`fat_prob` before following the model's nonfatal transition.

## Functions

Create starting matrices:

- `create_random_matrix()` infects an exact number of randomly selected cells.
- `create_center_matrix()` infects the center cell.
- `create_corner_matrix()` infects the four corner cells.
- `create_matrix()` is a backward-compatible wrapper around
  `create_random_matrix()`.

Run and summarize simulations:

- `simulate_sir()` runs one simulation until no infected cells remain.
- `simulate_many_runs()` repeats one scenario and returns average outcomes.
- `simulate_inf_seq()` sweeps across a sequence of infection probabilities.
- `multiple_run_heatmap()` counts how often each cell was ever infected across
  repeated runs.

## State Encoding

Simulation matrices use numeric state values:

| Value | State | Plain-English meaning |
| --- | --- | --- |
| `0` | susceptible | Healthy and able to become infected |
| `1` | infected | Currently infected and able to infect neighbors |
| `2` | recovered | No longer infected and treated as immune |
| `3` | deceased | Removed by optional mortality |

Initial matrices usually start with only susceptible and infected cells.
Infection can spread from any of the eight neighboring cells: horizontal,
vertical, or diagonal.

```r
initial <- create_random_matrix(
  row = 10,
  col = 10,
  start_infected = 8,
  seed = 94128
)

center_start <- create_center_matrix(row = 10, col = 10)
corner_start <- create_corner_matrix(row = 10, col = 10)
```

## Important Parameters

Most probabilities are written as numbers between `0` and `1`, where `0` means
"never" and `1` means "always."

| Concept | Argument | Meaning |
| --- | --- | --- |
| Infection probability | `prob_infect` | Chance that one infected neighbor infects a susceptible cell during a step. Multiple infected neighbors increase the overall infection chance. |
| Immunity probability | `imm_prob` | In `model = "SIRS"`, chance that an infected cell becomes recovered after its infectious step instead of becoming susceptible again. Ignored by `SIR` and `SIS`. |
| Fatality probability | `fat_prob` | Chance that an infected cell becomes deceased when `allow_death = TRUE`. Ignored when mortality is not enabled. |
| Starting matrix | `input_matrix` | The grid of starting states passed to simulation functions. Use the matrix helpers to create common starting patterns. |
| Seed | `seed` | Optional random seed for reproducible random matrices and simulation paths. |
| Full log | `full_log` | When `TRUE`, `simulate_sir()` includes one row per cell per step in the returned `full_log` data frame. |

## Reproducibility

Randomness appears when cells are chosen for a random starting matrix, when
infection spreads, when SIRS immunity is assigned, and when optional mortality
is enabled. Pass the same `seed` with the same inputs to reproduce the same
result.

```r
initial <- create_random_matrix(10, 10, start_infected = 8, seed = 94128)

result <- simulate_sir(
  prob_infect = 0.25,
  input_matrix = initial,
  model = "SIR",
  seed = 94128
)
```

For repeated-run helpers such as `simulate_many_runs()` and
`multiple_run_heatmap()`, one seed makes the whole sequence of runs
reproducible.

## Run One Simulation

```r
result <- simulate_sir(
  prob_infect = 0.25,
  input_matrix = initial,
  model = "SIR",
  seed = 94128,
  full_log = TRUE
)

result$steps
result$inf_prop
result$final_matrix
head(result$history)
```

## Reading Outputs

`simulate_sir()` returns a list with summary fields and, optionally, a detailed
log:

| Output | Meaning |
| --- | --- |
| `steps` | Number of simulation steps completed before no infected cells remained. Step `0` is the starting matrix and is included in `history` and `full_log`, not in this count. |
| `prob_infect` | Infection probability used for the run. |
| `inf_prop` | Proportion of cells that are recovered in the final matrix. In a basic `SIR` run without death this is often the final outbreak size; in `SIS`, `SIRS`, or mortality scenarios it is not the same as "ever infected." |
| `history` | Data frame with susceptible, infected, and recovered counts at each step, including step `0`. Deceased cells appear in `final_matrix` and `full_log`, but are not a separate `history` column. |
| `final_matrix` | Matrix of final cell states after the simulation stops. |
| `full_log` | Only returned when `full_log = TRUE`; a long data frame with one row per cell per step, including cell position, state, model settings, and state indicators. |

For "ever infected" questions, use `full_log` or
`multiple_run_heatmap()` instead of `inf_prop`.

## Full Logs For Cell-Level Analysis

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
names(logged_result$full_log)
```

The `full_log` table can be summarized into per-cell endpoints such as first
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

## More Simulation Workflows

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

Summarize how often each cell was infected at any point across repeated runs:

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

## Learn More

- The [package vignette](vignettes/sirsimulation.Rmd) walks through model
  choices, starting matrices, outputs, full logs, repeated runs, and heatmaps.
- The [test suite](tests/testthat/) contains executable examples of expected
  behavior, including matrix creation, infection spread, model transitions,
  reproducibility, and validation errors.

# SIRSsim Documentation

A documentation hub for an R package that simulates SIR, SIS, and SIRS disease-spread models on population matrices.

## What this project does

SIRSsim provides exported R functions for creating starting matrices, running single simulations, running repeated simulations, sweeping infection probabilities, and summarizing repeated infection patterns. The package supports SIR, SIS, and SIRS model options and optional mortality.

The package files document numeric state encoding, public functions, examples, manual pages, a vignette, and testthat coverage.

## Why it matters

The package gives users a small, inspectable simulation toolkit for experimenting with transition rules on a grid. Its documentation is organized around reproducible examples, state definitions, function references, and tests that check model behavior.

## Project contents

| Area | Files | What to read first |
|---|---|---|
| Package overview | [`../README.md`](../README.md), [`../README.Rmd`](../README.Rmd), [`../DESCRIPTION`](../DESCRIPTION) | Start here for the package purpose, public functions, state encoding, and examples. |
| Vignette | [`../vignettes/sirsimulation.Rmd`](../vignettes/sirsimulation.Rmd) | Walkthrough of matrices, models, mortality, history, full logs, repeated runs, and probability sweeps. |
| R source | [`../R/matrices.R`](../R/matrices.R), [`../R/simulation.R`](../R/simulation.R), [`../R/neighbors.R`](../R/neighbors.R), [`../R/plotting.R`](../R/plotting.R), [`../R/states.R`](../R/states.R) | Implementation files for matrix helpers, simulation logic, neighbor logic, plotting, and state constants. |
| Public API | [`../NAMESPACE`](../NAMESPACE), [`../man`](../man) | Exported functions and generated manual pages. |
| Tests | [`../tests/testthat`](../tests/testthat), [`../tests/testthat.R`](../tests/testthat.R) | Testthat checks for matrices, neighbors, simulations, full logs, mortality, and invalid inputs. |
| Assets | [`../assets/sirsim`](../assets/sirsim) | Images used by the package README or documentation. |

## How to read this documentation

1. Read the package `README.md` for the public function list and basic examples.
2. Read the vignette for a longer walkthrough of state encoding, model options, mortality, full logs, and repeated runs.
3. Use the files in `man/` for function-level documentation.
4. Review `tests/testthat/` when you want to see what behavior is checked by the package tests.

## Main source files used

<details>
<summary>Source files used to prepare this README</summary>

- `DESCRIPTION`
- `README.md`
- `README.Rmd`
- `NAMESPACE`
- `LICENSE`
- `LICENSE.md`
- `vignettes/sirsimulation.Rmd`
- `R/matrices.R`
- `R/simulation.R`
- `R/neighbors.R`
- `R/plotting.R`
- `R/states.R`
- `man/create_center_matrix.Rd`
- `man/create_corner_matrix.Rd`
- `man/create_matrix.Rd`
- `man/create_random_matrix.Rd`
- `man/multiple_run_heatmap.Rd`
- `man/simulate_inf_seq.Rd`
- `man/simulate_many_runs.Rd`
- `man/simulate_sir.Rd`
- `tests/testthat/test-matrices.R`
- `tests/testthat/test-neighbors.R`
- `tests/testthat/test-simulation.R`

</details>

## Known limits

- The package description frames the project as a simulator for disease-spread model options, not as a calibrated forecasting model.
- Simulation results depend on model choice, transition probabilities, starting matrix, random seed, and grid setup.
- Full logs are optional and must be requested with `full_log = TRUE`.
- The package is versioned as `0.0.0.9000` in `DESCRIPTION`.

---

<p align="center">
  <img src="./assets/signature.svg" alt="Austin Melendez signature" width="180" />
</p>

<p align="center">
  <strong>Austin Melendez</strong><br />
  <a href="https://austinmelendez.com">Portfolio</a> | <a href="https://github.com/austin-mel">GitHub</a> | <a href="https://github.com/austin-mel/sirmodelsimulation">Project repository</a>
</p>




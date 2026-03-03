# IPTW-BCH: Monte Carlo Evaluation of Combined Inverse Probability Weighting and BCH Method for Causal Inference with Latent Class Exposures

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Replication code for:

> **How Well Does Combining Inverse Probability Weighting with the BCH Method Work for Causal Inference with Latent Class Exposures? A Monte Carlo Evaluation**
> *Multivariate Behavioral Research* (under review)

---

## Overview

This repository contains all R code for a Monte Carlo simulation study evaluating the combined IPTW-BCH estimator for causal inference when the treatment variable is a latent class. The simulation crosses four sample sizes, three entropy conditions, and three confounding strengths (36 conditions × 500 replications = 18,000 total), comparing four estimators:

| Estimator | Classification Correction | Confounding Correction |
|-----------|--------------------------|------------------------|
| IPTW-BCH (oracle) | BCH weights | IPW from true class |
| IPTW-BCH (estimated) | BCH weights | IPW from modal class |
| BCH only | BCH weights | None |
| Classify-analyze | None | None |

---

## Requirements

### R (version ≥ 4.0)

```r
install.packages(c("MplusAutomation", "nnet", "MASS", "parallel",
                   "dplyr", "tidyr", "ggplot2", "hrbrthemes",
                   "patchwork", "knitr", "xtable"))
```

### Mplus (version ≥ 8)

A licensed installation of [Mplus](https://www.statmodel.com/) is required. The code uses `MplusAutomation` to generate input files and parse output programmatically.

---

## Simulation Design

| Factor | Levels | Values |
|--------|--------|--------|
| Sample size (*N*) | 4 | 250, 500, 1,000, 2,000 |
| Classification entropy | 3 | High (≈ 0.8), Medium (≈ 0.7), Low (≈ 0.6) |
| Confounding strength (γ) | 3 | None (0), Moderate (0.5), Strong (1.0) |
| **Total conditions** | **36** | 36 × 500 = 18,000 replications |

**Data generating process:**
- *K* = 3 latent classes (proportions ≈ 40%, 35%, 25%)
- *J* = 10 binary indicators
- 3 continuous confounders X = (X₁, X₂, X₃)
- Outcome: Y = α_c + 0.3X₁ + 0.2X₂ + 0.1X₃ + ε, ε ~ N(0,1)
- True ATEs: ATE(C₂ vs C₁) = 0.5, ATE(C₃ vs C₁) = 1.0

---

## File Structure

```
code/
├── 00_master.R            # Run entire pipeline from a single call
├── 01_dgp_functions.R     # Data generating process (DGP) functions
├── 02_iptw_functions.R    # Propensity score estimation and IPW weight computation
├── 03_mplus_syntax.R      # Mplus input file generators for all 4 estimators
├── 04_simulation_runner.R # Monte Carlo orchestration (single condition)
├── 05_results_analysis.R  # Performance metrics: bias, RMSE, coverage, SE ratio
├── 06_figure_generation.R # Publication-quality figures (500 DPI)
├── generate_tables.R      # LaTeX table generation
└── run_expanded_sim.R     # Full 36-condition simulation runner (parallelized)
```

### Script descriptions

**`00_master.R`**
Top-level orchestration. Sets paths, loads all helper scripts, and calls `run_expanded_sim.R` to execute the full simulation. Adjust `N_CORES` and `OUTPUT_DIR` at the top of this file.

**`01_dgp_functions.R`**
Functions for generating one simulation replication: confounders from a multivariate normal, latent class membership via multinomial logistic model (controlled by γ), binary indicators with class-specific item response probabilities (3 entropy conditions), and the continuous outcome.

**`02_iptw_functions.R`**
Propensity score estimation via multinomial logistic regression (`nnet::multinom`), stabilized IPW weight computation, 99th-percentile trimming, and weight diagnostic summary (mean, SD, max, CV, n_trimmed).

**`03_mplus_syntax.R`**
Functions that write Mplus `.inp` files for each of the four estimators:
- `write_step1_syntax()`: Unconditional LCA with BCH weight extraction (`SAVE = BCHWEIGHTS`)
- `write_step3_bch_iptw_syntax()`: BCH + IPTW (`TRAINING(BCH)` + `WEIGHT`)
- `write_step3_bch_only_syntax()`: BCH without IPTW
- `write_step3_classify_syntax()`: Classify-analyze (`KNOWNCLASS`)

**`04_simulation_runner.R`**
Runs one full replication: generates data → estimates LCA → applies label-switching correction → fits all four estimators → extracts ATE estimates, SEs, CIs, and weight diagnostics. Returns a single-row data frame.

**`05_results_analysis.R`**
Aggregates replications across conditions. Computes bias, RMSE, 95% CI coverage, SE ratio, and Monte Carlo standard errors (MCSE) for each metric.

**`06_figure_generation.R`**
Generates all manuscript figures using `ggplot2`:
- Figure 2: Bias by entropy × confounding, faceted
- Figure 3: Coverage rates across conditions
- Figure 4: RMSE by sample size and method
- Figure 5: Bias heatmap (entropy × confounding grid)
- Figure S1: SE ratio diagnostic

**`generate_tables.R`**
Generates LaTeX tables from the simulation output. Reads `simulation_output_v2/all_results.rds` and writes `.tex` files to `../tables/`.

**`run_expanded_sim.R`**
Parallelized runner for all 36 conditions. Uses `parallel::mclapply` (set `mc.cores` to available cores). Saves results as `simulation_output_v2/all_results.rds`.

---

## How to Run

### Full replication

```r
# In R, from the code/ directory:
source("00_master.R")
```

This will:
1. Run all 36 conditions × 500 replications in parallel
2. Save raw results to `simulation_output_v2/all_results.rds`
3. Compute performance metrics
4. Generate figures to `../figures/`
5. Generate LaTeX tables to `../tables/`

**Estimated runtime:** ~28 hours wall time with 9 cores (≈250 CPU-hours total)

### Quick test (3 replications per condition)

```r
source("01_dgp_functions.R")
source("02_iptw_functions.R")
source("03_mplus_syntax.R")
source("04_simulation_runner.R")

# Test one condition
result <- run_one_rep(rep_id = 1, N = 500, entropy = "high",
                      gamma = 0.5, output_dir = "simulation_output_test/")
```

### Reproduce figures and tables only

If `simulation_output_v2/all_results.rds` already exists:

```r
source("05_results_analysis.R")
source("06_figure_generation.R")
source("generate_tables.R")
```

---

## Key Mplus Syntax

The combined IPTW-BCH estimator uses Mplus Step 3 syntax:

```
VARIABLE:
  NAMES ARE U1-U10 Y IPW W1 W2 W3 MLC;
  USEVARIABLES ARE Y W1 W2 W3 IPW;
  CLASSES = C(3);
  TRAINING = W1 W2 W3(BCH);   ! BCH weights correct for classification error
  WEIGHT = IPW;               ! IPW weights adjust for confounding

ANALYSIS:
  TYPE = MIXTURE;
  ESTIMATOR = MLR;
  STARTS = 0;                 ! BCH weights fix class assignment; no random starts needed
```

Full syntax for all four estimators appears in Appendix E of the Supplementary Material.

---

## Output

After a complete run, the following files are produced:

```
simulation_output_v2/
└── all_results.rds      # 144,000-row data frame (36 conditions × 500 reps × 8 rows)

../figures/
├── Figure2_bias.png
├── Figure3_coverage.png
├── Figure4_rmse.png
├── Figure5_bias_heatmap.png
└── FigureS1_se_ratio.png

../tables/
├── table2_main_results_c2c1.tex
├── table2_main_results_c3c1.tex
├── table3_se_calibration.tex
└── table4_gamma0.tex
```

---

## Citation

If you use this code, please cite:

```bibtex
@article{iptw_bch_2025,
  author  = {Lu, Z{\'{e}}kai (Zachary)},
  title   = {How Well Does Combining Inverse Probability Weighting with the {BCH} Method
             Work for Causal Inference with Latent Class Exposures? {A} Monte Carlo Evaluation},
  journal = {Multivariate Behavioral Research},
  year    = {2025},
  note    = {Under review}
}
```

**Author:** Zékai (Zachary) Lu, Department of Sociology, McGill University

---

## License

Copyright (c) 2025 Zékai (Zachary) Lu

This project is licensed under the MIT License — see the [`LICENSE`](LICENSE) file for details. You are free to use, modify, and distribute this code for any purpose, including academic and commercial use, provided the original copyright notice is retained.

# ============================================================
# File: 05_results_analysis.R
# Purpose: Reusable functions for computing simulation
#          performance metrics with MCSEs
# Depends: 01_dgp_functions.R (for TRUE_ATE constants)
# ============================================================

library(dplyr)
library(tidyr)

# Source DGP constants
source("01_dgp_functions.R")

# ------------------------------------------------------------------
# True parameter values (from DGP)
# ------------------------------------------------------------------
TRUE_VALUES <- data.frame(
  class_pair = c("C2_vs_C1", "C3_vs_C1"),
  true_ate = c(TRUE_ATE_C2_vs_C1, TRUE_ATE_C3_vs_C1),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# load_results(): Load and validate simulation output
# ------------------------------------------------------------------
load_results <- function(results_path = "simulation_output_v2/all_results.rds") {

  if (!file.exists(results_path)) {
    stop("Results file not found: ", results_path,
         "\nRun the simulation first.")
  }

  results <- readRDS(results_path)

  required_cols <- c("rep_id", "N", "entropy", "gamma", "model",
                     "class_pair", "est_ate", "se", "ci_lower",
                     "ci_upper", "converged")

  missing <- setdiff(required_cols, names(results))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }

  results_conv <- results[results$converged == TRUE, ]

  cat(sprintf("Loaded %d total rows, %d converged (%.1f%%)\n",
              nrow(results), nrow(results_conv),
              100 * nrow(results_conv) / nrow(results)))

  results_conv
}

# ------------------------------------------------------------------
# compute_performance(): Performance metrics with MCSEs
# ------------------------------------------------------------------
compute_performance <- function(results) {

  results <- merge(results, TRUE_VALUES, by = "class_pair")

  perf <- results %>%
    group_by(N, entropy, gamma, model, class_pair, true_ate) %>%
    summarise(
      n_reps = n(),
      mean_est = mean(est_ate, na.rm = TRUE),
      bias = mean(est_ate, na.rm = TRUE) - first(true_ate),
      rel_bias_pct = 100 * (mean(est_ate, na.rm = TRUE) - first(true_ate)) /
                     first(true_ate),
      emp_se = sd(est_ate, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      se_ratio = mean(se, na.rm = TRUE) / sd(est_ate, na.rm = TRUE),
      mse = mean((est_ate - first(true_ate))^2, na.rm = TRUE),
      rmse = sqrt(mean((est_ate - first(true_ate))^2, na.rm = TRUE)),
      coverage = mean(ci_lower <= first(true_ate) &
                      ci_upper >= first(true_ate), na.rm = TRUE),
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      power = mean(ci_lower > 0 | ci_upper < 0, na.rm = TRUE),

      # MCSEs
      mcse_bias = sd(est_ate, na.rm = TRUE) / sqrt(n()),
      mcse_coverage = sqrt(
        mean(ci_lower <= first(true_ate) & ci_upper >= first(true_ate),
             na.rm = TRUE) *
        (1 - mean(ci_lower <= first(true_ate) & ci_upper >= first(true_ate),
                  na.rm = TRUE)) / n()
      ),

      # Weight diagnostics
      mean_w_mean = mean(w_mean, na.rm = TRUE),
      mean_w_max = mean(w_max, na.rm = TRUE),
      mean_w_cv = mean(w_cv, na.rm = TRUE),
      mean_w_trimmed = mean(w_n_trimmed, na.rm = TRUE),

      .groups = "drop"
    )

  perf
}

# ------------------------------------------------------------------
# compute_convergence_rates()
# ------------------------------------------------------------------
compute_convergence_rates <- function(results_all) {
  results_all %>%
    group_by(N, entropy, gamma, model) %>%
    summarise(
      total = n(),
      n_converged = sum(converged, na.rm = TRUE),  # Renamed to avoid collision with input column
      conv_rate = n_converged / total,              # Use the renamed column
      .groups = "drop"
    )
}

# ------------------------------------------------------------------
# run_analysis(): Master analysis pipeline
# ------------------------------------------------------------------
run_analysis <- function(results_path = "simulation_output_v2/all_results.rds") {

  results_all <- readRDS(results_path)
  results <- load_results(results_path)

  perf <- compute_performance(results)
  conv <- compute_convergence_rates(results_all)

  # Save
  out_dir <- dirname(results_path)
  saveRDS(perf, file.path(out_dir, "performance_metrics.rds"))
  write.csv(perf, file.path(out_dir, "performance_metrics.csv"),
            row.names = FALSE)

  cat("\nPerformance metrics saved.\n")
  invisible(perf)
}

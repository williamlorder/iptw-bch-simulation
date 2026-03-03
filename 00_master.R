# ============================================================
# File: 00_master.R
# Purpose: Master orchestration script — runs entire simulation
#          pipeline from data generation through results analysis
# Usage:  source("00_master.R")   OR   Rscript 00_master.R
# ============================================================

cat("============================================================\n")
cat("  IPTW-BCH Monte Carlo Simulation\n")
cat("  Causal Inference with Latent Class Exposures\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
N_CORES <- parallel::detectCores() - 1  # Leave 1 core free
if (N_CORES < 1) N_CORES <- 1

cat(sprintf("Detected %d cores, using %d for simulation.\n\n",
            parallel::detectCores(), N_CORES))

# ------------------------------------------------------------------
# Step 1: Source all component scripts
# ------------------------------------------------------------------
cat("--- Loading component scripts ---\n")
source("01_dgp_functions.R")
source("02_iptw_functions.R")
source("03_mplus_syntax.R")
source("04_simulation_runner.R")
source("05_results_analysis.R")
source("06_figure_generation.R")
cat("All scripts loaded.\n\n")

# ------------------------------------------------------------------
# Step 2: Verify DGP
# ------------------------------------------------------------------
cat("--- Verifying DGP ---\n")
compute_true_potential_outcomes()

test_dat <- generate_data(N = 1000, entropy = "high", gamma = 0.5, seed = 42)
cat(sprintf("\nTest data: N=%d, K=%d classes\n",
            nrow(test_dat$data), length(unique(test_dat$data$true_class))))
cat(sprintf("Class proportions: %.3f, %.3f, %.3f\n",
            test_dat$params$true_props[1],
            test_dat$params$true_props[2],
            test_dat$params$true_props[3]))

# Quick IPTW test
test_ps <- estimate_propensity_scores(test_dat$data, "true_class")
cat(sprintf("PS model converged: %s\n", test_ps$converged))
test_ipw <- compute_ipw_weights(test_ps$ps_matrix, test_dat$data$true_class)
cat(sprintf("IPW weight summary: mean=%.3f, sd=%.3f, max=%.3f\n",
            test_ipw$summary$mean, test_ipw$summary$sd, test_ipw$summary$max))
cat("\n")

# ------------------------------------------------------------------
# Step 3: Run Monte Carlo simulation
# ------------------------------------------------------------------
cat("--- Starting Monte Carlo simulation ---\n")
cat(sprintf("Conditions: %d\n", nrow(SIM_CONDITIONS)))
cat(sprintf("Replications per condition: %d\n", N_REPS))
cat(sprintf("Total replications: %d\n", nrow(SIM_CONDITIONS) * N_REPS))
cat(sprintf("Parallel cores: %d\n", N_CORES))

start_time <- Sys.time()

results <- run_simulation(
  conditions = SIM_CONDITIONS,
  n_reps = N_REPS,
  base_dir = "simulation_output",
  n_cores = N_CORES
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")
cat(sprintf("\nSimulation completed in %.1f hours.\n\n", as.numeric(elapsed)))

# ------------------------------------------------------------------
# Step 4: Analyze results and generate tables
# ------------------------------------------------------------------
cat("--- Analyzing results ---\n")
perf <- run_analysis("simulation_output/all_results.rds")

# ------------------------------------------------------------------
# Step 5: Generate figures
# ------------------------------------------------------------------
cat("\n--- Generating figures ---\n")
generate_all_figures("simulation_output/performance_metrics.rds")

# ------------------------------------------------------------------
# Done
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  Pipeline complete.\n")
cat(sprintf("  Results: simulation_output/\n"))
cat(sprintf("  Tables:  ../tables/\n"))
cat(sprintf("  Figures: ../figures/\n"))
cat("============================================================\n")

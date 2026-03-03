# ============================================================
# Run expanded simulation: 24 conditions × 500 reps
# 9-core parallel execution with real-time progress
# ============================================================
setwd("/Users/zekai/WorkPlace/LAB/MPLUS方法论论文/initial/scripts")

source("01_dgp_functions.R")
source("02_iptw_functions.R")
source("03_mplus_syntax.R")
source("04_simulation_runner.R")

library(parallel)

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
# Set PILOT_MODE = TRUE to test 3 reps per condition before full run
PILOT_MODE <- FALSE
N_REPS <- if (PILOT_MODE) 3 else 500
N_CORES <- 9
BATCH_SIZE <- if (PILOT_MODE) 3 else 50  # reps per batch for progress
OUTPUT_DIR <- if (PILOT_MODE) "simulation_pilot" else "simulation_output_v2"

# 36-condition factorial design: 4N × 3entropy × 3gamma
SIM_CONDITIONS <- expand.grid(
  N = c(250, 500, 1000, 2000),
  entropy = c("high", "medium", "low"),
  gamma = c(0, 0.5, 1.0),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# Header
# ------------------------------------------------------------------
cat("============================================================\n")
if (PILOT_MODE) {
  cat("PILOT TEST: 36 conditions × 3 reps (testing before full run)\n")
} else {
  cat("FULL SIMULATION: 36 conditions × 500 reps = 18,000 total\n")
}
cat("============================================================\n")
cat("Conditions:", nrow(SIM_CONDITIONS), "\n")
cat("Reps per condition:", N_REPS, "\n")
cat("Total replications:", nrow(SIM_CONDITIONS) * N_REPS, "\n")
cat("Parallel cores:", N_CORES, "\n")
cat("Estimators: BCH_IPTW_oracle, BCH_IPTW_est, BCH_only, Classify\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
sim_start <- Sys.time()

# ------------------------------------------------------------------
# Main loop: conditions sequential, reps parallel
# ------------------------------------------------------------------
for (cond_idx in 1:nrow(SIM_CONDITIONS)) {
  cond <- SIM_CONDITIONS[cond_idx, ]
  cond_label <- paste0("N", cond$N, "_", cond$entropy, "_g", cond$gamma)
  cond_dir <- file.path(OUTPUT_DIR, cond_label)
  dir.create(cond_dir, showWarnings = FALSE, recursive = TRUE)

  results_file <- file.path(cond_dir, "results.rds")

  # Resume support: check for existing results
  start_rep <- 1
  all_results <- list()

  if (file.exists(results_file)) {
    existing <- tryCatch(readRDS(results_file), error = function(e) NULL)
    if (!is.null(existing) && nrow(existing) > 0) {
      n_done <- length(unique(existing$rep_id))
      if (n_done >= N_REPS) {
        cat(sprintf("[%d/%d] %s: Already complete (%d reps). Skipping.\n",
                    cond_idx, nrow(SIM_CONDITIONS), cond_label, n_done))
        next
      }
      cat(sprintf("[%d/%d] %s: Resuming from rep %d\n",
                  cond_idx, nrow(SIM_CONDITIONS), cond_label, n_done + 1))
      start_rep <- n_done + 1
      # Rebuild results list from existing data
      for (rid in unique(existing$rep_id)) {
        all_results[[as.character(rid)]] <- existing[existing$rep_id == rid, ]
      }
    }
  }

  cat(sprintf("[%d/%d] %s: Running reps %d to %d...\n",
              cond_idx, nrow(SIM_CONDITIONS), cond_label, start_rep, N_REPS))
  t0 <- Sys.time()

  # Process in batches for progress reporting
  remaining_reps <- start_rep:N_REPS
  batch_starts <- seq(1, length(remaining_reps), by = BATCH_SIZE)

  for (b in batch_starts) {
    batch_indices <- b:min(b + BATCH_SIZE - 1, length(remaining_reps))
    batch_reps <- remaining_reps[batch_indices]

    # Run batch in parallel using mclapply (fork-based, macOS)
    batch_results <- mclapply(batch_reps, function(r) {
      tryCatch(
        run_single_replication(
          rep_id = r,
          N = cond$N,
          entropy = cond$entropy,
          gamma = cond$gamma,
          work_dir = cond_dir
        ),
        error = function(e) {
          make_failed_result(r, cond$N, cond$entropy, cond$gamma)
        }
      )
    }, mc.cores = min(N_CORES, length(batch_reps)))

    # Collect batch results
    for (i in seq_along(batch_reps)) {
      if (!is.null(batch_results[[i]])) {
        all_results[[as.character(batch_reps[i])]] <- batch_results[[i]]
      }
    }

    # Progress report
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    n_done_now <- length(all_results)
    n_converged <- sum(sapply(all_results, function(x) {
      any(x$converged, na.rm = TRUE)
    }))
    cat(sprintf("  Rep %d/%d (%.1f min elapsed, %d/%d converged)\n",
                max(batch_reps), N_REPS, elapsed, n_converged, n_done_now))

    # Checkpoint save
    checkpoint <- do.call(rbind, all_results)
    rownames(checkpoint) <- NULL
    saveRDS(checkpoint, results_file)
    cat(sprintf("  Checkpoint saved: %d reps\n", n_done_now))
  }

  # Condition summary
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  final_cond <- do.call(rbind, all_results)
  n_converged <- length(unique(
    final_cond$rep_id[final_cond$converged == TRUE]
  ))
  cat(sprintf("[%d/%d] %s: Complete. %d reps in %.1f min (%.0f%% converged)\n\n",
              cond_idx, nrow(SIM_CONDITIONS), cond_label,
              length(all_results), elapsed,
              100 * n_converged / length(all_results)))
}

# ------------------------------------------------------------------
# Combine all results
# ------------------------------------------------------------------
cat("Combining all results...\n")
all_files <- list.files(OUTPUT_DIR, pattern = "results\\.rds$",
                        recursive = TRUE, full.names = TRUE)
# Exclude the combined file itself
all_files <- all_files[!grepl("all_results", all_files)]

combined <- do.call(rbind, lapply(all_files, readRDS))
rownames(combined) <- NULL
saveRDS(combined, file.path(OUTPUT_DIR, "all_results.rds"))
write.csv(combined, file.path(OUTPUT_DIR, "all_results.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
total_elapsed <- round(as.numeric(difftime(Sys.time(), sim_start,
                                            units = "hours")), 1)
cat("\n============================================================\n")
cat("SIMULATION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Total rows: %d\n", nrow(combined)))
cat(sprintf("Unique reps: %d\n", length(unique(paste(
  combined$N, combined$entropy, combined$gamma, combined$rep_id
)))))
cat(sprintf("Overall convergence: %.1f%%\n",
            100 * mean(combined$converged, na.rm = TRUE)))
cat(sprintf("Total wall time: %.1f hours\n", total_elapsed))
cat(sprintf("Results saved: %s/all_results.rds\n", OUTPUT_DIR))
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

if (PILOT_MODE) {
  cat("\n*** PILOT TEST COMPLETE ***\n")
  cat("Review results above. If all conditions converged,\n")
  cat("set PILOT_MODE <- FALSE and re-run for full simulation.\n")
}

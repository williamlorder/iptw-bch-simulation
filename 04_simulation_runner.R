# ============================================================
# File: 04_simulation_runner.R
# Purpose: Monte Carlo simulation orchestration with
#          4 estimators, label-switching correction,
#          weight diagnostics, and result extraction
# Depends: 01_dgp_functions.R, 02_iptw_functions.R, 03_mplus_syntax.R
# ============================================================

library(MplusAutomation)
library(nnet)
library(parallel)

# Source helper scripts
source("01_dgp_functions.R")
source("02_iptw_functions.R")
source("03_mplus_syntax.R")

# ------------------------------------------------------------------
# Simulation design: 36-condition factorial grid
# ------------------------------------------------------------------
SIM_CONDITIONS <- expand.grid(
  N = c(250, 500, 1000, 2000),
  entropy = c("high", "medium", "low"),
  gamma = c(0, 0.5, 1.0),
  stringsAsFactors = FALSE
)

N_REPS <- 500  # Monte Carlo replications per condition
K <- 3         # Number of latent classes

# Model names for the 4 estimators
MODEL_NAMES <- c("BCH_IPTW_oracle", "BCH_IPTW_est", "BCH_only", "Classify")

# ------------------------------------------------------------------
# label_switch_greedy(): Resolve label switching
# ------------------------------------------------------------------
#' Correct label switching in LCA results
#'
#' For K=3, enumerates all 6 permutations. For general K, uses
#' greedy matching by minimizing squared differences to true means.
#'
#' @param estimated_means  Vector of K estimated class-specific means
#' @param true_means       Vector of K true class-specific means (ALPHA)
#' @return Integer vector: permutation mapping estimated -> true classes
label_switch_greedy <- function(estimated_means, true_means) {

  K <- length(true_means)
  if (length(estimated_means) != K) return(1:K)

  # For K=3, enumerate all 6 permutations directly
  if (K == 3) {
    perms <- list(
      c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
      c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
    )

    best_perm <- c(1, 2, 3)
    best_cost <- Inf

    for (perm in perms) {
      cost <- sum((estimated_means[perm] - true_means)^2)
      if (cost < best_cost) {
        best_cost <- cost
        best_perm <- perm
      }
    }
    return(best_perm)
  }

  # For general K: greedy matching
  remaining <- 1:K
  perm <- integer(K)

  for (k in 1:K) {
    diffs <- abs(estimated_means[remaining] - true_means[k])
    best_idx <- which.min(diffs)
    perm[k] <- remaining[best_idx]
    remaining <- remaining[-best_idx]
  }

  perm
}

# ------------------------------------------------------------------
# run_single_replication(): One MC replication with 4 estimators
# ------------------------------------------------------------------
#' Run a single Monte Carlo replication
#'
#' Estimators:
#'   1. BCH_IPTW_oracle: BCH + IPW from true-class PS (oracle)
#'   2. BCH_IPTW_est:    BCH + IPW from modal-class PS (realistic)
#'   3. BCH_only:        BCH weights, no IPW
#'   4. Classify:        Modal class, no corrections
#'
#' @param rep_id     Replication number
#' @param N          Sample size
#' @param entropy    "high" or "low"
#' @param gamma      Confounding strength
#' @param work_dir   Working directory for Mplus files
#' @return data.frame with results from all four models
run_single_replication <- function(rep_id, N, entropy, gamma, work_dir) {

  # Create unique subdirectory for this replication
  rep_dir <- file.path(work_dir, paste0("rep_", rep_id))
  dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)

  result <- tryCatch({

    # ---- 1. Generate data ----
    sim <- generate_data(N = N, entropy = entropy, gamma = gamma,
                         seed = rep_id * 10000 + N)
    dat <- sim$data

    # ---- 2. Compute oracle IPTW (true class) ----
    ps_oracle <- estimate_propensity_scores(dat, class_var = "true_class")

    if (is.null(ps_oracle$ps_matrix) || !ps_oracle$converged) {
      return(make_failed_result(rep_id, N, entropy, gamma))
    }

    ipw_oracle <- compute_ipw_weights(
      ps_matrix = ps_oracle$ps_matrix,
      assigned = dat$true_class,
      stabilize = TRUE, trim = 0.99
    )

    oracle_diag <- compute_weight_diagnostics(ipw_oracle)

    # Covariate balance diagnostic (stored but not blocking)
    balance_smd <- tryCatch({
      bal <- check_balance(dat, dat$true_class, ipw_oracle$weights)
      max(abs(bal$smd_weighted))
    }, error = function(e) NA_real_)

    # ---- 3. Prepare Mplus data (with oracle IPW) ----
    data_file <- file.path(rep_dir, "simdata.dat")
    prepare_mplus_data(dat, ipw_weights = ipw_oracle$weights,
                       file_path = data_file)

    # ---- 4. Step 1: LCA with BCH weight saving ----
    bch_output_file <- file.path(rep_dir, "bch_data.dat")

    step1_syntax <- generate_step1_syntax(
      data_file = "simdata.dat", K = K,
      output_file = "bch_data.dat", has_ipw = TRUE
    )

    step1_inp <- file.path(rep_dir, "step1.inp")
    step1_results <- write_and_run_mplus(step1_syntax, step1_inp)

    if (is.null(step1_results)) {
      return(make_failed_result(rep_id, N, entropy, gamma))
    }

    obs_entropy <- extract_entropy(step1_results)

    if (!file.exists(bch_output_file)) {
      return(make_failed_result(rep_id, N, entropy, gamma))
    }

    # ---- 5. Extract modal class, compute estimated IPTW ----
    bch_data <- tryCatch(
      read.table(bch_output_file, header = FALSE),
      error = function(e) NULL
    )

    # SAVEDATA output columns (SAVE = BCHWEIGHTS CPROB):
    # 1-10: U1-U10, 11: Y, 12: IPW, 13-15: BCHW1-BCHW3, 16-18: CPROB1-CPROB3, 19: MLC
    MLC_COL <- 19  # Most Likely Class column
    IPW_COL <- 12  # IPW weight column

    # Default: no estimated weights
    est_diag <- data.frame(
      w_mean = NA, w_sd = NA, w_max = NA, w_cv = NA, w_n_trimmed = NA
    )
    has_est_iptw <- FALSE

    if (!is.null(bch_data) && nrow(bch_data) == N && ncol(bch_data) >= MLC_COL) {
      mlc <- bch_data[, MLC_COL]

      # Handle Mplus missing value format
      if (is.character(mlc)) mlc <- as.numeric(mlc)
      mlc <- round(mlc)

      # Need at least 2 classes represented for PS estimation
      if (length(unique(na.omit(mlc))) >= 2) {
        dat$modal_class <- mlc

        ps_est <- estimate_propensity_scores(dat, class_var = "modal_class")

        if (!is.null(ps_est$ps_matrix) && ps_est$converged) {
          ipw_est <- compute_ipw_weights(
            ps_matrix = ps_est$ps_matrix,
            assigned = dat$modal_class,
            stabilize = TRUE, trim = 0.99
          )

          est_diag <- compute_weight_diagnostics(ipw_est)

          # Create modified BCH data with estimated IPW in IPW column
          bch_data_est <- bch_data
          bch_data_est[, IPW_COL] <- ipw_est$weights
          est_file <- file.path(rep_dir, "bch_data_est.dat")
          write.table(bch_data_est, est_file,
                      sep = "\t", row.names = FALSE, col.names = FALSE)
          has_est_iptw <- TRUE
        }
      }
    }

    # ---- 6. Run Step 3 models ----
    # 3a: BCH + oracle IPTW
    step3a_syntax <- generate_step3_bch_iptw_syntax(
      data_file = "bch_data.dat", K = K
    )
    step3a_inp <- file.path(rep_dir, "step3_bch_iptw.inp")
    step3a_results <- write_and_run_mplus(step3a_syntax, step3a_inp)

    # 3b: BCH only
    step3b_syntax <- generate_step3_bch_only_syntax(
      data_file = "bch_data.dat", K = K
    )
    step3b_inp <- file.path(rep_dir, "step3_bch_only.inp")
    step3b_results <- write_and_run_mplus(step3b_syntax, step3b_inp)

    # 3c: Classify-analyze
    step3c_syntax <- generate_classify_analyze_syntax(
      data_file = "bch_data.dat", K = K
    )
    step3c_inp <- file.path(rep_dir, "step3_classify.inp")
    step3c_results <- write_and_run_mplus(step3c_syntax, step3c_inp)

    # 3d: BCH + estimated IPTW
    step3d_results <- NULL
    if (has_est_iptw) {
      step3d_syntax <- generate_step3_bch_iptw_syntax(
        data_file = "bch_data_est.dat", K = K
      )
      step3d_inp <- file.path(rep_dir, "step3_bch_iptw_est.inp")
      step3d_results <- write_and_run_mplus(step3d_syntax, step3d_inp)
    }

    # ---- 7. Extract and align results ----
    na_diag <- data.frame(
      w_mean = NA, w_sd = NA, w_max = NA, w_cv = NA, w_n_trimmed = NA
    )

    models <- list(
      list(name = "BCH_IPTW_oracle", res = step3a_results, diag = oracle_diag),
      list(name = "BCH_IPTW_est",    res = step3d_results, diag = est_diag),
      list(name = "BCH_only",        res = step3b_results, diag = na_diag),
      list(name = "Classify",        res = step3c_results, diag = na_diag)
    )

    results_list <- list()

    for (model_info in models) {
      means_df <- extract_class_means(model_info$res, K)

      if (is.null(means_df) || nrow(means_df) < K) {
        results_list[[length(results_list) + 1]] <- data.frame(
          rep_id = rep_id, N = N, entropy = entropy, gamma = gamma,
          model = model_info$name,
          class_pair = c("C2_vs_C1", "C3_vs_C1"),
          est_ate = NA, se = NA,
          ci_lower = NA, ci_upper = NA,
          obs_entropy = obs_entropy,
          converged = FALSE,
          w_mean = model_info$diag$w_mean,
          w_sd = model_info$diag$w_sd,
          w_max = model_info$diag$w_max,
          w_cv = model_info$diag$w_cv,
          w_n_trimmed = model_info$diag$w_n_trimmed,
          max_smd_weighted = balance_smd
        )
        next
      }

      # Label switching correction
      est_means <- means_df$mean
      perm <- label_switch_greedy(est_means, ALPHA)
      means_df <- means_df[perm, ]
      means_df$class <- 1:K

      # Compute estimated ATEs
      ate_c2_c1 <- means_df$mean[2] - means_df$mean[1]
      ate_c3_c1 <- means_df$mean[3] - means_df$mean[1]

      # SE of ATE via delta method
      se_c2_c1 <- sqrt(means_df$se[2]^2 + means_df$se[1]^2)
      se_c3_c1 <- sqrt(means_df$se[3]^2 + means_df$se[1]^2)

      results_list[[length(results_list) + 1]] <- data.frame(
        rep_id = rep_id, N = N, entropy = entropy, gamma = gamma,
        model = model_info$name,
        class_pair = c("C2_vs_C1", "C3_vs_C1"),
        est_ate = c(ate_c2_c1, ate_c3_c1),
        se = c(se_c2_c1, se_c3_c1),
        ci_lower = c(ate_c2_c1 - 1.96 * se_c2_c1,
                     ate_c3_c1 - 1.96 * se_c3_c1),
        ci_upper = c(ate_c2_c1 + 1.96 * se_c2_c1,
                     ate_c3_c1 + 1.96 * se_c3_c1),
        obs_entropy = obs_entropy,
        converged = TRUE,
        w_mean = model_info$diag$w_mean,
        w_sd = model_info$diag$w_sd,
        w_max = model_info$diag$w_max,
        w_cv = model_info$diag$w_cv,
        w_n_trimmed = model_info$diag$w_n_trimmed,
        max_smd_weighted = balance_smd
      )
    }

    do.call(rbind, results_list)

  }, error = function(e) {
    make_failed_result(rep_id, N, entropy, gamma)
  })

  # Clean up Mplus temp files
  unlink(rep_dir, recursive = TRUE)

  result
}

# ------------------------------------------------------------------
# make_failed_result(): Helper for failed replications
# ------------------------------------------------------------------
make_failed_result <- function(rep_id, N, entropy, gamma) {
  data.frame(
    rep_id = rep_id, N = N, entropy = entropy, gamma = gamma,
    model = rep(MODEL_NAMES, each = 2),
    class_pair = rep(c("C2_vs_C1", "C3_vs_C1"), length(MODEL_NAMES)),
    est_ate = NA, se = NA,
    ci_lower = NA, ci_upper = NA,
    obs_entropy = NA,
    converged = FALSE,
    w_mean = NA, w_sd = NA, w_max = NA, w_cv = NA, w_n_trimmed = NA,
    max_smd_weighted = NA
  )
}

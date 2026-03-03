# ============================================================
# File: 02_iptw_functions.R
# Purpose: Propensity score estimation and IPW weight calculation
#          for causal inference with latent class exposures
# Input:  Simulated data from 01_dgp_functions.R
# Output: Functions for computing IPW weights
# ============================================================

library(nnet)  # for multinom()

# ------------------------------------------------------------------
# estimate_propensity_scores(): Multinomial logistic regression
# ------------------------------------------------------------------
#' Estimate propensity scores using multinomial logistic regression
#'
#' Fits P(C=k|X1,X2,X3) via multinomial logit.
#' Because true C is unobserved, we use the most likely class
#' from LCA as a proxy. This is a known limitation addressed
#' in the simulation evaluation.
#'
#' @param dat       data.frame with columns X1, X2, X3, assigned_class
#' @param class_var Character name of the class variable (default "assigned_class")
#' @return A list with:
#'   - ps_matrix: N x K matrix of propensity scores
#'   - model: the fitted multinom object
#'   - converged: logical
estimate_propensity_scores <- function(dat, class_var = "assigned_class") {

  dat$C_factor <- factor(dat[[class_var]])

  # Fit multinomial logistic regression
  # Suppress convergence messages
  fit <- tryCatch(
    multinom(C_factor ~ X1 + X2 + X3, data = dat, trace = FALSE,
             maxit = 200),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(ps_matrix = NULL, model = NULL, converged = FALSE))
  }

  # Predicted probabilities P(C=k|X)
  ps_matrix <- predict(fit, type = "probs")

  # Handle edge case: if only 2 classes predicted
  if (is.null(dim(ps_matrix))) {
    ps_matrix <- cbind(1 - ps_matrix, ps_matrix)
  }

  # Ensure columns are in order
  K <- length(unique(dat$C_factor))
  if (ncol(ps_matrix) != K) {
    warning("PS matrix dimension mismatch")
  }

  list(
    ps_matrix = ps_matrix,
    model = fit,
    converged = fit$convergence == 0
  )
}

# ------------------------------------------------------------------
# compute_ipw_weights(): Inverse probability weights
# ------------------------------------------------------------------
#' Compute inverse probability weights for causal identification
#'
#' For ATE estimation with multi-valued treatment:
#'   w_i = 1 / P(C = c_i | X_i)
#'
#' Stabilized weights:
#'   w_i^{stab} = P(C = c_i) / P(C = c_i | X_i)
#'
#' @param ps_matrix  N x K matrix of propensity scores
#' @param assigned   Vector of assigned class memberships (integer 1:K)
#' @param stabilize  Logical: use stabilized weights? (default TRUE)
#' @param trim       Numeric: percentile for weight trimming (default 0.99)
#' @return Named list with:
#'   - weights: vector of IPW weights
#'   - summary: descriptive statistics of weights
compute_ipw_weights <- function(ps_matrix, assigned,
                                stabilize = TRUE, trim = 0.99) {

  N <- length(assigned)
  K <- ncol(ps_matrix)

  # Extract P(C = c_i | X_i) for each individual's assigned class
  ps_treated <- sapply(1:N, function(i) ps_matrix[i, assigned[i]])

  # Floor propensity scores to avoid extreme weights
  ps_treated <- pmax(ps_treated, 0.01)

  if (stabilize) {
    # Marginal class probabilities
    marginal <- table(assigned) / N
    numerator <- sapply(1:N, function(i) marginal[as.character(assigned[i])])
    weights <- as.numeric(numerator) / ps_treated
  } else {
    weights <- 1 / ps_treated
  }

  # Trim extreme weights at the specified percentile
  upper <- quantile(weights, trim, na.rm = TRUE)
  n_trimmed <- sum(weights > upper)
  weights <- pmin(weights, upper)

  # Normalize weights to sum to N (for interpretability)
  weights <- weights * N / sum(weights)

  list(
    weights = weights,
    n_trimmed = n_trimmed,
    summary = data.frame(
      mean = mean(weights),
      sd = sd(weights),
      min = min(weights),
      q25 = quantile(weights, 0.25),
      median = median(weights),
      q75 = quantile(weights, 0.75),
      max = max(weights),
      cv = sd(weights) / mean(weights)
    )
  )
}

# ------------------------------------------------------------------
# compute_ipw_for_true_class(): Oracle weights (for benchmarking)
# ------------------------------------------------------------------
#' Compute IPW weights using the TRUE class membership
#' This serves as a benchmark: what IPTW would achieve with
#' perfect class assignment (no measurement error)
#'
#' @param dat  data.frame with X1-X3 and true_class
#' @return IPW weight vector
compute_ipw_oracle <- function(dat) {

  ps_result <- estimate_propensity_scores(dat, class_var = "true_class")

  if (!ps_result$converged || is.null(ps_result$ps_matrix)) {
    return(rep(1, nrow(dat)))
  }

  ipw_result <- compute_ipw_weights(
    ps_matrix = ps_result$ps_matrix,
    assigned = dat$true_class,
    stabilize = TRUE, trim = 0.99
  )

  ipw_result$weights
}

# ------------------------------------------------------------------
# compute_weight_diagnostics(): Summary statistics for IPW weights
# ------------------------------------------------------------------
#' Extract weight diagnostics from an IPW result
#'
#' @param ipw_result  Output from compute_ipw_weights()
#' @return data.frame with w_mean, w_sd, w_max, w_cv, w_n_trimmed
compute_weight_diagnostics <- function(ipw_result) {
  data.frame(
    w_mean = ipw_result$summary$mean,
    w_sd = ipw_result$summary$sd,
    w_max = ipw_result$summary$max,
    w_cv = ipw_result$summary$cv,
    w_n_trimmed = ipw_result$n_trimmed
  )
}

# ------------------------------------------------------------------
# check_balance(): Covariate balance after weighting
# ------------------------------------------------------------------
#' Check covariate balance across latent classes after IPW weighting
#'
#' Computes standardized mean differences (SMD) for each covariate
#' between each class pair, before and after weighting
#'
#' @param dat       data.frame with X1-X3
#' @param assigned  class assignments
#' @param weights   IPW weights
#' @return data.frame of balance statistics
check_balance <- function(dat, assigned, weights) {

  covariates <- c("X1", "X2", "X3")
  K <- length(unique(assigned))
  results <- list()

  for (cov in covariates) {
    x <- dat[[cov]]
    # Unweighted means and SDs by class
    for (c1 in 1:(K - 1)) {
      for (c2 in (c1 + 1):K) {
        idx1 <- assigned == c1
        idx2 <- assigned == c2

        # Unweighted SMD
        m1_uw <- mean(x[idx1])
        m2_uw <- mean(x[idx2])
        s_pool <- sqrt((var(x[idx1]) + var(x[idx2])) / 2)
        smd_uw <- (m1_uw - m2_uw) / s_pool

        # Weighted SMD
        w1 <- weights[idx1]
        w2 <- weights[idx2]
        m1_w <- weighted.mean(x[idx1], w1)
        m2_w <- weighted.mean(x[idx2], w2)
        v1_w <- sum(w1 * (x[idx1] - m1_w)^2) / sum(w1)
        v2_w <- sum(w2 * (x[idx2] - m2_w)^2) / sum(w2)
        s_pool_w <- sqrt((v1_w + v2_w) / 2)
        smd_w <- (m1_w - m2_w) / s_pool_w

        results[[length(results) + 1]] <- data.frame(
          covariate = cov,
          comparison = paste0("C", c1, " vs C", c2),
          smd_unweighted = round(smd_uw, 4),
          smd_weighted = round(smd_w, 4)
        )
      }
    }
  }

  do.call(rbind, results)
}

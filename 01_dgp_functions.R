# ============================================================
# File: 01_dgp_functions.R
# Purpose: Data Generating Process (DGP) for Monte Carlo simulation
#          of IPTW-BCH causal inference with latent class exposures
# Output: Functions for generating simulated data
# ============================================================

library(MASS)

# ------------------------------------------------------------------
# Item response probability matrices for binary indicators U1-U10
# Rows = latent classes (K=3), Columns = indicators (J=10)
# ------------------------------------------------------------------

# High entropy condition (~0.8): well-separated classes
ITEM_PROBS_HIGH <- matrix(c(
  # U1   U2   U3   U4   U5   U6   U7   U8   U9   U10
  0.85, 0.80, 0.75, 0.10, 0.15, 0.10, 0.15, 0.10, 0.50, 0.45,  # Class 1
  0.15, 0.20, 0.10, 0.80, 0.85, 0.75, 0.15, 0.10, 0.50, 0.55,  # Class 2
  0.10, 0.15, 0.10, 0.15, 0.10, 0.15, 0.85, 0.80, 0.75, 0.70   # Class 3
), nrow = 3, ncol = 10, byrow = TRUE)

# Medium entropy condition (~0.7): moderate class separation
ITEM_PROBS_MEDIUM <- matrix(c(
  # U1   U2   U3   U4   U5   U6   U7   U8   U9   U10
  0.78, 0.73, 0.68, 0.18, 0.22, 0.18, 0.22, 0.18, 0.50, 0.45,  # Class 1
  0.22, 0.27, 0.18, 0.73, 0.78, 0.68, 0.22, 0.18, 0.50, 0.55,  # Class 2
  0.18, 0.22, 0.18, 0.22, 0.18, 0.22, 0.78, 0.73, 0.68, 0.63   # Class 3
), nrow = 3, ncol = 10, byrow = TRUE)

# Low entropy condition (~0.6): overlapping classes
ITEM_PROBS_LOW <- matrix(c(
  # U1   U2   U3   U4   U5   U6   U7   U8   U9   U10
  0.70, 0.65, 0.60, 0.25, 0.30, 0.25, 0.30, 0.25, 0.50, 0.45,  # Class 1
  0.30, 0.35, 0.25, 0.65, 0.70, 0.60, 0.30, 0.25, 0.50, 0.55,  # Class 2
  0.25, 0.30, 0.25, 0.30, 0.25, 0.30, 0.70, 0.65, 0.60, 0.55   # Class 3
), nrow = 3, ncol = 10, byrow = TRUE)

# ------------------------------------------------------------------
# True causal parameters
# ------------------------------------------------------------------
TRUE_ATE_C2_vs_C1 <- 0.5   # True ATE of Class 2 relative to Class 1
TRUE_ATE_C3_vs_C1 <- 1.0   # True ATE of Class 3 relative to Class 1

# Class-specific intercepts for outcome Y
ALPHA <- c(0.0, 0.5, 1.0)  # alpha_1, alpha_2, alpha_3

# Direct confounding effects of X on Y
BETA_X <- c(0.3, 0.2, 0.1) # beta_1 (X1), beta_2 (X2), beta_3 (X3)

# Residual SD of Y
SIGMA_Y <- 1.0

# ------------------------------------------------------------------
# generate_data(): Main DGP function
# ------------------------------------------------------------------
#' Generate one simulated dataset
#'
#' @param N        Sample size
#' @param entropy  Character: "high" or "low"
#' @param gamma    Confounding strength: coefficient of X on class membership
#' @param seed     Random seed (optional)
#' @return A list with:
#'   - data: data.frame with X1-X3, U1-U10, Y, true_class, ipw
#'   - true_ate: named vector of true ATEs
#'   - params: list of DGP parameters
generate_data <- function(N, entropy = "high", gamma = 0.5, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Select item response probabilities
  item_probs <- switch(entropy,
    "high"   = ITEM_PROBS_HIGH,
    "medium" = ITEM_PROBS_MEDIUM,
    "low"    = ITEM_PROBS_LOW,
    stop("entropy must be 'high', 'medium', or 'low'")
  )

  # ---- Step 1: Generate confounders X1, X2, X3 ----
  # Correlated confounders with moderate correlation
  Sigma_X <- matrix(c(
    1.0, 0.3, 0.2,
    0.3, 1.0, 0.25,
    0.2, 0.25, 1.0
  ), nrow = 3)
  X <- mvrnorm(N, mu = c(0, 0, 0), Sigma = Sigma_X)
  colnames(X) <- c("X1", "X2", "X3")

  # ---- Step 2: Generate true latent class C ----
  # Multinomial logistic model: P(C=k|X) = exp(eta_k) / sum(exp(eta_j))
  # Reference class: Class 1 (eta_1 = 0)

  # Class membership depends on confounders
  # gamma controls confounding strength
  eta1 <- rep(0, N)                                          # Reference

  eta2 <- -0.2 + gamma * X[, 1] + 0.3 * gamma * X[, 2]      # Class 2
  eta3 <- -0.5 + 0.5 * gamma * X[, 1] + gamma * X[, 3]      # Class 3

  # Softmax to get class probabilities
  eta_mat <- cbind(eta1, eta2, eta3)
  max_eta <- apply(eta_mat, 1, max)
  exp_eta <- exp(eta_mat - max_eta)  # numerical stability
  class_probs <- exp_eta / rowSums(exp_eta)

  # Draw true class membership
  true_class <- apply(class_probs, 1, function(p) {
    sample(1:3, size = 1, prob = p)
  })

  # ---- Step 3: Generate binary indicators U1-U10 ----
  U <- matrix(0, nrow = N, ncol = 10)
  colnames(U) <- paste0("U", 1:10)

  for (i in 1:N) {
    k <- true_class[i]
    U[i, ] <- rbinom(10, size = 1, prob = item_probs[k, ])
  }

  # ---- Step 4: Generate continuous outcome Y ----
  # Y = alpha_c + beta_1*X1 + beta_2*X2 + beta_3*X3 + epsilon
  # The true ATE is alpha_2 - alpha_1 = 0.5, alpha_3 - alpha_1 = 1.0
  Y <- ALPHA[true_class] +
       BETA_X[1] * X[, 1] +
       BETA_X[2] * X[, 2] +
       BETA_X[3] * X[, 3] +
       rnorm(N, 0, SIGMA_Y)

  # ---- Assemble dataset ----
  dat <- data.frame(
    X1 = X[, 1], X2 = X[, 2], X3 = X[, 3],
    U, Y = Y,
    true_class = true_class
  )

  # Store true class proportions and probabilities
  true_props <- table(true_class) / N

  list(
    data = dat,
    true_ate = c(C2_vs_C1 = TRUE_ATE_C2_vs_C1,
                 C3_vs_C1 = TRUE_ATE_C3_vs_C1),
    class_probs = class_probs,
    params = list(
      N = N, entropy = entropy, gamma = gamma,
      item_probs = item_probs,
      alpha = ALPHA, beta_x = BETA_X, sigma_y = SIGMA_Y,
      true_props = true_props
    )
  )
}

# ------------------------------------------------------------------
# compute_true_potential_outcomes(): For verification
# ------------------------------------------------------------------
#' Compute true average potential outcomes E[Y(c)] for each class
#' Under the DGP, the potential outcome Y(c) = alpha_c + beta'X + epsilon
#' So E[Y(c)] = alpha_c + beta' * E[X] = alpha_c (since E[X]=0)
#' And ATE(c vs 1) = alpha_c - alpha_1 = alpha_c
#' This holds because the confounding is fully captured by X
compute_true_potential_outcomes <- function() {
  cat("True potential outcomes (population level):\n")
  cat(sprintf("  E[Y(1)] = %.2f\n", ALPHA[1]))
  cat(sprintf("  E[Y(2)] = %.2f\n", ALPHA[2]))
  cat(sprintf("  E[Y(3)] = %.2f\n", ALPHA[3]))
  cat(sprintf("  ATE(2 vs 1) = %.2f\n", ALPHA[2] - ALPHA[1]))
  cat(sprintf("  ATE(3 vs 1) = %.2f\n", ALPHA[3] - ALPHA[1]))
}

# ============================================================
# File: 03_mplus_syntax.R
# Purpose: Generate Mplus syntax for LCA with BCH 3-step approach
#          and IPTW-BCH combined estimation
# Depends: MplusAutomation
# ============================================================

library(MplusAutomation)

# ------------------------------------------------------------------
# prepare_mplus_data(): Write data file for Mplus
# ------------------------------------------------------------------
#' Prepare data for Mplus analysis
#'
#' Writes a tab-delimited .dat file with indicators, outcome,
#' and optional IPW weights. Mplus requires specific formatting.
#'
#' @param dat        data.frame from DGP
#' @param ipw_weights  Numeric vector of IPW weights (or NULL)
#' @param file_path  Path for the .dat file
#' @return Path to the written data file
prepare_mplus_data <- function(dat, ipw_weights = NULL, file_path) {

  # Select variables for Mplus
  mplus_dat <- dat[, c(paste0("U", 1:10), "Y")]

  # Add IPW weights if provided
  if (!is.null(ipw_weights)) {
    mplus_dat$IPW <- ipw_weights
  }

  # Replace NA with -999 (Mplus missing value code)
  mplus_dat[is.na(mplus_dat)] <- -999

  # Write tab-delimited file (no header, no row names)
  write.table(mplus_dat, file = file_path,
              sep = "\t", row.names = FALSE, col.names = FALSE,
              na = "-999")

  file_path
}

# ------------------------------------------------------------------
# Step 1: Unconditional LCA with BCH weight saving
# ------------------------------------------------------------------
#' Generate Mplus syntax for Step 1: Unconditional LCA
#'
#' Estimates a K-class LCA using binary indicators U1-U10.
#' Saves BCH weights and class probabilities for Step 3.
#'
#' @param data_file   Path to the input data file
#' @param K           Number of latent classes
#' @param output_file Path for BCH weight output data
#' @param has_ipw     Logical: does data contain IPW weight column?
#' @return Character string of Mplus input syntax
generate_step1_syntax <- function(data_file, K = 3, output_file,
                                  has_ipw = FALSE) {

  var_names <- paste0("U", 1:10)
  if (has_ipw) {
    all_vars <- c(var_names, "Y", "IPW")
  } else {
    all_vars <- c(var_names, "Y")
  }

  syntax <- paste0(
    "TITLE: Step 1 - Unconditional LCA with BCH weights;\n",
    "\n",
    "DATA:\n",
    "  FILE IS ", data_file, ";\n",
    "\n",
    "VARIABLE:\n",
    "  NAMES ARE ", paste(all_vars, collapse = " "), ";\n",
    "  CATEGORICAL ARE ", paste(var_names, collapse = " "), ";\n",
    "  CLASSES = C(", K, ");\n",
    "  USEVARIABLES ARE ", paste(var_names, collapse = " "), ";\n",
    if (has_ipw) {
      paste0("  AUXILIARY = Y IPW;\n")
    } else {
      paste0("  AUXILIARY = Y;\n")
    },
    "\n",
    "ANALYSIS:\n",
    "  TYPE = MIXTURE;\n",
    "  ESTIMATOR = MLR;\n",
    "  STARTS = 100 20;\n",
    "  STITERATIONS = 20;\n",
    "  PROCESS = 2;\n",
    "\n",
    "OUTPUT:\n",
    "  TECH1 TECH8 TECH11 TECH14;\n",
    "\n",
    "SAVEDATA:\n",
    "  FILE IS ", output_file, ";\n",
    "  SAVE = BCHWEIGHTS CPROB;\n"
  )

  syntax
}

# ------------------------------------------------------------------
# Step 3a: BCH + IPTW (proposed method)
# ------------------------------------------------------------------
#' Generate Mplus syntax for Step 3: BCH with IPTW weights
#'
#' Uses BCH weights (from Step 1) as TRAINING data for the
#' latent class variable, while simultaneously applying IPTW
#' sampling weights to adjust for confounding.
#'
#' @param data_file    Path to the saved BCH weight data file
#' @param K            Number of latent classes
#' @param has_ipw      Logical: data contains IPW column
#' @return Character string of Mplus input syntax
generate_step3_bch_iptw_syntax <- function(data_file, K = 3) {

  # Variable names from SAVEDATA output (SAVE = BCHWEIGHTS CPROB):
  # U1-U10, Y, IPW, BCHW1-BCHW3, CPROB1-CPROB3, C
  bch_vars <- paste0("BCHW", 1:K)
  cprob_vars <- paste0("CPROB", 1:K)

  # Break NAMES across multiple lines to stay within Mplus 90-char limit
  names_line <- paste0(
    "  NAMES ARE U1 U2 U3 U4 U5\n",
    "  U6 U7 U8 U9 U10\n",
    "  Y IPW ", paste(bch_vars, collapse = " "), "\n",
    "  ", paste(cprob_vars, collapse = " "), " MLC;\n"
  )

  syntax <- paste0(
    "TITLE: Step 3 - BCH with IPTW;\n",
    "\n",
    "DATA:\n",
    "  FILE IS ", data_file, ";\n",
    "\n",
    "VARIABLE:\n",
    names_line,
    "  USEVARIABLES ARE Y\n",
    "  ", paste(bch_vars, collapse = " "), " IPW;\n",
    "  CLASSES = C(", K, ");\n",
    "  TRAINING = ", paste(bch_vars, collapse = " "),
    "(BCH);\n",
    "  WEIGHT = IPW;\n",
    "\n",
    "ANALYSIS:\n",
    "  TYPE = MIXTURE;\n",
    "  ESTIMATOR = MLR;\n",
    "  STARTS = 0;\n",
    "\n",
    "MODEL:\n",
    "  %OVERALL%\n",
    "  [Y];\n",
    "  Y;\n",
    "\n",
    paste(sapply(1:(K - 1), function(k) {
      paste0("  %C#", k, "%\n  [Y];\n  Y;\n")
    }), collapse = "\n"),
    "\n",
    "OUTPUT:\n",
    "  CINTERVAL;\n"
  )

  syntax
}

# ------------------------------------------------------------------
# Step 3b: BCH only (naive, no IPTW)
# ------------------------------------------------------------------
#' Generate Mplus syntax for Step 3: BCH without IPTW
#'
#' Uses BCH weights to correct for classification uncertainty
#' but does NOT apply IPTW to adjust for confounding.
#' This is the naive comparison model.
#'
#' @param data_file    Path to the saved BCH weight data file
#' @param K            Number of latent classes
#' @return Character string of Mplus input syntax
generate_step3_bch_only_syntax <- function(data_file, K = 3) {

  bch_vars <- paste0("BCHW", 1:K)
  cprob_vars <- paste0("CPROB", 1:K)

  names_line <- paste0(
    "  NAMES ARE U1 U2 U3 U4 U5\n",
    "  U6 U7 U8 U9 U10\n",
    "  Y IPW ", paste(bch_vars, collapse = " "), "\n",
    "  ", paste(cprob_vars, collapse = " "), " MLC;\n"
  )

  syntax <- paste0(
    "TITLE: Step 3 - BCH only (no IPTW);\n",
    "\n",
    "DATA:\n",
    "  FILE IS ", data_file, ";\n",
    "\n",
    "VARIABLE:\n",
    names_line,
    "  USEVARIABLES ARE Y\n",
    "  ", paste(bch_vars, collapse = " "), ";\n",
    "  CLASSES = C(", K, ");\n",
    "  TRAINING = ", paste(bch_vars, collapse = " "),
    "(BCH);\n",
    "\n",
    "ANALYSIS:\n",
    "  TYPE = MIXTURE;\n",
    "  ESTIMATOR = MLR;\n",
    "  STARTS = 0;\n",
    "\n",
    "MODEL:\n",
    "  %OVERALL%\n",
    "  [Y];\n",
    "  Y;\n",
    "\n",
    paste(sapply(1:(K - 1), function(k) {
      paste0("  %C#", k, "%\n  [Y];\n  Y;\n")
    }), collapse = "\n"),
    "\n",
    "OUTPUT:\n",
    "  CINTERVAL;\n"
  )

  syntax
}

# ------------------------------------------------------------------
# Classify-analyze (most naive approach)
# ------------------------------------------------------------------
#' Generate Mplus syntax for classify-analyze approach
#'
#' Simply uses the most likely class assignment from Step 1
#' and regresses Y on class dummies. No correction for
#' classification uncertainty or confounding.
#'
#' @param data_file  Path to saved data (with MLC column)
#' @param K          Number of latent classes
#' @return Character string of Mplus input syntax
generate_classify_analyze_syntax <- function(data_file, K = 3) {

  bch_vars <- paste0("BCHW", 1:K)
  cprob_vars <- paste0("CPROB", 1:K)

  names_line <- paste0(
    "  NAMES ARE U1 U2 U3 U4 U5\n",
    "  U6 U7 U8 U9 U10\n",
    "  Y IPW ", paste(bch_vars, collapse = " "), "\n",
    "  ", paste(cprob_vars, collapse = " "), " MLC;\n"
  )

  syntax <- paste0(
    "TITLE: Classify-analyze;\n",
    "\n",
    "DATA:\n",
    "  FILE IS ", data_file, ";\n",
    "\n",
    "VARIABLE:\n",
    names_line,
    "  USEVARIABLES ARE Y MLC;\n",
    "  CLASSES = C(", K, ");\n",
    "  KNOWNCLASS = C(MLC);\n",
    "\n",
    "ANALYSIS:\n",
    "  TYPE = MIXTURE;\n",
    "  ESTIMATOR = MLR;\n",
    "  STARTS = 0;\n",
    "\n",
    "MODEL:\n",
    "  %OVERALL%\n",
    "  [Y];\n",
    "  Y;\n",
    "\n",
    paste(sapply(1:(K - 1), function(k) {
      paste0("  %C#", k, "%\n  [Y];\n  Y;\n")
    }), collapse = "\n"),
    "\n",
    "OUTPUT:\n",
    "  CINTERVAL;\n"
  )

  syntax
}

# ------------------------------------------------------------------
# write_and_run_mplus(): Utility to write syntax and run
# ------------------------------------------------------------------
#' Write Mplus input file and run the model
#'
#' @param syntax    Mplus syntax string
#' @param inp_file  Path for .inp file
#' @param quiet     Suppress Mplus output? (default TRUE)
#' @return readModels() output or NULL on failure
write_and_run_mplus <- function(syntax, inp_file, quiet = TRUE) {

  # Write input file
  writeLines(syntax, inp_file)

  # Run Mplus
  tryCatch({
    runModels(inp_file, showOutput = !quiet, logFile = NULL)
    out_file <- sub("\\.inp$", ".out", inp_file)
    if (file.exists(out_file)) {
      readModels(out_file, quiet = quiet)
    } else {
      NULL
    }
  }, error = function(e) {
    if (!quiet) message("Mplus error: ", e$message)
    NULL
  })
}

# ------------------------------------------------------------------
# extract_class_means(): Extract class-specific means from results
# ------------------------------------------------------------------
#' Extract class-specific Y means from Mplus output
#'
#' @param results  readModels() output
#' @param K        Number of classes
#' @return data.frame with columns: class, mean, se, ci_lower, ci_upper
extract_class_means <- function(results, K = 3) {

  if (is.null(results)) return(NULL)

  params <- tryCatch(
    results$parameters$unstandardized,
    error = function(e) NULL
  )

  if (is.null(params)) return(NULL)

  # Filter for Y means (intercepts)
  y_means <- params[params$paramHeader == "Means" &
                     params$param == "Y", ]

  if (nrow(y_means) == 0) {
    # Try alternative header
    y_means <- params[params$paramHeader == "Intercepts" &
                      params$param == "Y", ]
  }

  if (nrow(y_means) == 0) return(NULL)

  # Extract class-specific values
  out <- data.frame(
    class = 1:min(nrow(y_means), K),
    mean = y_means$est[1:min(nrow(y_means), K)],
    se = y_means$se[1:min(nrow(y_means), K)]
  )

  # Compute 95% CI
  out$ci_lower <- out$mean - 1.96 * out$se
  out$ci_upper <- out$mean + 1.96 * out$se

  out
}

# ------------------------------------------------------------------
# extract_entropy(): Extract entropy from LCA Step 1
# ------------------------------------------------------------------
#' Extract entropy from Mplus LCA output
#'
#' @param results  readModels() output from Step 1
#' @return Numeric entropy value or NA
extract_entropy <- function(results) {
  if (is.null(results)) return(NA)
  tryCatch(
    results$summaries$Entropy,
    error = function(e) NA
  )
}

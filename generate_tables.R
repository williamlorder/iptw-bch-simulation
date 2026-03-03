# ============================================================
# File: generate_tables.R
# Purpose: Compute performance metrics and generate LaTeX tables
#          for the expanded 36-condition simulation
# Input:  simulation_output_v2/all_results.rds
# Output: ../tables/*.tex, simulation_output_v2/performance_metrics.rds
# ============================================================

library(dplyr)
library(tidyr)

# Source DGP constants (loads MASS before dplyr::select is used)
source("01_dgp_functions.R")

# Ensure dplyr::select takes priority over MASS::select
select <- dplyr::select

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
RESULTS_PATH <- "simulation_output_v2/all_results.rds"
TABLE_DIR <- "../tables"
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# True ATEs
TRUE_VALUES <- data.frame(
  class_pair = c("C2_vs_C1", "C3_vs_C1"),
  true_ate = c(TRUE_ATE_C2_vs_C1, TRUE_ATE_C3_vs_C1),
  stringsAsFactors = FALSE
)

# Model display names and order
MODEL_ORDER <- c("BCH_IPTW_oracle", "BCH_IPTW_est", "BCH_only", "Classify")
MODEL_LABELS_TEX <- c(
  BCH_IPTW_oracle = "BCH+IPTW (Oracle)",
  BCH_IPTW_est    = "BCH+IPTW (Est.)",
  BCH_only        = "BCH Only",
  Classify        = "Classify-Analyze"
)

ENTROPY_ORDER <- c("high", "medium", "low")
ENTROPY_LABELS <- c(high = "High", medium = "Medium", low = "Low")

# ------------------------------------------------------------------
# Load and validate results
# ------------------------------------------------------------------
cat("Loading results from:", RESULTS_PATH, "\n")
df_all <- readRDS(RESULTS_PATH)
cat(sprintf("Total rows: %d\n", nrow(df_all)))
cat(sprintf("Columns: %s\n", paste(names(df_all), collapse = ", ")))

# Convergence summary
cat(sprintf("Overall convergence: %.1f%%\n",
            100 * mean(df_all$converged, na.rm = TRUE)))

# Filter to converged
df <- df_all[df_all$converged == TRUE, ]
cat(sprintf("Converged rows: %d (%.1f%%)\n",
            nrow(df), 100 * nrow(df) / nrow(df_all)))

# ------------------------------------------------------------------
# Compute performance metrics with MCSEs
# ------------------------------------------------------------------
cat("\nComputing performance metrics...\n")

perf <- df %>%
  merge(TRUE_VALUES, by = "class_pair") %>%
  group_by(N, entropy, gamma, model, class_pair, true_ate) %>%
  summarise(
    n_reps = n(),

    # Bias
    mean_est = mean(est_ate, na.rm = TRUE),
    bias = mean(est_ate, na.rm = TRUE) - first(true_ate),
    rel_bias_pct = 100 * (mean(est_ate, na.rm = TRUE) - first(true_ate)) /
                   first(true_ate),

    # Empirical SE and model-based SE
    emp_se = sd(est_ate, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    se_ratio = mean(se, na.rm = TRUE) / sd(est_ate, na.rm = TRUE),

    # MSE and RMSE
    mse = mean((est_ate - first(true_ate))^2, na.rm = TRUE),
    rmse = sqrt(mean((est_ate - first(true_ate))^2, na.rm = TRUE)),

    # Coverage
    coverage = mean(ci_lower <= first(true_ate) &
                    ci_upper >= first(true_ate), na.rm = TRUE),
    mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),

    # Power (reject H0: ATE=0)
    power = mean(ci_lower > 0 | ci_upper < 0, na.rm = TRUE),

    # Type I error: reject when true ATE != 0
    # (meaningful only for gamma=0 in some interpretations)
    reject_rate = mean(ci_lower > 0 | ci_upper < 0, na.rm = TRUE),

    # Monte Carlo Standard Errors
    mcse_bias = sd(est_ate, na.rm = TRUE) / sqrt(n()),
    mcse_coverage = sqrt(
      mean(ci_lower <= first(true_ate) & ci_upper >= first(true_ate),
           na.rm = TRUE) *
      (1 - mean(ci_lower <= first(true_ate) & ci_upper >= first(true_ate),
                na.rm = TRUE)) / n()
    ),
    # Weight diagnostics (mean across reps for oracle/est models)
    mean_w_mean = mean(w_mean, na.rm = TRUE),
    mean_w_max = mean(w_max, na.rm = TRUE),
    mean_w_cv = mean(w_cv, na.rm = TRUE),
    mean_w_trimmed = mean(w_n_trimmed, na.rm = TRUE),

    .groups = "drop"
  )

# Factor ordering
perf$model <- factor(perf$model, levels = MODEL_ORDER)
perf$entropy <- factor(perf$entropy, levels = ENTROPY_ORDER)

# Save
saveRDS(perf, "simulation_output_v2/performance_metrics.rds")
write.csv(perf, "simulation_output_v2/performance_metrics.csv",
          row.names = FALSE)
cat("Saved performance metrics.\n")

# ------------------------------------------------------------------
# Table 2: Main results (4 estimators)
# ------------------------------------------------------------------
generate_main_table <- function(perf, cp, true_val, label_suffix,
                                output_file) {
  tab <- perf %>%
    filter(class_pair == cp) %>%
    select(N, entropy, gamma, model, bias, mcse_bias, rmse, coverage,
           mcse_coverage) %>%
    pivot_wider(
      names_from = model,
      values_from = c(bias, mcse_bias, rmse, coverage, mcse_coverage),
      names_sep = "_"
    ) %>%
    arrange(entropy, gamma, N)

  lines <- c()
  lines <- c(lines, "\\begin{table}[htbp]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, "\\singlespacing")
  lines <- c(lines, sprintf(
    "\\caption{Monte Carlo Simulation Results: $\\widehat{\\text{ATE}}(C_%s)$, True = %.1f}",
    label_suffix, true_val
  ))
  lines <- c(lines, sprintf("\\label{tab:sim_results_%s}",
                             gsub("_", "", cp)))
  lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
  lines <- c(lines, "\\footnotesize")
  lines <- c(lines, "\\begin{tabular}{llr rrr rrr rrr rrr}")
  lines <- c(lines, "\\toprule")
  lines <- c(lines, paste0(
    " & & & \\multicolumn{3}{c}{BCH+IPTW (Oracle)} & ",
    "\\multicolumn{3}{c}{BCH+IPTW (Est.)} & ",
    "\\multicolumn{3}{c}{BCH Only} & ",
    "\\multicolumn{3}{c}{Classify-Analyze} \\\\"
  ))
  lines <- c(lines, paste0(
    "\\cmidrule(lr){4-6} \\cmidrule(lr){7-9} ",
    "\\cmidrule(lr){10-12} \\cmidrule(lr){13-15}"
  ))
  lines <- c(lines, paste0(
    "Entropy & $\\gamma$ & $N$ & ",
    "Bias & RMSE & Cov. & ",
    "Bias & RMSE & Cov. & ",
    "Bias & RMSE & Cov. & ",
    "Bias & RMSE & Cov. \\\\"
  ))
  lines <- c(lines, "\\midrule")

  prev_entropy <- ""
  for (i in 1:nrow(tab)) {
    row <- tab[i, ]
    ent_label <- ENTROPY_LABELS[as.character(row$entropy)]

    if (ent_label != prev_entropy && prev_entropy != "") {
      lines <- c(lines, "\\addlinespace")
    }
    prev_entropy <- ent_label

    line <- sprintf(
      paste0(
        "%s & %.1f & %d & ",
        "$%.3f$ & $%.3f$ & $%.0f$ & ",
        "$%.3f$ & $%.3f$ & $%.0f$ & ",
        "$%.3f$ & $%.3f$ & $%.0f$ & ",
        "$%.3f$ & $%.3f$ & $%.0f$ \\\\"
      ),
      ent_label, row$gamma, row$N,
      row$bias_BCH_IPTW_oracle, row$rmse_BCH_IPTW_oracle,
      row$coverage_BCH_IPTW_oracle * 100,
      row$bias_BCH_IPTW_est, row$rmse_BCH_IPTW_est,
      row$coverage_BCH_IPTW_est * 100,
      row$bias_BCH_only, row$rmse_BCH_only,
      row$coverage_BCH_only * 100,
      row$bias_Classify, row$rmse_Classify,
      row$coverage_Classify * 100
    )
    lines <- c(lines, line)
  }

  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}%")
  lines <- c(lines, "}% end resizebox")
  lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
  lines <- c(lines, "\\begin{tablenotes}[flushleft]")
  lines <- c(lines, "\\small")
  lines <- c(lines, sprintf(paste0(
    "\\item \\textit{Note.} Bias = $\\bar{\\widehat{\\text{ATE}}} - ",
    "\\text{ATE}_{\\text{true}}$; RMSE = root mean squared error; ",
    "Cov. = 95\\%% confidence interval coverage rate (\\%%). ",
    "Results based on 500 Monte Carlo replications per condition. ",
    "BCH+IPTW (Oracle) = BCH correction with IPW from true-class ",
    "propensity scores; BCH+IPTW (Est.) = BCH correction with IPW from ",
    "modal-class propensity scores; BCH Only = BCH correction without ",
    "confounding adjustment; Classify-Analyze = modal class assignment ",
    "without corrections. True ATE = %.1f. ",
    "Entropy levels: High $\\approx$ 0.8, Medium $\\approx$ 0.7, ",
    "Low $\\approx$ 0.6."
  ), true_val))
  lines <- c(lines, "\\end{tablenotes}%")
  lines <- c(lines, "}% end resizebox")
  lines <- c(lines, "\\end{table}")

  writeLines(lines, output_file)
  cat("Wrote:", output_file, "\n")
}

generate_main_table(perf, "C2_vs_C1", 0.5, "2$ vs $C_1",
                    file.path(TABLE_DIR, "table2_main_results_c2c1.tex"))
generate_main_table(perf, "C3_vs_C1", 1.0, "3$ vs $C_1",
                    file.path(TABLE_DIR, "table2_main_results_c3c1.tex"))

# ------------------------------------------------------------------
# Table 3: SE Calibration and Power (4 estimators)
# ------------------------------------------------------------------
tab_se <- perf %>%
  filter(class_pair == "C2_vs_C1") %>%
  select(N, entropy, gamma, model, se_ratio, power) %>%
  pivot_wider(
    names_from = model,
    values_from = c(se_ratio, power),
    names_sep = "_"
  ) %>%
  arrange(entropy, gamma, N)

lines <- c()
lines <- c(lines, "\\begin{table}[htbp]")
lines <- c(lines, "\\centering")
lines <- c(lines, "\\singlespacing")
lines <- c(lines, "\\caption{Standard Error Calibration and Statistical Power}")
lines <- c(lines, "\\label{tab:se_calibration}")
lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
lines <- c(lines, "\\footnotesize")
lines <- c(lines, "\\begin{tabular}{llr rr rr rr rr}")
lines <- c(lines, "\\toprule")
lines <- c(lines, paste0(
  " & & & \\multicolumn{2}{c}{BCH+IPTW (Oracle)} & ",
  "\\multicolumn{2}{c}{BCH+IPTW (Est.)} & ",
  "\\multicolumn{2}{c}{BCH Only} & ",
  "\\multicolumn{2}{c}{Classify-Analyze} \\\\"
))
lines <- c(lines, paste0(
  "\\cmidrule(lr){4-5} \\cmidrule(lr){6-7} ",
  "\\cmidrule(lr){8-9} \\cmidrule(lr){10-11}"
))
lines <- c(lines, paste0(
  "Entropy & $\\gamma$ & $N$ & ",
  "SE Ratio & Power & SE Ratio & Power & ",
  "SE Ratio & Power & SE Ratio & Power \\\\"
))
lines <- c(lines, "\\midrule")

prev_ent <- ""
for (i in 1:nrow(tab_se)) {
  row <- tab_se[i, ]
  ent_label <- ENTROPY_LABELS[as.character(row$entropy)]
  if (ent_label != prev_ent && prev_ent != "") {
    lines <- c(lines, "\\addlinespace")
  }
  prev_ent <- ent_label

  line <- sprintf(
    paste0(
      "%s & %.1f & %d & ",
      "$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & ",
      "$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\"
    ),
    ent_label, row$gamma, row$N,
    row$se_ratio_BCH_IPTW_oracle, row$power_BCH_IPTW_oracle,
    row$se_ratio_BCH_IPTW_est, row$power_BCH_IPTW_est,
    row$se_ratio_BCH_only, row$power_BCH_only,
    row$se_ratio_Classify, row$power_Classify
  )
  lines <- c(lines, line)
}

lines <- c(lines, "\\bottomrule")
lines <- c(lines, "\\end{tabular}%")
lines <- c(lines, "}% end resizebox")
lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
lines <- c(lines, "\\begin{tablenotes}[flushleft]")
lines <- c(lines, "\\small")
lines <- c(lines, paste0(
  "\\item \\textit{Note.} SE Ratio = mean model-based SE / empirical SE ",
  "(values near 1.0 indicate well-calibrated standard errors). ",
  "Power = proportion of replications rejecting $H_0$: ATE = 0 at ",
  "$\\alpha = .05$. Results shown for ATE($C_2$ vs $C_1$) = 0.5. ",
  "500 replications per condition."
))
lines <- c(lines, "\\end{tablenotes}%")
lines <- c(lines, "}% end resizebox")
lines <- c(lines, "\\end{table}")

writeLines(lines, file.path(TABLE_DIR, "table3_se_calibration.tex"))
cat("Wrote:", file.path(TABLE_DIR, "table3_se_calibration.tex"), "\n")

# ------------------------------------------------------------------
# Table 4: Type I Error Rate (gamma = 0 conditions)
# ------------------------------------------------------------------
cat("\nComputing Type I error rates (gamma = 0)...\n")

# Under gamma=0, no confounding → BCH_only, BCH_IPTW_oracle, BCH_IPTW_est
# should all be unbiased. Type I error = rejection rate when true ATE != 0.
# But our true ATEs are 0.5 and 1.0 (never 0), so "Type I error" isn't
# directly applicable. Instead, we report coverage and rejection rates
# under the no-confounding condition to show all methods work when gamma=0.
# The key question is: does IPTW hurt when it's unnecessary?

gamma0 <- perf %>%
  filter(gamma == 0, class_pair == "C2_vs_C1") %>%
  select(N, entropy, model, bias, rmse, coverage, se_ratio,
         mcse_bias, mcse_coverage) %>%
  arrange(entropy, N, model)

# Pivot for display
tab_gamma0 <- gamma0 %>%
  pivot_wider(
    names_from = model,
    values_from = c(bias, coverage),
    names_sep = "_"
  ) %>%
  arrange(entropy, N)

lines <- c()
lines <- c(lines, "\\begin{table}[htbp]")
lines <- c(lines, "\\centering")
lines <- c(lines, "\\singlespacing")
lines <- c(lines, "\\caption{Performance Under No Confounding ($\\gamma = 0$): ATE($C_2$ vs $C_1$) = 0.5}")
lines <- c(lines, "\\label{tab:gamma0}")
lines <- c(lines, "\\resizebox{0.9\\textwidth}{!}{%")
lines <- c(lines, "\\footnotesize")
lines <- c(lines, "\\begin{tabular}{lr rr rr rr rr}")
lines <- c(lines, "\\toprule")
lines <- c(lines, paste0(
  " & & \\multicolumn{2}{c}{BCH+IPTW (Oracle)} & ",
  "\\multicolumn{2}{c}{BCH+IPTW (Est.)} & ",
  "\\multicolumn{2}{c}{BCH Only} & ",
  "\\multicolumn{2}{c}{Classify-Analyze} \\\\"
))
lines <- c(lines, paste0(
  "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} ",
  "\\cmidrule(lr){7-8} \\cmidrule(lr){9-10}"
))
lines <- c(lines, paste0(
  "Entropy & $N$ & Bias & Cov. & Bias & Cov. & ",
  "Bias & Cov. & Bias & Cov. \\\\"
))
lines <- c(lines, "\\midrule")

prev_ent <- ""
for (i in 1:nrow(tab_gamma0)) {
  row <- tab_gamma0[i, ]
  ent_label <- ENTROPY_LABELS[as.character(row$entropy)]
  if (ent_label != prev_ent && prev_ent != "") {
    lines <- c(lines, "\\addlinespace")
  }
  prev_ent <- ent_label

  line <- sprintf(
    paste0(
      "%s & %d & $%.3f$ & $%.0f$ & $%.3f$ & $%.0f$ & ",
      "$%.3f$ & $%.0f$ & $%.3f$ & $%.0f$ \\\\"
    ),
    ent_label, row$N,
    row$bias_BCH_IPTW_oracle, row$coverage_BCH_IPTW_oracle * 100,
    row$bias_BCH_IPTW_est, row$coverage_BCH_IPTW_est * 100,
    row$bias_BCH_only, row$coverage_BCH_only * 100,
    row$bias_Classify, row$coverage_Classify * 100
  )
  lines <- c(lines, line)
}

lines <- c(lines, "\\bottomrule")
lines <- c(lines, "\\end{tabular}%")
lines <- c(lines, "}% end resizebox")
lines <- c(lines, "\\resizebox{0.9\\textwidth}{!}{%")
lines <- c(lines, "\\begin{tablenotes}[flushleft]")
lines <- c(lines, "\\small")
lines <- c(lines, paste0(
  "\\item \\textit{Note.} Results under no confounding ($\\gamma = 0$). ",
  "When confounding is absent, IPTW should neither help nor harm. ",
  "Bias = $\\bar{\\widehat{\\text{ATE}}} - \\text{ATE}_{\\text{true}}$; ",
  "Cov. = 95\\% CI coverage (\\%). 500 replications per condition."
))
lines <- c(lines, "\\end{tablenotes}%")
lines <- c(lines, "}% end resizebox")
lines <- c(lines, "\\end{table}")

writeLines(lines, file.path(TABLE_DIR, "table4_gamma0.tex"))
cat("Wrote:", file.path(TABLE_DIR, "table4_gamma0.tex"), "\n")

# ------------------------------------------------------------------
# Table S1: Convergence Rates
# ------------------------------------------------------------------
cat("\nGenerating convergence table...\n")

conv_rates <- df_all %>%
  group_by(N, entropy, gamma, model) %>%
  summarise(
    total = n(),
    converged = sum(converged, na.rm = TRUE),
    conv_rate = mean(converged, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(conv_pct = sprintf("%.1f", 100 * conv_rate))

conv_tab <- conv_rates %>%
  select(N, entropy, gamma, model, conv_pct) %>%
  pivot_wider(names_from = model, values_from = conv_pct) %>%
  mutate(entropy = factor(entropy, levels = ENTROPY_ORDER)) %>%
  arrange(entropy, gamma, N)

lines <- c()
lines <- c(lines, "\\begin{table}[htbp]")
lines <- c(lines, "\\centering")
lines <- c(lines, "\\singlespacing")
lines <- c(lines, "\\caption{Convergence Rates by Simulation Condition (\\%)}")
lines <- c(lines, "\\label{tab:convergence}")
lines <- c(lines, "\\footnotesize")
lines <- c(lines, "\\begin{tabular}{llr rrrr}")
lines <- c(lines, "\\toprule")
lines <- c(lines, paste0(
  "Entropy & $\\gamma$ & $N$ & ",
  "BCH+IPTW (Or.) & BCH+IPTW (Est.) & BCH Only & Classify \\\\"
))
lines <- c(lines, "\\midrule")

prev_ent <- ""
for (i in 1:nrow(conv_tab)) {
  row <- conv_tab[i, ]
  ent_label <- ENTROPY_LABELS[as.character(row$entropy)]
  if (ent_label != prev_ent && prev_ent != "") {
    lines <- c(lines, "\\addlinespace")
  }
  prev_ent <- ent_label

  line <- sprintf(
    "%s & %.1f & %d & %s\\%% & %s\\%% & %s\\%% & %s\\%% \\\\",
    ent_label, row$gamma, row$N,
    row$BCH_IPTW_oracle, row$BCH_IPTW_est,
    row$BCH_only, row$Classify
  )
  lines <- c(lines, line)
}

lines <- c(lines, "\\bottomrule")
lines <- c(lines, "\\end{tabular}")
lines <- c(lines, "\\end{table}")

writeLines(lines, file.path(TABLE_DIR, "tableS1_convergence.tex"))
cat("Wrote:", file.path(TABLE_DIR, "tableS1_convergence.tex"), "\n")

# ------------------------------------------------------------------
# Table S2: Weight Diagnostics
# ------------------------------------------------------------------
cat("\nGenerating weight diagnostics table...\n")

weight_diag <- perf %>%
  filter(class_pair == "C2_vs_C1",
         model %in% c("BCH_IPTW_oracle", "BCH_IPTW_est")) %>%
  select(N, entropy, gamma, model, mean_w_mean, mean_w_max,
         mean_w_cv, mean_w_trimmed) %>%
  pivot_wider(
    names_from = model,
    values_from = c(mean_w_mean, mean_w_max, mean_w_cv, mean_w_trimmed),
    names_sep = "_"
  ) %>%
  arrange(entropy, gamma, N)

lines <- c()
lines <- c(lines, "\\begin{table}[htbp]")
lines <- c(lines, "\\centering")
lines <- c(lines, "\\singlespacing")
lines <- c(lines, "\\caption{IPW Weight Diagnostics by Condition}")
lines <- c(lines, "\\label{tab:weight_diag}")
lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
lines <- c(lines, "\\footnotesize")
lines <- c(lines, "\\begin{tabular}{llr rrrr rrrr}")
lines <- c(lines, "\\toprule")
lines <- c(lines, paste0(
  " & & & \\multicolumn{4}{c}{Oracle PS Weights} & ",
  "\\multicolumn{4}{c}{Estimated PS Weights} \\\\"
))
lines <- c(lines, "\\cmidrule(lr){4-7} \\cmidrule(lr){8-11}")
lines <- c(lines, paste0(
  "Entropy & $\\gamma$ & $N$ & ",
  "Mean & Max & CV & Trim & Mean & Max & CV & Trim \\\\"
))
lines <- c(lines, "\\midrule")

prev_ent <- ""
for (i in 1:nrow(weight_diag)) {
  row <- weight_diag[i, ]
  ent_label <- ENTROPY_LABELS[as.character(row$entropy)]
  if (ent_label != prev_ent && prev_ent != "") {
    lines <- c(lines, "\\addlinespace")
  }
  prev_ent <- ent_label

  line <- sprintf(
    paste0(
      "%s & %.1f & %d & ",
      "$%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & ",
      "$%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ \\\\"
    ),
    ent_label, row$gamma, row$N,
    row$mean_w_mean_BCH_IPTW_oracle, row$mean_w_max_BCH_IPTW_oracle,
    row$mean_w_cv_BCH_IPTW_oracle, row$mean_w_trimmed_BCH_IPTW_oracle,
    row$mean_w_mean_BCH_IPTW_est, row$mean_w_max_BCH_IPTW_est,
    row$mean_w_cv_BCH_IPTW_est, row$mean_w_trimmed_BCH_IPTW_est
  )
  lines <- c(lines, line)
}

lines <- c(lines, "\\bottomrule")
lines <- c(lines, "\\end{tabular}%")
lines <- c(lines, "}% end resizebox")
lines <- c(lines, "\\resizebox{\\textwidth}{!}{%")
lines <- c(lines, "\\begin{tablenotes}[flushleft]")
lines <- c(lines, "\\small")
lines <- c(lines, paste0(
  "\\item \\textit{Note.} Mean = average stabilized IPW weight; ",
  "Max = average maximum weight across replications; ",
  "CV = coefficient of variation; Trim = average number of weights ",
  "trimmed at the 99th percentile. Oracle PS = propensity scores ",
  "estimated from true class membership; Estimated PS = propensity scores ",
  "estimated from modal (most likely) class assignment. ",
  "Values averaged across 500 replications."
))
lines <- c(lines, "\\end{tablenotes}%")
lines <- c(lines, "}% end resizebox")
lines <- c(lines, "\\end{table}")

writeLines(lines, file.path(TABLE_DIR, "tableS2_weight_diagnostics.tex"))
cat("Wrote:", file.path(TABLE_DIR, "tableS2_weight_diagnostics.tex"), "\n")

# ------------------------------------------------------------------
# Console Summary
# ------------------------------------------------------------------
cat("\n====================================================\n")
cat("          SIMULATION RESULTS SUMMARY\n")
cat("====================================================\n\n")

for (cp in c("C2_vs_C1", "C3_vs_C1")) {
  true_val <- TRUE_VALUES$true_ate[TRUE_VALUES$class_pair == cp]
  cat(sprintf("--- %s (True ATE = %.1f) ---\n\n", cp, true_val))

  sub <- perf %>%
    filter(class_pair == cp) %>%
    group_by(model) %>%
    summarise(
      avg_bias = mean(bias),
      max_abs_bias = max(abs(bias)),
      avg_rmse = mean(rmse),
      avg_coverage = mean(coverage),
      avg_se_ratio = mean(se_ratio),
      .groups = "drop"
    )

  for (i in 1:nrow(sub)) {
    cat(sprintf("  %-20s Avg Bias: %+.4f  Max|Bias|: %.4f  RMSE: %.4f  Cov: %.3f  SE Ratio: %.3f\n",
                as.character(sub$model[i]), sub$avg_bias[i],
                sub$max_abs_bias[i], sub$avg_rmse[i],
                sub$avg_coverage[i], sub$avg_se_ratio[i]))
  }
  cat("\n")
}

# Key comparisons
cat("--- Key Findings ---\n\n")
cat("Under strong confounding (gamma=1.0), C2_vs_C1:\n")
strong <- perf %>% filter(gamma == 1.0, class_pair == "C2_vs_C1")
for (m in MODEL_ORDER) {
  avg <- mean(abs(strong$bias[strong$model == m]))
  cat(sprintf("  %-20s avg |bias|: %.4f\n", m, avg))
}

cat("\nOracle vs Estimated IPTW comparison (gamma=1.0, C2_vs_C1):\n")
oracle_est <- perf %>%
  filter(gamma == 1.0, class_pair == "C2_vs_C1",
         model %in% c("BCH_IPTW_oracle", "BCH_IPTW_est")) %>%
  group_by(model) %>%
  summarise(avg_bias = mean(bias), avg_rmse = mean(rmse),
            avg_cov = mean(coverage), .groups = "drop")
print(as.data.frame(oracle_est))

cat("\nAll tables generated!\n")

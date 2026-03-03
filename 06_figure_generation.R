# ============================================================
# File: 06_figure_generation.R
# Purpose: Generate publication-quality figures for manuscript
#          — 36 conditions, 4 estimators, 500 DPI
# Input:  simulation_output_v2/performance_metrics.rds
# Output: ../figures/*.png (500 DPI)
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# ------------------------------------------------------------------
# Theme and palette (journal-quality, colorblind-friendly)
# ------------------------------------------------------------------
THEME_PUB <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
  )

# 4-model palette (colorblind-friendly)
MODEL_COLORS <- c(
  "BCH_IPTW_oracle" = "#0072B2",   # Blue
  "BCH_IPTW_est"    = "#56B4E9",   # Light blue
  "BCH_only"        = "#D55E00",   # Vermillion
  "Classify"        = "#009E73"    # Bluish green
)

MODEL_LABELS <- c(
  "BCH_IPTW_oracle" = "BCH+IPTW (Oracle)",
  "BCH_IPTW_est"    = "BCH+IPTW (Estimated)",
  "BCH_only"        = "BCH Only",
  "Classify"        = "Classify-Analyze"
)

MODEL_SHAPES <- c(
  "BCH_IPTW_oracle" = 16,  # Filled circle
  "BCH_IPTW_est"    = 1,   # Open circle
  "BCH_only"        = 17,  # Filled triangle
  "Classify"        = 15   # Filled square
)

MODEL_LINETYPES <- c(
  "BCH_IPTW_oracle" = "solid",
  "BCH_IPTW_est"    = "dashed",
  "BCH_only"        = "solid",
  "Classify"        = "solid"
)

# Labels for parsed facets (use ~ for spaces in expression context)
ENTROPY_LABELS <- c(
  "high"   = "High~Entropy",
  "medium" = "Medium~Entropy",
  "low"    = "Low~Entropy"
)

MODEL_ORDER <- c("BCH_IPTW_oracle", "BCH_IPTW_est", "BCH_only", "Classify")

FIG_DIR <- "../figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Helper: prepare data for plotting
# ------------------------------------------------------------------
prep_data <- function(perf) {
  perf$model <- factor(perf$model, levels = MODEL_ORDER)
  perf$entropy_label <- factor(
    ENTROPY_LABELS[as.character(perf$entropy)],
    levels = c("High~Entropy", "Medium~Entropy", "Low~Entropy")
  )
  perf$gamma_label <- paste0("gamma == ", perf$gamma)
  perf$ate_label <- ifelse(
    perf$class_pair == "C2_vs_C1",
    "ATE(C[2]~vs~C[1]) == 0.5",
    "ATE(C[3]~vs~C[1]) == 1.0"
  )
  perf
}

# ------------------------------------------------------------------
# Figure 1: DAG
# ------------------------------------------------------------------
generate_dag_figure <- function(output_file = file.path(FIG_DIR, "Figure1_dag.tex")) {
  dag_tex <- '\\documentclass[tikz,border=10pt]{standalone}
\\usepackage{tikz}
\\usetikzlibrary{arrows.meta, positioning}

\\begin{document}
\\begin{tikzpicture}[
    node distance=2.5cm,
    every node/.style={font=\\sffamily},
    latent/.style={draw, circle, minimum size=1.2cm, thick},
    observed/.style={draw, rectangle, minimum size=1cm, thick},
    arrow/.style={-{Stealth[length=3mm]}, thick},
    dashedarrow/.style={-{Stealth[length=3mm]}, thick, dashed}
]

% Nodes
\\node[observed] (X) {$\\mathbf{X}$};
\\node[latent, right=of X] (C) {$C^*$};
\\node[observed, right=of C] (Y) {$Y$};
\\node[observed, below=1.5cm of C] (U) {$\\mathbf{U}$};

% Arrows
\\draw[arrow] (X) -- node[above, font=\\small] {confounding} (C);
\\draw[arrow] (C) -- node[above, font=\\small] {causal} (Y);
\\draw[arrow] (X) to[bend left=30] node[above, font=\\small] {direct} (Y);
\\draw[arrow] (C) -- node[right, font=\\small] {measurement} (U);

% Labels
\\node[below=0.3cm of X, font=\\footnotesize, text=gray] {Confounders};
\\node[below=0.3cm of Y, font=\\footnotesize, text=gray] {Outcome};
\\node[below=0.3cm of U, font=\\footnotesize, text=gray] {LCA Indicators};
\\node[above=0.3cm of C, font=\\footnotesize, text=gray] {Latent Class};

\\end{tikzpicture}
\\end{document}'

  writeLines(dag_tex, output_file)
  cat("Wrote DAG TikZ source:", output_file, "\n")
}

# ------------------------------------------------------------------
# Figure 2: Bias comparison (main figure)
# ------------------------------------------------------------------
plot_bias_comparison <- function(perf,
                                 output_file = file.path(FIG_DIR, "Figure2_bias.png")) {
  perf <- prep_data(perf)

  p <- ggplot(perf, aes(x = factor(N), y = bias,
                         color = model, shape = model,
                         linetype = model, group = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5), linewidth = 0.5) +
    facet_grid(ate_label ~ entropy_label + gamma_label,
               labeller = label_parsed) +
    scale_color_manual(values = MODEL_COLORS, labels = MODEL_LABELS,
                       name = "Method") +
    scale_shape_manual(values = MODEL_SHAPES, labels = MODEL_LABELS,
                       name = "Method") +
    scale_linetype_manual(values = MODEL_LINETYPES, labels = MODEL_LABELS,
                          name = "Method") +
    labs(x = "Sample Size (N)", y = "Bias",
         title = "Bias in ATE Estimation Across Simulation Conditions") +
    THEME_PUB +
    theme(legend.key.width = unit(1.5, "cm"))

  ggsave(output_file, p, width = 14, height = 7, dpi = 500)
  cat("Saved:", output_file, "\n")
}

# ------------------------------------------------------------------
# Figure 3: 95% CI Coverage
# ------------------------------------------------------------------
plot_coverage <- function(perf,
                           output_file = file.path(FIG_DIR, "Figure3_coverage.png")) {
  perf <- prep_data(perf)

  p <- ggplot(perf, aes(x = factor(N), y = coverage,
                         color = model, shape = model,
                         linetype = model, group = model)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5), linewidth = 0.5) +
    facet_grid(ate_label ~ entropy_label + gamma_label,
               labeller = label_parsed) +
    scale_color_manual(values = MODEL_COLORS, labels = MODEL_LABELS,
                       name = "Method") +
    scale_shape_manual(values = MODEL_SHAPES, labels = MODEL_LABELS,
                       name = "Method") +
    scale_linetype_manual(values = MODEL_LINETYPES, labels = MODEL_LABELS,
                          name = "Method") +
    scale_y_continuous(limits = c(0.5, 1.0),
                       breaks = seq(0.5, 1.0, 0.05)) +
    labs(x = "Sample Size (N)", y = "95% CI Coverage Rate",
         title = "Confidence Interval Coverage Across Conditions") +
    THEME_PUB +
    theme(legend.key.width = unit(1.5, "cm"))

  ggsave(output_file, p, width = 14, height = 7, dpi = 500)
  cat("Saved:", output_file, "\n")
}

# ------------------------------------------------------------------
# Figure 4: RMSE comparison
# ------------------------------------------------------------------
plot_rmse <- function(perf,
                       output_file = file.path(FIG_DIR, "Figure4_rmse.png")) {
  perf <- prep_data(perf)

  p <- ggplot(perf, aes(x = factor(N), y = rmse,
                         color = model, shape = model,
                         linetype = model, group = model)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5), linewidth = 0.5) +
    facet_grid(ate_label ~ entropy_label + gamma_label,
               labeller = label_parsed) +
    scale_color_manual(values = MODEL_COLORS, labels = MODEL_LABELS,
                       name = "Method") +
    scale_shape_manual(values = MODEL_SHAPES, labels = MODEL_LABELS,
                       name = "Method") +
    scale_linetype_manual(values = MODEL_LINETYPES, labels = MODEL_LABELS,
                          name = "Method") +
    labs(x = "Sample Size (N)", y = "RMSE",
         title = "Root Mean Squared Error Across Conditions") +
    THEME_PUB +
    theme(legend.key.width = unit(1.5, "cm"))

  ggsave(output_file, p, width = 14, height = 7, dpi = 500)
  cat("Saved:", output_file, "\n")
}

# ------------------------------------------------------------------
# Figure 5: Bias heatmap (entropy x confounding interaction)
# ------------------------------------------------------------------
plot_bias_heatmap <- function(perf,
                               output_file = file.path(FIG_DIR, "Figure5_bias_heatmap.png")) {

  heat_data <- perf %>%
    filter(class_pair == "C2_vs_C1") %>%
    group_by(entropy, gamma, model) %>%
    summarise(avg_abs_bias = mean(abs(bias)), .groups = "drop") %>%
    mutate(
      model = factor(model, levels = MODEL_ORDER),
      entropy_label = factor(
        ENTROPY_LABELS[as.character(entropy)],
        levels = c("High~Entropy", "Medium~Entropy", "Low~Entropy")
      ),
      gamma_label = paste0("gamma = ", gamma)
    )

  p <- ggplot(heat_data, aes(x = gamma_label, y = entropy_label,
                              fill = avg_abs_bias)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", avg_abs_bias)),
              color = "white", fontface = "bold", size = 3.5) +
    facet_wrap(~ model, labeller = labeller(model = MODEL_LABELS),
               nrow = 1) +
    scale_fill_gradient(low = "#2166AC", high = "#B2182B",
                        name = "Mean |Bias|") +
    labs(x = "Confounding Strength",
         y = "Classification Entropy",
         title = "Average Absolute Bias: Entropy x Confounding Interaction",
         subtitle = "ATE(C2 vs C1), averaged across sample sizes") +
    THEME_PUB +
    theme(panel.grid = element_blank())

  ggsave(output_file, p, width = 12, height = 4.5, dpi = 500)
  cat("Saved:", output_file, "\n")
}

# ------------------------------------------------------------------
# Figure S1: SE Ratio diagnostic
# ------------------------------------------------------------------
plot_se_ratio <- function(perf,
                           output_file = file.path(FIG_DIR, "FigureS1_se_ratio.png")) {
  perf <- prep_data(perf)

  p <- ggplot(perf %>% filter(class_pair == "C2_vs_C1"),
              aes(x = factor(N), y = se_ratio,
                  color = model, shape = model,
                  linetype = model, group = model)) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5), linewidth = 0.5) +
    facet_grid(gamma_label ~ entropy_label,
               labeller = label_parsed) +
    scale_color_manual(values = MODEL_COLORS, labels = MODEL_LABELS,
                       name = "Method") +
    scale_shape_manual(values = MODEL_SHAPES, labels = MODEL_LABELS,
                       name = "Method") +
    scale_linetype_manual(values = MODEL_LINETYPES, labels = MODEL_LABELS,
                          name = "Method") +
    scale_y_continuous(limits = c(0.5, 1.5)) +
    labs(x = "Sample Size (N)",
         y = "SE Ratio (Model SE / Empirical SE)",
         title = "Standard Error Calibration",
         subtitle = "Values near 1.0 indicate well-calibrated model-based SEs") +
    THEME_PUB +
    theme(legend.key.width = unit(1.5, "cm"))

  ggsave(output_file, p, width = 10, height = 8, dpi = 500)
  cat("Saved:", output_file, "\n")
}

# ------------------------------------------------------------------
# Figure S2: Weight distribution comparison
# ------------------------------------------------------------------
plot_weight_diagnostics <- function(perf,
                                     output_file = file.path(FIG_DIR, "FigureS2_weights.png")) {

  weight_data <- perf %>%
    filter(class_pair == "C2_vs_C1",
           model %in% c("BCH_IPTW_oracle", "BCH_IPTW_est")) %>%
    select(N, entropy, gamma, model, mean_w_mean, mean_w_max, mean_w_cv) %>%
    pivot_longer(cols = starts_with("mean_w_"),
                 names_to = "metric", values_to = "value") %>%
    mutate(
      model = factor(model, levels = MODEL_ORDER),
      entropy_label = factor(
        ENTROPY_LABELS[as.character(entropy)],
        levels = c("High~Entropy", "Medium~Entropy", "Low~Entropy")
      ),
      metric_label = case_when(
        metric == "mean_w_mean" ~ "Mean~Weight",
        metric == "mean_w_max"  ~ "Max~Weight",
        metric == "mean_w_cv"   ~ "CV~of~Weights"
      ),
      metric_label = factor(metric_label,
                             levels = c("Mean~Weight", "Max~Weight",
                                       "CV~of~Weights"))
    )

  p <- ggplot(weight_data,
              aes(x = factor(N), y = value,
                  color = model, shape = model, group = model)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
    geom_line(position = position_dodge(width = 0.4), linewidth = 0.5) +
    facet_grid(metric_label ~ entropy_label + paste0("gamma == ", gamma),
               scales = "free_y", labeller = label_parsed) +
    scale_color_manual(
      values = c("BCH_IPTW_oracle" = "#0072B2",
                 "BCH_IPTW_est" = "#56B4E9"),
      labels = c("BCH_IPTW_oracle" = "Oracle PS",
                 "BCH_IPTW_est" = "Estimated PS"),
      name = "PS Source"
    ) +
    scale_shape_manual(
      values = c("BCH_IPTW_oracle" = 16, "BCH_IPTW_est" = 1),
      labels = c("BCH_IPTW_oracle" = "Oracle PS",
                 "BCH_IPTW_est" = "Estimated PS"),
      name = "PS Source"
    ) +
    labs(x = "Sample Size (N)",
         y = "Value (averaged across replications)",
         title = "IPW Weight Diagnostics: Oracle vs. Estimated PS") +
    THEME_PUB

  ggsave(output_file, p, width = 14, height = 8, dpi = 500)
  cat("Saved:", output_file, "\n")
}

# ------------------------------------------------------------------
# generate_all_figures(): Master function
# ------------------------------------------------------------------
generate_all_figures <- function(perf_path = "simulation_output_v2/performance_metrics.rds") {

  if (!file.exists(perf_path)) {
    stop("Performance metrics not found. Run generate_tables.R first.")
  }

  perf <- readRDS(perf_path)

  cat("Generating figures...\n\n")

  generate_dag_figure()
  plot_bias_comparison(perf)
  plot_coverage(perf)
  plot_rmse(perf)
  plot_bias_heatmap(perf)
  plot_se_ratio(perf)
  plot_weight_diagnostics(perf)

  cat("\nAll figures saved to:", FIG_DIR, "\n")
}

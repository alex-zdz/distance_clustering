# =============================================================================
# Summary & Visualization Script – Variable Dimension Gaussian Mixture Clustering
# =============================================================================
# Purpose: Summarize simulation results and produce two key plots:
#          1. ARI vs. true labels (boxplots)
#          2. ARI vs. true labels with SALSO consensus benchmarks (Binder & VI)
# Date: March 2026
# Author: Alexander

# ---- 1. User-configurable parameters ----
chosen_D <- 10L          # ← CHANGE THIS to any dimension you ran (2, 5, 10, 20, ...)
results_base_dir <- "results/simstudy/variableD"

# ---- 2. Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(cowplot)      # for combining plots if needed
})

# ---- 3. Construct paths and load data ----
dir_D <- file.path(results_base_dir, sprintf("D%02d", chosen_D))

if (!dir.exists(dir_D)) {
  stop(sprintf("Directory not found: %s\nPlease check chosen_D and results_base_dir.", dir_D))
}

cat("Summarizing results for dimension D =", chosen_D, "from:\n", dir_D, "\n\n")

files <- list.files(dir_D, pattern = "^sim_.*\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No simulation result files (*.rds starting with 'sim_') found in the directory.")
}

cat("Found", length(files), "result files.\n")

all_results <- map_dfr(files, readRDS) |>
  mutate(
    scenario_label = case_when(
      scenario == "scen1" ~ "Scenario 1 (moderate overlap)",
      scenario == "scen2" ~ "Scenario 2 (well separated)",
      TRUE ~ scenario
    ),
    clusterwise_label = if_else(clusterwise, "Yes (per-component)", "No (global mixture)")
  )

cat("Loaded", nrow(all_results), "configurations.\n")

# ---- 4. Load posterior objects and true labels to compute consensus ARI ----
scen1_post <- readRDS(file.path(dir_D, "scen1_post.rds"))
scen2_post <- readRDS(file.path(dir_D, "scen2_post.rds"))

scen1_true <- readRDS(file.path(dir_D, "scen1_data.rds"))$cluster_true
scen2_true <- readRDS(file.path(dir_D, "scen2_data.rds"))$cluster_true

# Compute SALSO consensus once per scenario
binder_scen1 <- salso(scen1_post$clustering_matrix, loss = binder())
vi_scen1     <- salso(scen1_post$clustering_matrix, loss = VI())
ari_binder_scen1 <- adjustedRandIndex(binder_scen1, scen1_true)
ari_vi_scen1     <- adjustedRandIndex(vi_scen1,     scen1_true)

binder_scen2 <- salso(scen2_post$clustering_matrix, loss = binder())
vi_scen2     <- salso(scen2_post$clustering_matrix, loss = VI())
ari_binder_scen2 <- adjustedRandIndex(binder_scen2, scen2_true)
ari_vi_scen2     <- adjustedRandIndex(vi_scen2,     scen2_true)

# Attach to every row
all_results <- all_results |>
  mutate(
    ARI_Binder_consensus = if_else(scenario == "scen1", ari_binder_scen1, ari_binder_scen2),
    ARI_VI_consensus     = if_else(scenario == "scen1", ari_vi_scen1,     ari_vi_scen2)
  )


View(all_results)


cat("Consensus ARI vs. truth:\n",
    "  Scenario 1: Binder =", round(ari_binder_scen1, 3),
    " | VI =", round(ari_vi_scen1, 3), "\n",
    "  Scenario 2: Binder =", round(ari_binder_scen2, 3),
    " | VI =", round(ari_vi_scen2, 3), "\n\n")

# ---- 5. Prepare reference lines for plotting ----
ref_lines <- all_results |>
  distinct(scenario_label, ARI_Binder_consensus, ARI_VI_consensus) |>
  pivot_longer(
    cols = c(ARI_Binder_consensus, ARI_VI_consensus),
    names_to = "consensus_type",
    values_to = "consensus_ARI"
  ) |>
  mutate(
    consensus_type = case_when(
      consensus_type == "ARI_Binder_consensus" ~ "Binder consensus",
      consensus_type == "ARI_VI_consensus"     ~ "VI consensus"
    )
  )

# ---- 6. Plot 1 + 2 combined: ARI boxplots with consensus reference lines ----

p <- all_results |>
  ggplot(aes(x = factor(n_sweet), y = ARI, colour = clusterwise_label)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.22, dodge.width = 0.75),
             size = 2.6, alpha = 0.9) +
  
  # Dashed reference lines (same color, non-black)
  geom_hline(
    data = ref_lines,
    aes(yintercept = consensus_ARI, linetype = consensus_type),
    colour = "#636EFA",          # soft blue – change to preferred hex if desired
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  facet_wrap(~ scenario_label, ncol = 1, labeller = label_wrap_gen(multi_line = TRUE)) +
  
  scale_colour_viridis_d(option = "C", name = "Clusterwise Wasserstein") +
  scale_linetype_manual(
    values = c("Binder consensus" = "dashed", "VI consensus" = "dashed"),
    name = "Posterior consensus benchmark"
  ) +
  
  labs(
    title = sprintf("Clustering Performance – Dimension D = %d", chosen_D),
    subtitle = "Boxplots: ARI vs. true labels | Dashed lines: ARI vs. SALSO consensus (Binder & VI)",
    x = "Number of sweetening iterations (n_sweet)",
    y = "Adjusted Rand Index (ARI)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

# ---- 7. Display the plot ----
print(p)

# Optional: save high-resolution version
ggsave(sprintf("summary_D%02d_ari_with_consensus.pdf", chosen_D),
       plot = p, width = 10, height = 9, dpi = 300)
cat("\nPlot saved as: summary_D%02d_ari_with_consensus.pdf\n", chosen_D)
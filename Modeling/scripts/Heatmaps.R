# 05_Summary_Figs



library(tidyverse)
library(ggplot2)

# Assuming 'all_groups_cor_df' is already generated from your loop:
# It should contain columns: Spectral_Metric, Plant_Metric, Structural_Group, Correlation, p_val, Significance

# 1. Filter or clean metric names for better axis labels if desired
heatmap_dat <- all_groups_cor_df %>%
  mutate(
    # Add an asterisk or label for significance
    Sig_Label = case_when(
      p_val < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Optional: clean up metric names for neat plotting
    Plant_Metric = str_remove(Plant_Metric, "_TopDown|_Layered")
  )

# 2. Build the heatmap plot faceted by Structural Group
summary_heatmap <- ggplot(heatmap_dat, aes(x = Spectral_Metric, y = Plant_Metric, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # Add asterisks for significant correlations
  geom_text(aes(label = Sig_Label), color = "black", size = 4, vjust = 0.7) +
  # Diverging color scale (e.g., Cool-Warm or PiYG/PRGn)
  scale_fill_gradient2(
    low = "#21908CFF", 
    mid = "white", 
    high = "#440154FF", 
    midpoint = 0, 
    limits = c(-1, 1),
    name = "Pearson's r"
  ) +
  facet_wrap(~ Structural_Group, ncol = 3) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(
    title = "Summary of Spectral-Plant Correlations Across Structural Groups",
    subtitle = "Asterisks (*) indicate p < 0.05",
    x = "EMIT Spectral Metric",
    y = "NEON Plant Metric"
  )

print(summary_heatmap)

# 3. Save the high-res figure
ggsave("Modeling/figures/summary_correlation_heatmap.jpg", summary_heatmap, width = 14, height = 8, dpi = 300)



# -------------------------------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)

# Load data and cluster assignments
dat <- read.csv("Modeling/data_work/Joined_Results.csv")

# Re-run site clustering to get Structural_Group assignments if not already in workspace
selected_cluster_vars <- c("pct_visible_tree_cover", "mean_canopy_height_site", "visible_stems_per_ha")
site_cluster_dat <- dat %>%
  dplyr::select(siteID, all_of(selected_cluster_vars)) %>%
  drop_na()

scaled_metrics <- scale(site_cluster_dat[, selected_cluster_vars])
set.seed(42)
k_fit <- kmeans(scaled_metrics, centers = 3, nstart = 25)

site_cluster_dat <- site_cluster_dat %>%
  mutate(
    Structural_Group = case_when(
      visible_stems_per_ha > 600 ~ "High Stem Density",
      pct_visible_tree_cover < 40 ~ "Sparse Cover",
      TRUE ~ "Dense Canopy"
    ),
    Structural_Group = factor(
      Structural_Group, 
      levels = c("Sparse Cover", "Dense Canopy", "High Stem Density")
    )
  )

dat_analysis <- dat %>%
  left_join(site_cluster_dat %>% dplyr::select(siteID, Structural_Group), by = "siteID")

spec_vars_use <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise",  
                   "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise")

topdown_vars <- c("TD_Richness_TopDown", "TD_Shannon_Eff_TopDown", "TD_Simpson_Eff_TopDown", "TD_Beta_Bray_TopDown",  
                  "TD_Beta_Dispersion_TopDown", "PD_Faith_TopDown", "PD_Beta_Sorensen_TopDown", "PD_Beta_Dispersion_TopDown",
                  "FD_FDis_m1_TopDown", "FD_RaoQ_m1_TopDown", "FD_Beta_Dispersion_m1_TopDown", "FD_Beta_Turnover_m1_TopDown")

layered_vars <- c("TD_Richness_Layered", "TD_Shannon_Eff_Layered", "TD_Simpson_Eff_Layered", "TD_Beta_Bray_Layered",
                  "TD_Beta_Dispersion_Layered", "PD_Faith_Layered", "PD_Beta_Sorensen_Layered", "PD_Beta_Dispersion_Layered",
                  "FD_FDis_m1_Layered", "FD_RaoQ_m1_Layered", "FD_Beta_Dispersion_m1_Layered", "FD_Beta_Turnover_m1_Layered")

# Function to compute correlation matrices for a given subset
get_cor_df <- function(df, group_label) {
  cor_res <- psych::corr.test(df[, spec_vars_use], df[, topdown_vars], method = "pearson", adjust = "none")
  
  as.data.frame(as.table(cor_res$r)) %>%
    rename(Spectral_Metric = Var1, Plant_Metric = Var2, Correlation = Freq) %>%
    left_join(
      as.data.frame(as.table(cor_res$p)) %>% rename(Spectral_Metric = Var1, Plant_Metric = Var2, p_val = Freq),
      by = c("Spectral_Metric", "Plant_Metric")
    ) %>%
    mutate(
      Dataset = group_label,
      Significance = if_else(p_val < 0.05, "*", "")
    )
}

# 1. All Sites correlation
all_sites_df <- get_cor_df(dat_analysis, "All Sites")

# 2. Group-specific correlations
group_df <- dat_analysis %>%
  drop_na(Structural_Group) %>%
  filter(Structural_Group != "High Stem Density") %>%
  group_split(Structural_Group) %>%
  map_dfr(~ get_cor_df(.x, unique(.x$Structural_Group)))

# Combine all into one plotting dataframe
heatmap_data_topdown <- bind_rows(all_sites_df, group_df) %>%
  mutate(Dataset = factor(Dataset, levels = c("All Sites", "Sparse Cover", "Dense Canopy", "High Stem Density")))

# Plot Heatmap
topdown_heatmap <- ggplot(heatmap_data_topdown, aes(x = Spectral_Metric, y = Plant_Metric, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Significance), color = "black", size = 4, vjust = 0.7) +
  scale_fill_gradient2(low = "#21908CFF", mid = "white", high = "#440154FF", midpoint = 0, limit = c(-1, 1), name = "Pearson (r)") +
  facet_wrap(~ Dataset, ncol = 3) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "TopDown NEON vs. EMIT Spectral Metric Correlations",
    subtitle = "(* p < 0.05)",
    x = "EMIT Spectral Metric",
    y = "NEON TopDown Metric"
  )

print(topdown_heatmap)
ggsave("Modeling/figures/TopDown_correlation_heatmap.jpg", topdown_heatmap, width = 12, height = 10, dpi = 300)

# 04_woodyCanopyCover: site data plotting

library(tidyverse)
library(ggplot2)
library(psych)
library(ggpubr)
library(cluster)
library(ggrepel)
library(plotly)

# Load data 
dat <- read.csv("Modeling/data_work/Joined_Results.csv")

spec_vars <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise",
               "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise",
               "plot_sum_squares", "plot_beta_dispersion", "plot_avg_pairwise")

site_vars <- c("n_plots", "total_visible_stems", "visible_stems_per_ha", 
               "avg_visible_canopy_per_plot", "pct_visible_tree_cover",
               "mean_canopy_height_site", "max_canopy_height_site", 
               "mean_base_crown_height_site", "mean_elevation")


# -------------------------------------------------------------------------------------
# Global corr test
global_cor_res <- psych::corr.test(
  dat[, spec_vars], 
  dat[, site_vars], 
  method = "pearson", 
  adjust = "none"
)

# -------------------------------------------------------------------------------------
# gt sig
all_pairs_long <- as.data.frame(as.table(global_cor_res$r)) %>%
  rename(Spectral_Metric = Var1, Plant_Metric = Var2, Correlation = Freq) %>%
  left_join(
    as.data.frame(as.table(global_cor_res$p)) %>%
      rename(Spectral_Metric = Var1, Plant_Metric = Var2, p_val = Freq),
    by = c("Spectral_Metric", "Plant_Metric")
  )

significant_pairs <- all_pairs_long %>% filter(p_val < 0.05)

# -------------------------------------------------------------------------------------
# plot em
plot_ready_all <- dat %>%
  dplyr::select(siteID, pct_visible_tree_cover, all_of(spec_vars), all_of(site_vars)) %>%
  pivot_longer(cols = all_of(spec_vars), names_to = "Spectral_Metric", values_to = "Spectral_Value") %>%
  pivot_longer(cols = all_of(site_vars), names_to = "Plant_Metric", values_to = "Plant_Value")

plot_ready_all_merged <- plot_ready_all %>%
  left_join(all_pairs_long, by = c("Spectral_Metric", "Plant_Metric"))

# full gallary
all_pairs_plot <- ggplot(plot_ready_all_merged, aes(x = Spectral_Value, y = Plant_Value)) +
  geom_point(color = "black", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "grey90", linewidth = 0.6, linetype = "solid") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           color = "black", size = 3, label.x.npc = "left") +
  facet_wrap(Plant_Metric ~ Spectral_Metric, scales = "free", ncol = 8) +
  scale_color_viridis_c(option = "mako", direction = -1) + 
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right", 
    legend.key.height = unit(1, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  ) +
  labs(
    title = "All EMIT vs. Top-Down NEON Pairs",
    x = "EMIT Spectral Metric Value",
    y = "NEON Top-Down Metric Value",
    color = "% Tree Cover"
  )

#print(all_pairs_plot)

# plot sigs
plot_ready_sig <- plot_ready_all %>%
  inner_join(significant_pairs, by = c("Spectral_Metric", "Plant_Metric"))

num_sig_plots <- nrow(significant_pairs)
calculated_cols <- if(num_sig_plots <= 3) num_sig_plots else 3

sig_topdown_plot <- ggplot(plot_ready_sig, aes(x = Spectral_Value, y = Plant_Value)) +
  geom_point(color = "gray", alpha = 0.8, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "grey90", linewidth = 0.8, linetype = "solid") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           color = "black", size = 3.5, label.x.npc = "left") +
  facet_wrap(Plant_Metric ~ Spectral_Metric, scales = "free", ncol = calculated_cols) +
  scale_color_viridis_c(option = "mako", direction = -1) + 
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right", 
    legend.key.height = unit(1, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  ) +
  labs(
    title = "Significant EMIT vs. Top-Down NEON Metrics",
    subtitle = paste("Pearson Correlation | p < 0.05 | Total Significant Pairs:", num_sig_plots),
    x = "EMIT Spectral Metric Value",
    y = "NEON Top-Down Metric Value",
    color = "% Tree Cover"
  )

print(sig_topdown_plot)
ggsave("Modeling/figures/significant_woody_pairs.jpg", sig_topdown_plot, width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------------
# site variable correlations and choice
site_cor_res <- psych::corr.test(
  dat[, site_vars], 
  method = "pearson", 
  adjust = "holm"
)

print(round(site_cor_res$r, 2))

corrplot::corrplot(
  site_cor_res$r, 
  method = "color", 
  type = "upper", 
  order = "hclust", 
  addCoef.col = "black", 
  tl.col = "black", 
  tl.srt = 45, 
  sig.level = 0.05, 
  insig = "blank",
  diag = FALSE,
  title = "Correlation Matrix of Site Structural & Topographic Variables",
  mar = c(0,0,1,0)
)

selected_cluster_vars <- c(
  "pct_visible_tree_cover", 
  "mean_canopy_height_site", 
  "visible_stems_per_ha"
)

# -------------------------------------------------------------------------------------
# look at clustering based on woody canopy cover
site_cluster_dat <- dat %>%
  dplyr::select(siteID, all_of(selected_cluster_vars)) %>%
  drop_na()

scaled_metrics <- scale(site_cluster_dat[, selected_cluster_vars])

# k means clustering for 3
set.seed(42)
k_fit <- kmeans(scaled_metrics, centers = 3, nstart = 25)

# 4. Assign clusters and structural group labels
site_cluster_dat <- site_cluster_dat %>%
  mutate(
    Cluster_ID = factor(k_fit$cluster),
    Structural_Group = case_when(
      visible_stems_per_ha > 600 ~ "High Stem Density",
      pct_visible_tree_cover < 40 ~ "Open / Sparse Cover",
      TRUE ~ "Closed / Dense Canopy"
    ),
    Structural_Group = factor(
      Structural_Group, 
      levels = c("Open / Sparse Cover", "Closed / Dense Canopy", "High Stem Density")
    )
  )

# 3d plot
p_3d <- plot_ly(
  site_cluster_dat,
  x = ~pct_visible_tree_cover,
  y = ~visible_stems_per_ha,
  z = ~mean_canopy_height_site,
  color = ~Structural_Group,
  colors = c("#440154FF", "#21908CFF", "#FDE725FF"),
  text = ~siteID,
  type = "scatter3d",
  mode = "text+markers",
  marker = list(size = 5)
) %>%
  plotly::layout(
    title = "3D Site Structural Clustering (k = 3)",
    scene = list(
      xaxis = list(title = "Percent Tree Cover (%)"),
      yaxis = list(title = "Stems per Hectare"),
      zaxis = list(title = "Mean Canopy Height (m)")
    )
  )

print(p_3d)

# -------------------------------------------------------------------------------------
# run layered and topdown analysis with these groups 
# TopDown
dat_analysis <- dat %>%
  left_join(
    site_cluster_dat %>% dplyr::select(siteID, Structural_Group), 
    by = "siteID"
  ) %>%
  drop_na(Structural_Group)

spec_vars_use <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise",  
                   "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise")

topdown_vars <- c("TD_Richness_TopDown", "TD_Shannon_Eff_TopDown", "TD_Simpson_Eff_TopDown", "TD_Beta_Bray_TopDown",  
                  "TD_Beta_Dispersion_TopDown", "PD_Faith_TopDown", "PD_Beta_Sorensen_TopDown", "PD_Beta_Dispersion_TopDown",
                  "FD_FDis_m1_TopDown", "FD_RaoQ_m1_TopDown", "FD_Beta_Dispersion_m1_TopDown", "FD_Beta_Turnover_m1_TopDown",
                  "FD_FDis_m2_TopDown", "FD_RaoQ_m2_TopDown", "FD_Beta_Dispersion_m2_TopDown", "FD_Beta_Turnover_m2_TopDown")


# group correlations
run_group_correlations <- function(df, group_name) {
  subset_df <- df %>% filter(Structural_Group == group_name)
  
  # Correlation test
  cor_res <- psych::corr.test(
    subset_df[, spec_vars_use], 
    subset_df[, topdown_vars], 
    method = "pearson", 
    adjust = "none"
  )
  
  # Format to long table
  as.data.frame(as.table(cor_res$r)) %>%
    rename(Spectral_Metric = Var1, Plant_Metric = Var2, Correlation = Freq) %>%
    left_join(
      as.data.frame(as.table(cor_res$p)) %>%
        rename(Spectral_Metric = Var1, Plant_Metric = Var2, p_val = Freq),
      by = c("Spectral_Metric", "Plant_Metric")
    ) %>%
    mutate(
      Structural_Group = group_name,
      R_squared = Correlation^2,
      Significance = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
}

groups <- unique(dat_analysis$Structural_Group)
all_groups_cor_list <- map(groups, ~ run_group_correlations(dat_analysis, .x))
all_groups_cor_df <- bind_rows(all_groups_cor_list)


# sig plots
plot_ready_all <- dat_analysis %>%
  dplyr::select(siteID, Structural_Group, pct_visible_tree_cover, all_of(spec_vars_use), all_of(topdown_vars)) %>%
  pivot_longer(cols = all_of(spec_vars_use), names_to = "Spectral_Metric", values_to = "Spectral_Value") %>%
  pivot_longer(cols = all_of(topdown_vars), names_to = "Plant_Metric", values_to = "Plant_Value")

plot_ready_merged <- plot_ready_all %>%
  inner_join(all_groups_cor_df %>% filter(p_val < 0.05), by = c("Structural_Group", "Spectral_Metric", "Plant_Metric"))

# facet sig plot
sig_clustered_topdown_plot <- ggplot(plot_ready_merged, aes(x = Spectral_Value, y = Plant_Value, color = Structural_Group)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "grey90", linewidth = 0.6, linetype = "solid") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           color = "black", size = 3, label.x.npc = "left") +
  facet_wrap(Structural_Group ~ Plant_Metric + Spectral_Metric, scales = "free", ncol = 4) +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) + 
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom", 
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "EMIT vs. Top-Down NEON Metrics Split by Structural Group",
    subtitle = "Significant Correlations (p < 0.05) evaluated within each structural cluster",
    x = "EMIT Spectral Metric Value",
    y = "NEON Top-Down Metric Value",
    color = "Structural Group"
  )

print(sig_clustered_topdown_plot)
ggsave("Modeling/figures/significant_topdown_by_structural_group.jpg", sig_clustered_topdown_plot, width = 14, height = 12, dpi = 300)


# group-specific tables
topdown_taxonomic <- topdown_vars[str_detect(topdown_vars, "^TD_")]
topdown_phylogenetic <- topdown_vars[str_detect(topdown_vars, "^PD_")]
topdown_functional <- topdown_vars[str_detect(topdown_vars, "^FD_")]

write.csv(all_groups_cor_df %>% filter(Plant_Metric %in% topdown_taxonomic), 
          "Modeling/data_out/topDown_taxonomic_correlations_by_group.csv", row.names = FALSE)
write.csv(all_groups_cor_df %>% filter(Plant_Metric %in% topdown_phylogenetic), 
          "Modeling/data_out/topDown_phylogenetic_correlations_by_group.csv", row.names = FALSE)
write.csv(all_groups_cor_df %>% filter(Plant_Metric %in% topdown_functional), 
          "Modeling/data_out/topDown_functional_correlations_by_group.csv", row.names = FALSE)

# -------------------------------------------------------------------------------------
# run layered and topdown analysis with these groups 
# Layered
dat_analysis <- dat %>%
  left_join(
    site_cluster_dat %>% dplyr::select(siteID, Structural_Group), 
    by = "siteID"
  ) %>%
  drop_na(Structural_Group)

spec_vars_use <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise",  
                   "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise")

layered_vars <- c("TD_Richness_Layered", "TD_Shannon_Eff_Layered", "TD_Simpson_Eff_Layered", "TD_Beta_Bray_Layered",
                  "TD_Beta_Dispersion_Layered", "PD_Faith_Layered", "PD_Beta_Sorensen_Layered", "PD_Beta_Dispersion_Layered",
                  "FD_FDis_m1_Layered", "FD_RaoQ_m1_Layered", "FD_Beta_Dispersion_m1_Layered", "FD_Beta_Turnover_m1_Layered",
                  "FD_FDis_m2_Layered", "FD_RaoQ_m2_Layered", "FD_Beta_Dispersion_m2_Layered", "FD_Beta_Turnover_m2_Layered")

# group correlations
run_group_correlations <- function(df, group_name) {
  subset_df <- df %>% filter(Structural_Group == group_name)
  
  # Correlation test
  cor_res <- psych::corr.test(
    subset_df[, spec_vars_use], 
    subset_df[, layered_vars], 
    method = "pearson", 
    adjust = "none"
  )
  
  # Format to long table
  as.data.frame(as.table(cor_res$r)) %>%
    rename(Spectral_Metric = Var1, Plant_Metric = Var2, Correlation = Freq) %>%
    left_join(
      as.data.frame(as.table(cor_res$p)) %>%
        rename(Spectral_Metric = Var1, Plant_Metric = Var2, p_val = Freq),
      by = c("Spectral_Metric", "Plant_Metric")
    ) %>%
    mutate(
      Structural_Group = group_name,
      R_squared = Correlation^2,
      Significance = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
}

groups <- unique(dat_analysis$Structural_Group)
all_groups_cor_list <- map(groups, ~ run_group_correlations(dat_analysis, .x))
all_groups_cor_df <- bind_rows(all_groups_cor_list)


# sig plots
plot_ready_all <- dat_analysis %>%
  dplyr::select(siteID, Structural_Group, pct_visible_tree_cover, all_of(spec_vars_use), all_of(layered_vars)) %>%
  pivot_longer(cols = all_of(spec_vars_use), names_to = "Spectral_Metric", values_to = "Spectral_Value") %>%
  pivot_longer(cols = all_of(layered_vars), names_to = "Plant_Metric", values_to = "Plant_Value")

plot_ready_merged <- plot_ready_all %>%
  inner_join(all_groups_cor_df %>% filter(p_val < 0.05), by = c("Structural_Group", "Spectral_Metric", "Plant_Metric"))

# facet sig plot
sig_clustered_topdown_plot <- ggplot(plot_ready_merged, aes(x = Spectral_Value, y = Plant_Value, color = Structural_Group)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "grey90", linewidth = 0.6, linetype = "solid") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           color = "black", size = 3, label.x.npc = "left") +
  facet_wrap(Structural_Group ~ Plant_Metric + Spectral_Metric, scales = "free", ncol = 4) +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) + 
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom", 
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "EMIT vs. Layered NEON Metrics Split by Structural Group",
    subtitle = "Significant Correlations (p < 0.05) evaluated within each structural cluster",
    x = "EMIT Spectral Metric Value",
    y = "NEON Top-Down Metric Value",
    color = "Structural Group"
  )

print(sig_clustered_topdown_plot)
ggsave("Modeling/figures/significant_layered_by_structural_group.jpg", sig_clustered_topdown_plot, width = 14, height = 12, dpi = 300)


# group-specific tables
layered_taxonomic <- layered_vars[str_detect(layered_vars, "^TD_")]
layered_phylogenetic <- layered_vars[str_detect(layered_vars, "^PD_")]
layered_functional <- layered_vars[str_detect(layered_vars, "^FD_")]

write.csv(all_groups_cor_df %>% filter(Plant_Metric %in% layered_taxonomic), 
          "Modeling/data_out/layered_taxonomic_correlations_by_group.csv", row.names = FALSE)
write.csv(all_groups_cor_df %>% filter(Plant_Metric %in% layered_phylogenetic), 
          "Modeling/data_out/layered_phylogenetic_correlations_by_group.csv", row.names = FALSE)
write.csv(all_groups_cor_df %>% filter(Plant_Metric %in% layered_functional), 
          "Modeling/data_out/layered_functional_correlations_by_group.csv", row.names = FALSE)


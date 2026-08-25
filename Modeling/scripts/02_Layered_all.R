# 02_Layered-Analysis

library(tidyverse)
library(ggplot2)
library(psych)
library(ggpubr)

# load
dat <- read.csv("Modeling/data_work/Joined_Results.csv")

# Define spectral variables
spec_vars <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise", "gamma_sum_squares",
               "gamma_dispersion", "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise" ,
               "plot_sum_squares", "plot_beta_dispersion", "plot_avg_pairwise")

# Define top-down plant variables
layered_vars <- c("TD_Richness_Layered", "TD_Shannon_Eff_Layered", "TD_Simpson_Eff_Layered", "TD_Beta_Bray_Layered",
                  "TD_Beta_Dispersion_Layered", "PD_Faith_Layered", "PD_Beta_Sorensen_Layered", "PD_Beta_Dispersion_Layered",
                  "FD_FDis_m1_Layered", "FD_RaoQ_m1_Layered", "FD_Beta_Dispersion_m1_Layered", "FD_Beta_Turnover_m1_Layered",
                  "FD_FDis_m2_Layered", "FD_RaoQ_m2_Layered", "FD_Beta_Dispersion_m2_Layered", "FD_Beta_Turnover_m2_Layered")

# -------------------------------------------------------------------------------------
# Global corr test
global_cor_res <- psych::corr.test(
  dat[, spec_vars_use], 
  dat[, layered_vars], 
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
  dplyr::select(siteID, pct_visible_tree_cover, all_of(spec_vars_use), all_of(layered_vars)) %>%
  pivot_longer(cols = all_of(spec_vars_use), names_to = "Spectral_Metric", values_to = "Spectral_Value") %>%
  pivot_longer(cols = all_of(layered_vars), names_to = "Plant_Metric", values_to = "Plant_Value")

plot_ready_all_merged <- plot_ready_all %>%
  left_join(all_pairs_long, by = c("Spectral_Metric", "Plant_Metric"))

# full gallary
all_pairs_plot <- ggplot(plot_ready_all_merged, aes(x = Spectral_Value, y = Plant_Value)) +
  geom_point(aes(color = pct_visible_tree_cover), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "grey90", linewidth = 0.6, linetype = "solid") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           color = "black", size = 3, label.x.npc = "left") +
  facet_wrap(Plant_Metric ~ Spectral_Metric, scales = "free", ncol = 16) +
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

print(all_pairs_plot)

# plot sigs
plot_ready_sig <- plot_ready_all %>%
  inner_join(significant_pairs, by = c("Spectral_Metric", "Plant_Metric"))

num_sig_plots <- nrow(significant_pairs)
calculated_cols <- if(num_sig_plots <= 3) num_sig_plots else 3

sig_topdown_plot <- ggplot(plot_ready_sig, aes(x = Spectral_Value, y = Plant_Value)) +
  geom_point(aes(color = pct_visible_tree_cover), alpha = 0.8, size = 2.5) +
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
ggsave("Modeling/figures/significant_layered_pairs.jpg", sig_topdown_plot, width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------------
# table of correlations w/ pvalues
layered_taxonomic <- layered_vars[str_detect(layered_vars, "^TD_")]
layered_phylogenetic <- layered_vars[str_detect(layered_vars, "^PD_")]
layered_functional <- layered_vars[str_detect(layered_vars, "^FD_")]

# add R2 column
all_pairs_table <- all_pairs_long %>%
  mutate(
    R_squared = Correlation^2,
    Significance = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  ) %>%
  arrange(p_val, desc(abs(Correlation)))

# filter each table
taxonomic_table <- all_pairs_table %>% 
  filter(Plant_Metric %in% layered_taxonomic)

phylogenetic_table <- all_pairs_table %>% 
  filter(Plant_Metric %in% layered_phylogenetic)

functional_table <- all_pairs_table %>% 
  filter(Plant_Metric %in% layered_functional)

write.csv(taxonomic_table, "Modeling/data_out/layered_taxonomic_correlations.csv", row.names = FALSE)
write.csv(phylogenetic_table, "Modeling/data_out/layered_phylogenetic_correlations.csv", row.names = FALSE)
write.csv(functional_table, "Modeling/data_out/layered_functional_correlations.csv", row.names = FALSE)


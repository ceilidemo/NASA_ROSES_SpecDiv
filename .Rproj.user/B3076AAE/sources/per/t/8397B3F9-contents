# 01_variableEDA

library(tidyverse)
library(ggplot2)

# -------------------------------------------------------------------------------------
# load results final 
dat <- read.csv("Modeling/data_work/Joined_Results.csv")

# group variables...
spec_vars <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise", "gamma_sum_squares",
  "gamma_dispersion", "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise" )

layered_vars <- c("TD_Richness_Layered", "TD_Shannon_Eff_Layered", "TD_Simpson_Eff_Layered", "TD_Beta_Bray_Layered",
  "TD_Beta_Dispersion_Layered", "PD_Faith_Layered", "PD_Beta_Sorensen_Layered", "PD_Beta_Dispersion_Layered",
  "FD_FDis_m1_Layered", "FD_RaoQ_m1_Layered", "FD_Beta_Dispersion_m1_Layered", "FD_Beta_Turnover_m1_Layered",
  "FD_FDis_m2_Layered", "FD_RaoQ_m2_Layered", "FD_Beta_Dispersion_m2_Layered", "FD_Beta_Turnover_m2_Layered")

topdown_vars <- c("TD_Richness_TopDown", "TD_Shannon_Eff_TopDown", "TD_Simpson_Eff_TopDown", "TD_Beta_Bray_TopDown", 
         "TD_Beta_Dispersion_TopDown", "PD_Faith_TopDown", "PD_Beta_Sorensen_TopDown", "PD_Beta_Dispersion_TopDown",
         "FD_FDis_m1_TopDown", "FD_RaoQ_m1_TopDown", "FD_Beta_Dispersion_m1_TopDown", "FD_Beta_Turnover_m1_TopDown",
         "FD_FDis_m2_TopDown", "FD_RaoQ_m2_TopDown", "FD_Beta_Dispersion_m2_TopDown", "FD_Beta_Turnover_m2_TopDown")

site_vars <- c("n_plots", "total_visible_stems", "visible_stems_per_ha", 
               "avg_visible_canopy_per_plot", "pct_visible_tree_cover",
               "mean_canopy_height_site", "max_canopy_height_site", 
               "mean_base_crown_height_site", "mean_elevation")
# -------------------------------------------------------------------------------------
# spec variable correlations
spec_cor_matrix <- cor(dat %>% select(all_of(spec_vars)), use = "pairwise.complete.obs")
spec_cor_long <- as.data.frame(as.table(spec_cor_matrix)) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation = Freq)

spec_vars_use <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise", 
                       "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise")

# -------------------------------------------------------------------------------------
# layerd vs topdown plant
base_metrics <- str_remove(layered_vars, "_Layered")

matching_corrs <- map_dfr(base_metrics, function(m) {
  col_layer <- paste0(m, "_Layered")
  col_top   <- paste0(m, "_TopDown")
  
  cor_val <- cor(dat[[col_layer]], dat[[col_top]], use = "pairwise.complete.obs")
  
  tibble(
    Metric = m,
    Correlation = cor_val
  )
})

print(matching_corrs)

ggplot(matching_corrs, aes(x = reorder(Metric, Correlation), y = Correlation)) +
  geom_col(fill = "steelblue", width = 0.6) +
  coord_flip() +
  ylim(0, 1) +
  labs(
    title = "Correlations between Layered and TopDown Plant Metrics",
    x = "Plant Metric",
    y = "Pearson Correlation (r)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------------------------------------------------------------------------------------
# layered plant vs spec
plant_spec_layered_matrix <- cor(
  dat %>% select(all_of(spec_vars_use)), 
  dat %>% select(all_of(layered_vars)), 
  use = "pairwise.complete.obs"
)

plant_spec_layered_long <- as.data.frame(as.table(plant_spec_layered_matrix)) %>%
  rename(Spectral_Var = Var1, Plant_Metric = Var2, Correlation = Freq) %>%
  mutate(Approach = "Layered")

# -------------------------------------------------------------------------------------
# topdown plant vs spec
plant_spec_topdown_matrix <- cor(
  dat %>% select(all_of(spec_vars_use)), 
  dat %>% select(all_of(topdown_vars)), 
  use = "pairwise.complete.obs"
)

plant_spec_topdown_long <- as.data.frame(as.table(plant_spec_topdown_matrix)) %>%
  rename(Spectral_Var = Var1, Plant_Metric = Var2, Correlation = Freq) %>%
  mutate(Approach = "TopDown")

# -------------------------------------------------------------------------------------
# all plant vs spec
plant_spec_combined <- bind_rows(plant_spec_layered_long, plant_spec_topdown_long)

plant_spec_combined <- plant_spec_combined %>%
  mutate(Base_Metric = str_remove(Plant_Metric, "_Layered|_TopDown"))

print(head(plant_spec_combined %>% arrange(desc(abs(Correlation))), 15))





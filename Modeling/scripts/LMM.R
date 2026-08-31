library(tidyverse)
library(ggplot2)
library(ggeffects)
library(broom)
library(ggpubr)

# Scale and prep data
dat_model <- dat %>%
  mutate(
    scaled_tree_cover = as.numeric(scale(pct_visible_tree_cover)),
    scaled_canopy_height = as.numeric(scale(mean_canopy_height_site)),
    scaled_stems = as.numeric(scale(visible_stems_per_ha)),
    scaled_sum_squares = as.numeric(scale(sum_squares))
  )

spec_vars_use <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise",
                   "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise")

layered_vars <- c("PD_Faith_Layered", "PD_Beta_Sorensen_Layered", "PD_Beta_Dispersion_Layered",
                  "FD_FDis_m1_Layered", "FD_RaoQ_m1_Layered", "FD_Beta_Dispersion_m1_Layered", "FD_Beta_Turnover_m1_Layered")

topdown_vars <- c("PD_Faith_TopDown", "PD_Beta_Sorensen_TopDown", "PD_Beta_Dispersion_TopDown",
                  "FD_FDis_m1_TopDown", "FD_RaoQ_m1_TopDown", "FD_Beta_Dispersion_m1_TopDown", "FD_Beta_Turnover_m1_TopDown")


# Fit models and pull sum stats
run_model_summary <- function(spec, plant) {
  formula_str <- paste0(plant, " ~ ", spec, " * scaled_tree_cover + scaled_canopy_height")
  
  model <- tryCatch({
    lm(as.formula(formula_str), data = dat_model)
  }, error = function(e) { NULL })
  
  if (is.null(model)) return(NULL)
  
  glance_res <- broom::glance(model)
  model_p_val <- glance_res$p.value
  
  if (!is.na(model_p_val) && model_p_val < 0.05) {
    tibble(
      Plant_Metric = plant,
      Spectral_Metric = spec,
      R_squared = round(glance_res$r.squared, 3),
      Adj_R_squared = round(glance_res$adj.r.squared, 3),
      Model_P_Value = signif(model_p_val, 3)
    )
  } else {
    NULL
  }
}

# Sum table
model_summary_table <- expand_grid(spec = spec_vars_use, plant = topdown_vars) %>%
  pmap_dfr(~ run_model_summary(.x, .y))

if(nrow(model_summary_table) == 0) {
  stop("No significant models found matching criteria.")
}

# Save CSV backup
write.csv(model_summary_table, "Modeling/data_out/significant_lmm_model_summary_layered.csv", row.names = FALSE)

# Display
display_table_df <- model_summary_table %>%
  rename(
    "Plant Metric" = Plant_Metric,
    "Spectral Metric" = Spectral_Metric,
    "R²" = R_squared,
    "Adj R²" = Adj_R_squared,
    "Model p-value" = Model_P_Value
  )

# Table plot
table_plot <- ggtexttable(
  display_table_df, 
  rows = NULL, 
  theme = ttheme(
    padding = unit(c(6, 6), "mm"),
    colnames.style = colnames_style(fill = "gray80", face = "bold", size = 11),
    tbody.style = tbody_style(size = 10)
  )
)

print(table_plot)

# Save as JPG
num_rows <- nrow(display_table_df)
fig_height <- max(4, num_rows * 0.4)

ggsave(
  filename = "Modeling/figures/significant_lmm_summary_table_Layered.jpg", 
  plot = table_plot, 
  width = 10, 
  height = fig_height, 
  dpi = 300
)

print(paste0("Successfully saved clean summary table to Modeling/figures/significant_lmm_summary_table_Layered.jpg"))


# --------------------------------------------------------------------------------------------------------


dat_model <- dat %>%
  mutate(
    scaled_tree_cover = as.numeric(scale(pct_visible_tree_cover)),
    scaled_canopy_height = as.numeric(scale(mean_canopy_height_site)),
    scaled_stems = as.numeric(scale(visible_stems_per_ha)),
    scaled_sum_squares = as.numeric(scale(sum_squares))
  )

spec_vars_use <- c("sum_squares", "Beta_dispersion", "Beta_avg_pairwise",
                   "beta_agg_sum_squares", "beta_agg_dispersion", "Beta_agg_avg_pairwise")

topdown_vars <- c("PD_Faith_TopDown", "PD_Beta_Sorensen_TopDown", "PD_Beta_Dispersion_TopDown",
                  "FD_FDis_m1_TopDown", "FD_RaoQ_m1_TopDown", "FD_Beta_Dispersion_m1_TopDown", "FD_Beta_Turnover_m1_TopDown")

layered_vars <- c("PD_Faith_Layered", "PD_Beta_Sorensen_Layered", "PD_Beta_Dispersion_Layered",
                  "FD_FDis_m1_Layered", "FD_RaoQ_m1_Layered", "FD_Beta_Dispersion_m1_Layered", "FD_Beta_Turnover_m1_Layered")

all_plant_vars <- c(
  setNames(topdown_vars, rep("TopDown", length(topdown_vars))),
  setNames(layered_vars, rep("Layered", length(layered_vars)))
)

#  Fit models and get sum stats
run_model_summary <- function(spec, plant, extraction_type) {
  formula_str <- paste0(plant, " ~ ", spec, " * scaled_tree_cover + scaled_canopy_height")
  
  model <- tryCatch({
    lm(as.formula(formula_str), data = dat_model)
  }, error = function(e) { NULL })
  
  if (is.null(model)) return(NULL)
  
  glance_res <- broom::glance(model)
  model_p_val <- glance_res$p.value
  
  if (!is.na(model_p_val) && model_p_val < 0.05) {
    sig_stars <- case_when(
      model_p_val < 0.001 ~ "***",
      model_p_val < 0.01  ~ "**",
      model_p_val < 0.05  ~ "*",
      TRUE ~ ""
    )
    
    tibble(
      Plant_Metric = plant,
      Metric_Type = extraction_type,
      Spectral_Metric = spec,
      R_squared = glance_res$r.squared,
      Adj_R_squared = glance_res$adj.r.squared,
      Model_P_Value = model_p_val,
      Significance = sig_stars
    )
  } else {
    NULL
  }
}

# Do it for all combs
model_summary_table <- imap_dfr(all_plant_vars, function(plant_var, type_label) {
  expand_grid(spec = spec_vars_use, plant = plant_var) %>%
    pmap_dfr(~ run_model_summary(.x, .y, type_label))
})

if(nrow(model_summary_table) == 0) {
  stop("No significant models found matching criteria.")
}

write.csv(model_summary_table, "Modeling/data_out/significant_lmm_model_summary_combined.csv", row.names = FALSE)

# heatmap go time
heatmap_plot <- ggplot(model_summary_table, aes(x = Spectral_Metric, y = Plant_Metric, fill = Adj_R_squared)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Significance), color = "black", size = 4.5, vjust = 0.7) +
  facet_wrap(~ Metric_Type, scales = "free_y", ncol = 1) +
  scale_fill_gradient2(
    low = "white",
    mid = "#21908CFF",
    high = "#440154FF",
    midpoint = max(model_summary_table$Adj_R_squared, na.rm = TRUE) / 2,
    name = "Adjusted R²"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Significant LMM Model Performance Heatmap (TopDown & Layered)",
    subtitle = "Color = Adjusted R² | Stars indicate p-value (* p < 0.05, ** p < 0.01, *** p < 0.001)",
    x = "EMIT Spectral Metric",
    y = "NEON Plant Metric"
  )

print(heatmap_plot)

# Save as JPG
ggsave(
  filename = "Modeling/figures/significant_lmm_HEAT_HOT.jpg", 
  plot = heatmap_plot, 
  width = 10, 
  height = 10, 
  dpi = 300
)



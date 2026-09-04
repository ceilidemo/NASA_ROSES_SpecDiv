library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Load data
emit <- read.csv("Modeling/data_in/EMIT_qDTM_Beta_2024_full.csv")
neon <- read.csv("Modeling/data_in/NEON_TopDown_qDTM_Results_2024.csv")
site <- read.csv("Modeling/data_in/Site_info_2024.csv")

# 2. Pivot EMIT data wide (selecting the first file per site)
emit_wide <- emit %>%
  group_by(siteID) %>%
  filter(file_id == first(file_id)) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = siteID,
    names_from = variant,
    values_from = c(n_comm, M_beta, Ht_beta, qDT_beta, qDTM_beta),
    names_glue = "{.value}_{variant}"
  )

# 3. Merge datasets and clean NAs in spectral variables
model_data <- site %>%
  left_join(emit_wide, by = "siteID") %>%
  left_join(neon, by = "siteID") %>%
  drop_na(matches("^(n_comm|M_beta|Ht_beta|qDT_beta|qDTM_beta)_"))

# 4. Define predictor columns (BrightVectorNorm ONLY) and response columns (NEON Beta)
predictor_cols <- grep("BrightVectorNorm$", names(emit_wide), value = TRUE)
response_cols <- grep("^Beta_", names(neon), value = TRUE)

pred_data <- model_data %>% select(all_of(predictor_cols)) %>% select(where(is.numeric))
resp_data <- model_data %>% select(all_of(response_cols)) %>% select(where(is.numeric))

# 5. Compute cross-correlations and p-values
cross_cor <- cor(pred_data, resp_data, use = "pairwise.complete.obs")

get_p_values <- function(x, y) {
  p_mat <- matrix(NA, ncol = ncol(y), nrow = ncol(x))
  rownames(p_mat) <- colnames(x)
  colnames(p_mat) <- colnames(y)
  for (i in colnames(x)) {
    for (j in colnames(y)) {
      test <- cor.test(x[[i]], y[[j]], use = "complete.obs")
      p_mat[i, j] <- test$p.value
    }
  }
  return(p_mat)
}

p_values <- get_p_values(pred_data, resp_data)

# 6. Define explicit readable label mappings for both Predictors and Responses 
# so equivalent qDTM components (Magnitude/Dispersion, Evenness, Synthesis Diversity) align clearly.
spectral_labels <- c(
  "n_comm_BrightVectorNorm"    = "Spectral Pixel Count (n_comm)",
  "M_beta_BrightVectorNorm"    = "Spectral Dispersion / Magnitude (M)",
  "Ht_beta_BrightVectorNorm"   = "Spectral Entropy / Evenness (Ht)",
  "qDT_beta_BrightVectorNorm"  = "Spectral Effective Diversity (qDT)",
  "qDTM_beta_BrightVectorNorm" = "Spectral Synthesis Diversity (qDTM)"
)

response_labels <- c(
  "Beta_FD_qDTM"        = "Functional Synthesis Diversity (qDTM)",
  "Beta_PD_qDTM"        = "Phylogenetic Synthesis Diversity (qDTM)",
  "Beta_Dispersion_FD"  = "Functional Dispersion / Magnitude (M)",
  "Beta_Dispersion_PD"  = "Phylogenetic Dispersion / Magnitude (M)",
  "Beta_Evenness_FD"    = "Functional Evenness (Ht)",
  "Beta_Evenness_PD"    = "Phylogenetic Evenness (Ht)"
)

# 7. Reshape for plotting with significance filtering and comprehensive label re-coding
cor_df <- cross_cor %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Predictor") %>%
  pivot_longer(-Predictor, names_to = "Response", values_to = "Correlation")

pval_df <- p_values %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Predictor") %>%
  pivot_longer(-Predictor, names_to = "Response", values_to = "p_value")

plot_df <- left_join(cor_df, pval_df, by = c("Predictor", "Response")) %>%
  mutate(
    Predictor = recode(Predictor, !!!spectral_labels),
    Response  = recode(Response, !!!response_labels),
    Correlation_sig = ifelse(p_value < 0.05, Correlation, NA),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# 8. Plot Heatmap with cleanly structured component-matched labels
ggplot(plot_df, aes(x = Response, y = Predictor, fill = Correlation_sig)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Significance), color = "black", size = 3) +
  scale_fill_gradient2(low = "#2b5c8f", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), na.value = "grey95") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank()
  ) +
  labs(
    title = "EMIT Spectral Diversity vs. NEON Plant Beta Diversity Components",
    subtitle = "Non-significant pairs (p >= 0.05) greyed out; (* p<0.05, ** p<0.01, *** p<0.001)",
    x = "NEON Plant Beta Response Metric",
    y = "EMIT Spectral Predictor Metric",
    fill = "Pearson r"
  )

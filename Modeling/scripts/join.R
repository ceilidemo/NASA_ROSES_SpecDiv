# Data visualization modeling
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

emit <- read.csv("Modeling/data_in/EMIT_qDTM_Beta_2024_USE.csv")
neon <- read.csv("Modeling/data_in/NEON_Plants_qDTM_2024.csv")
gedi <- read.csv("Modeling/data_in/GEDI_2024_250m_q1_ALL_normalized_model_predictors.csv")

emit_pro<- emit %>%
  filter(variant == "BrightVectorNorm") %>%
  group_by(siteID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

#==========================================================
## EMIT and NEON
library(dplyr)
library(ggplot2)
library(broom)

# 1. Prepare EMIT predictor data
emit_predictors <- emit %>%
  filter(variant == "BrightVectorNorm") %>%
  group_by(siteID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# 2. Get all unscaled NEON beta metrics
beta_cols <- grep("beta_unscaled", names(neon), value = TRUE)

# 3. Loop through each unscaled metric and plot with cleaned labels removing _unscaled
for (n_metric in beta_cols) {
  
  # Determine the matching EMIT variable based on the metric name
  e_var <- NULL
  if (grepl("qDTM", n_metric)) {
    e_var <- "qDTM"
  } else if (grepl("Et", n_metric)) {
    e_var <- ifelse("qEt" %in% names(emit_predictors), "qEt", "Et")
  } else if (grepl("prime", n_metric)) {
    e_var <- "M_prime"
  } else if (grepl("M", n_metric) && !grepl("prime", n_metric)) {
    e_var <- "M"
  }
  
  if (!is.null(e_var) && e_var %in% names(emit_predictors)) {
    sub_df <- neon %>%
      select(siteID, all_of(n_metric)) %>%
      inner_join(emit_predictors %>% select(siteID, all_of(e_var)), by = "siteID") %>%
      rename(neon_val = all_of(n_metric), emit_val = all_of(e_var)) %>%
      filter(!is.na(neon_val), !is.na(emit_val))
    
    if (nrow(sub_df) > 5) {
      fit <- lm(neon_val ~ emit_val, data = sub_df)
      r2 <- glance(fit)$r.squared
      p_val <- tidy(fit)$p.value[2]
      
      # Clean string by removing "_unscaled"
      clean_n_metric <- sub("_unscaled", "", n_metric)
      
      p <- ggplot(sub_df, aes(x = emit_val, y = neon_val)) +
        geom_point(color = "gray40", size = 2.5) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "gray70", alpha = 0.3, se = TRUE, linewidth = 0.9) +
        theme_minimal() +
        labs(
          title = sprintf("NEON.%s ~ EMIT.%s", clean_n_metric, e_var),
          subtitle = sprintf("R² = %.2f, p = %.3f", r2, p_val),
          x = sprintf("EMIT.%s", e_var),
          y = sprintf("NEON.%s", clean_n_metric)
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          axis.text = element_text(size = 9)
        )
      
      print(p)
    }
  }
}

#=================================================================
# Model incorporating Gedi LASSO BABY

library(dplyr)
library(tidyr)
library(glmnet)
library(broom)

# 1. Prepare EMIT predictor data (site-level means)
emit_predictors <- emit %>%
  filter(variant == "BrightVectorNorm") %>%
  group_by(siteID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# 2. Prepare GEDI predictor data (site-level)
gedi_predictors <- gedi

# 3. Get all NEON beta metrics (FD and PD, scaled and unscaled)
beta_cols <- grep("beta", names(neon), value = TRUE)

# 4. Merge all datasets by siteID
combined_data <- neon %>%
  select(siteID, all_of(beta_cols)) %>%
  inner_join(emit_predictors, by = "siteID", suffix = c("", "_emit")) %>%
  inner_join(gedi_predictors, by = "siteID", suffix = c("", "_gedi"))

# Define predictor feature matrix columns (all numeric EMIT and GEDI columns)
predictor_cols <- setdiff(names(combined_data), c("siteID", beta_cols))

# 5. Loop through every beta metric and run 10-fold CV LASSO
lasso_results_list <- list()

set.seed(42)
for (target_metric in beta_cols) {
  
  # Subset complete cases for this specific target to maximize available data
  modeling_df <- combined_data %>%
    select(siteID, all_of(target_metric), all_of(predictor_cols)) %>%
    na.omit()
  
  if (nrow(modeling_df) > 10) {
    x_mat <- as.matrix(modeling_df %>% select(all_of(predictor_cols)))
    y_vec <- modeling_df[[target_metric]]
    
    if (length(unique(y_vec)) > 1) {
      # Fit cross-validated LASSO
      cv_lasso <- cv.glmnet(x_mat, y_vec, alpha = 1, nfolds = 10)
      
      # Extract non-zero coefficients at optimal lambda.min
      lasso_coefs <- coef(cv_lasso, s = "lambda.min")
      coef_df <- data.frame(
        neon_metric = target_metric,
        variable = rownames(lasso_coefs),
        coefficient = as.vector(lasso_coefs)
      ) %>%
        filter(variable != "(Intercept)", coefficient != 0)
      
      if (nrow(coef_df) > 0) {
        lasso_results_list[[target_metric]] <- coef_df
      }
    }
  }
}

# Combine all selected features across all beta metrics into a single summary table
all_lasso_summary <- bind_rows(lasso_results_list)
print(all_lasso_summary)


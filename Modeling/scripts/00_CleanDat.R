# 00_CleanDat: Bring plant and spec div together

library(tidyverse)
library(dplyr)

# -------------------------------------------------------------------------------------
# Load the EMIT results for specdiv
emit_2024 <- read.csv("EMIT_SpecDiv/data_out/results_2024/EMIT_SpecDiv_2024_new.csv") %>% 
  mutate(year = 2024)
emit_2023 <- read.csv("EMIT_SpecDiv/data_out/results_2023/EMIT_SpecDiv_2023_new.csv") %>% 
  mutate(year = 2023)

emit_data <- bind_rows(emit_2023, emit_2024) %>% 
  mutate(
    siteID = str_extract(site, "[A-Z]{4}"),
    tile_num = if_else(is.na(str_extract(site, "\\d+$")), "1", str_extract(site, "\\d+$"))
  ) %>% 
  dplyr::select(siteID, year, tile_num, sum_squares, Beta_dispersion, Beta_avg_pairwise,
                gamma_sum_squares, gamma_dispersion, beta_agg_sum_squares, beta_agg_dispersion, Beta_agg_avg_pairwise,
                plot_sum_squares, plot_beta_dispersion, plot_avg_pairwise)

# -------------------------------------------------------------------------------------
# Load the NEON plant div results
neon_layered <- read.csv("NEON_PlantDiv/data_out/NEON_Layered_FullPlantDiv_2024.csv")
neon_topdown <- read.csv("NEON_PlantDiv/data_out/NEON_TopDown_FullPlantDiv_2024.csv")

NEON_data <- full_join(neon_layered, neon_topdown, by = "siteID", suffix = c("_Layered", "_TopDown")) %>%
  mutate(year = 2024)

# and the site info
site_info <- read.csv("NEON_PlantDiv/data_out/Site_info_2024.csv")

# -------------------------------------------------------------------------------------
# join everything up
dat <- NEON_data %>%
  left_join(emit_data, by = c("siteID", "year")) %>%
  drop_na(sum_squares) %>%
  left_join(site_info, by = "siteID") %>%
  mutate(across(any_of(c("total_visible_stems", "pct_visible_tree_cover", "visible_stems_per_ha")), 
                ~replace_na(., 0))) %>%  
  mutate(site_year = paste0(siteID, "_", year))

# keep only one tile per each (recent yr)
set.seed(42) 

dat_final <- dat %>%  
  group_by(siteID) %>%  
  filter(year == max(year)) %>%  
  slice_sample(n = 1) %>%  
  ungroup()

write.csv(dat_final, "Modeling/data_work/Joined_Results.csv",   row.names = FALSE)

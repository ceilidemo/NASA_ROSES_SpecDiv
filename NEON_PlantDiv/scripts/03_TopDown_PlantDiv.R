## 03_TopDown_PlantDiv: "Top-down" view plant diversity calcs

library(dplyr)
library(tidyr)
library(sf)
library(purrr)
library(vegan)
library(sf)

# Get community and phylo data
load("data_work/community_2024.RData")
load("data_work/phylo.RData")
load("data_out/layered_results_2024.RData")

# clipping function dealing with crown dynamics etc 
clip_overlapping_crowns <- function(df, plot_bbox) {
  if (nrow(df) == 0) return(NULL)
  
  visible_geoms <- list()
  covered_area <- st_sfc(crs = st_crs(df))
  
  for (i in 1:nrow(df)) {
    current_tree <- df[i, ]
    geom <- st_make_valid(current_tree$geometry)
    geom <- st_intersection(geom, plot_bbox)
    
    if (length(covered_area) == 0) {
      visible_part <- geom
    } else {
      mask <- st_make_valid(st_union(covered_area))
      visible_part <- st_difference(geom, mask)
    }
    
    if (length(visible_part) > 0 && !st_is_empty(visible_part)) {
      current_tree$geometry <- visible_part
      current_tree$visible_area <- as.numeric(st_area(visible_part))
      
      if (current_tree$visible_area > 0.01) {
        visible_geoms[[length(visible_geoms) + 1]] <- current_tree
        covered_area <- st_make_valid(st_union(covered_area, visible_part))
      }
    }
  }
  return(if(length(visible_geoms) > 0) bind_rows(visible_geoms) else NULL)
}

# set up spatial enviro
plot_bbox <- st_polygon(list(matrix(c(-10,-10, 10,-10, 10,10, -10,10, -10,-10), ncol=2, byrow=TRUE))) %>% 
  st_sfc()

# map thee trees
canopy_polys <- structure_all$vst_mappingandtagging %>%
  distinct(individualID, .keep_all = TRUE) %>%
  inner_join(trees_df %>% dplyr::select(individualID, height, area_m2), by = "individualID") %>%
  mutate(theta = stemAzimuth * pi / 180, 
         x = stemDistance * sin(theta), 
         y = stemDistance * cos(theta),
         rad = sqrt(area_m2/pi)) %>%
  filter(!is.na(x)) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_buffer(dist = .$rad)

# run that func
canopy_visible <- canopy_polys %>% 
  arrange(plotID, desc(height)) %>% 
  group_by(plotID) %>% 
  group_split() %>% 
  map_dfr(~clip_overlapping_crowns(.x, plot_bbox))

# calc gap fraction (what parts are not covered by trees)
gap_summary <- canopy_visible %>% 
  st_drop_geometry() %>% 
  group_by(plotID) %>% 
  summarize(tree_area = sum(visible_area), .groups = "drop") %>%
  left_join(shrubs_df %>% group_by(plotID) %>% summarize(s_area = sum(area_m2), .groups = "drop"), by = "plotID") %>%
  mutate(gap_frac = pmax(0, (400 - tree_area - replace_na(s_area, 0)) / 400))

# using that logic now we need new community for top-down
# get rid of abiotic: 
abiotic_ids <- c("lichen", "moss", "fungi", "biocrustMoss", "biocrustLichen", 
                 "biocrustLightCyanobacteria", "biocrustDarkCyanobacteria", 
                 "standingDead", "standingDeadWoody", "standingDeadHerbaceous", 
                 "litter", "wood", "soil", "rock", "scat", "bare ground", "water", "otherNonVascular")

# get comm
topdown_ready <- bind_rows(
  canopy_visible %>% st_drop_geometry() %>% dplyr::select(siteID, plotID, taxonID, area_m2 = visible_area),
  shrubs_df %>% dplyr::select(siteID, plotID, taxonID, area_m2),
  ground_df %>% left_join(gap_summary, by = "plotID") %>% 
    mutate(area_m2 = (mean_pct/100) * 400 * replace_na(gap_frac, 1))
) %>%
  filter(!taxonID %in% abiotic_ids) %>%
  left_join(taxon_lookup %>% dplyr::select(taxonID, phyloName), by = "taxonID") %>%
  mutate(taxonID = coalesce(phyloName, taxonID)) %>%
  inner_join(plot_metadata %>% filter(plotType == "distributed") %>% dplyr::select(plotID), by = "plotID") %>%
  group_by(siteID, plotID, taxonID) %>%
  summarize(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop")

# calc div
topdown_td <- topdown_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_t_div(.x))

topdown_pd <- topdown_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_p_div(.x, tre$scenario.3))

# join and save
final_topdown_results <- topdown_td %>%
  left_join(topdown_pd, by = "siteID") %>%
  mutate(across(where(is.numeric), ~round(., 4)))

write.csv(final_topdown_results, "data_out/NEON_TopDown_PlantDiv_Results_2024.csv", row.names = FALSE)

save(topdown_ready, final_topdown_results, file = "data_out/topdown_results_2024.RData")

###################################################
###################################################
# MAPPING CHECK #

site_to_check <- "MOAB" 

site_map_data <- canopy_visible %>% 
  filter(grepl(site_to_check, plotID))

ggplot() +
  geom_sf(data = plot_bbox, fill = "gray95", color = "black", linetype = "dashed") +
  geom_sf(data = site_map_data, aes(fill = taxonID), alpha = 0.8, color = "white", size = 0.1) +
  facet_wrap(~plotID, ncol = 8) + 
  scale_fill_viridis_d(option = "plasma") + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = paste("'Top-Down' View of plots in:", site_to_check),
    subtitle = "Visible mapped canopy in each plot w/species ID",
    fill = "Taxon ID",
    x = NULL, y = NULL
  )

######################################
#### export canopy metrics / site info

plot_height_metrics <- structure_all$vst_apparentindividual %>%
  filter(!is.na(height)) %>%
  group_by(plotID) %>%
  summarise(
    mean_canopy_height = mean(height, na.rm = TRUE),
    max_canopy_height  = max(height, na.rm = TRUE),
    sd_canopy_height   = sd(height, na.rm = TRUE),
    mean_base_crown_ht = mean(baseCrownHeight, na.rm = TRUE),
    .groups = "drop"
  )

# Get plot-level elevation from perplotperyear
plot_elevation <- structure_all$vst_perplotperyear %>%
  select(plotID, elevation) %>%
  distinct(plotID, .keep_all = TRUE)

# num distributed plots per site
site_plots_count <- plot_metadata %>%
  filter(plotType == "distributed") %>%
  mutate(siteID = substr(plotID, 1, 4)) %>%  
  group_by(siteID) %>%
  summarize(n_plots = n(), .groups = "drop")

# Merge spatial with canopy
canopy_shown_plots <- canopy_visible %>%
  st_drop_geometry() %>%
  inner_join(plot_metadata %>% filter(plotType == "distributed") %>% select(plotID), by = "plotID") %>%
  left_join(plot_height_metrics, by = "plotID") %>%
  left_join(plot_elevation, by = "plotID")

# aggregate to site-level
site_canopy_metrics <- canopy_shown_plots %>%
  group_by(siteID) %>%
  summarize(
    total_visible_stems         = n_distinct(individualID),
    total_site_canopy_area      = sum(visible_area, na.rm = TRUE),
    mean_canopy_height_site     = mean(mean_canopy_height, na.rm = TRUE),
    max_canopy_height_site      = max(max_canopy_height, na.rm = TRUE),
    canopy_height_sd_site       = mean(sd_canopy_height, na.rm = TRUE),
    mean_base_crown_height_site = mean(mean_base_crown_ht, na.rm = TRUE),
    mean_elevation              = mean(elevation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(site_plots_count, by = "siteID") %>%
  mutate(
    avg_visible_canopy_per_plot = total_site_canopy_area / n_plots,
    pct_visible_tree_cover      = (avg_visible_canopy_per_plot / 400) * 100,
    visible_stems_per_ha        = total_visible_stems / (n_plots * 0.04)
  )

# clean and write
final_site_info <- site_plots_count %>%
  left_join(site_canopy_metrics, by = c("siteID", "n_plots")) %>%
  mutate(across(
    c(total_visible_stems, avg_visible_canopy_per_plot, pct_visible_tree_cover, 
      visible_stems_per_ha, mean_canopy_height_site, max_canopy_height_site, 
      canopy_height_sd_site, mean_base_crown_height_site, mean_elevation),
    ~ replace_na(., 0)
  )) %>%
  mutate(across(where(is.numeric), ~ round(., 5)))

write.csv(final_site_info, "data_out/Site_info_2024.csv", row.names = FALSE)


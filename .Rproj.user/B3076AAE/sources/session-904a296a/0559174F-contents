## 04_FuncDiv : Functional Diversity of Plants across NEON sites
# Ceili DeMarais

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(betapart)
library(ape)
library(FD)
library(vegan)

load("data_in/foliar_data_2024.RData")
load("data_work/community_2024.RData")
load("data_out/layered_results_2024.RData")
load("data_out/topdown_results_2024.RData")

## Build NEON Trait Table
# NEON taxonID -> phyloName used in community data
phylo_to_neon <- taxon_lookup %>%
  select(taxonID, phyloName) %>%
  distinct() %>%
  filter(!is.na(phyloName))

# Site-specific traits with cross-site fallback, keyed by community_taxonID
trait_cols <- c("LMA_gm2", "N_pct", "C_pct", "chlorophyll_mgg",
                "lignin_pct", "cellulose_pct", "EWT_g_cm2", 
                "d13C", "d15N", "carotenoids_mgg")

# attach community_taxonID to both site and cross-site tables
neon_traits_site_keyed <- neon_traits_site %>%
  left_join(phylo_to_neon, by = "taxonID") %>%
  mutate(community_taxonID = coalesce(phyloName, taxonID))

neon_traits_cross_keyed <- neon_traits_fullCross %>%
  left_join(phylo_to_neon, by = "taxonID") %>%
  mutate(community_taxonID = coalesce(phyloName, taxonID))


## Build foliar trait matrices

# trait pools for models
traits_m1 <- neon_traits_cross_keyed %>%
  filter(!is.na(LMA_gm2), !is.na(N_pct), !is.na(C_pct)) %>%
  distinct(community_taxonID, .keep_all = TRUE) %>%
  select(community_taxonID, LMA = LMA_gm2, N = N_pct, C = C_pct) %>%
  column_to_rownames("community_taxonID")

traits_m2 <- neon_traits_cross_keyed %>%
  drop_na(LMA_gm2, N_pct, C_pct, chlorophyll_mgg,
          lignin_pct, cellulose_pct, EWT_g_cm2) %>%
  distinct(community_taxonID, .keep_all = TRUE) %>%
  select(community_taxonID,
         LMA       = LMA_gm2,
         N         = N_pct,
         C         = C_pct,
         Chl       = chlorophyll_mgg,
         Lignin    = lignin_pct,
         Cellulose = cellulose_pct,
         Water     = EWT_g_cm2) %>%
  column_to_rownames("community_taxonID")

# Filter the data to sites with coverage
layered_keyed <- layered_ready %>%
  left_join(phylo_to_neon, by = "taxonID") %>%
  mutate(community_taxonID = coalesce(phyloName, taxonID))

sites_to_keep <- layered_keyed %>%
  group_by(siteID) %>%
  summarize(
    pct_area_covered = 100 * sum(area_m2[community_taxonID %in% rownames(traits_m1)], na.rm = TRUE) / sum(area_m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(pct_area_covered >= 80.0) %>%
  pull(siteID)

cat("Keeping", length(sites_to_keep), "sites that have 80% area coverage for foliar traits.\n")

# -------------------------------------------------------------------------------------
# Model 1: LMA + N + C 
comm_m1_raw <- layered_keyed %>%
  filter(siteID %in% sites_to_keep, community_taxonID %in% rownames(traits_m1)) %>%
  group_by(siteID, community_taxonID) %>%
  summarize(total_area = sum(area_m2), .groups = "drop")

comm_m1 <- comm_m1_raw %>%
  pivot_wider(names_from = community_taxonID, values_from = total_area, values_fill = 0) %>%
  column_to_rownames("siteID")

shared_m1 <- intersect(rownames(traits_m1), colnames(comm_m1))
traits_m1_final <- traits_m1[shared_m1, ]
comm_m1         <- comm_m1[, shared_m1]

# -------------------------------------------------------------------------------------
# Model 2: all 7 traits

comm_m2_raw <- layered_keyed %>%
  filter(siteID %in% sites_to_keep, community_taxonID %in% rownames(traits_m2)) %>%
  group_by(siteID, community_taxonID) %>%
  summarize(total_area = sum(area_m2), .groups = "drop")

comm_m2 <- comm_m2_raw %>%
  pivot_wider(names_from = community_taxonID, values_from = total_area, values_fill = 0) %>%
  column_to_rownames("siteID")

shared_m2 <- intersect(rownames(traits_m2), colnames(comm_m2))
traits_m2_final <- traits_m2[shared_m2, ]
comm_m2         <- comm_m2[, shared_m2]

# -------------------------------------------------------------------------------------
# Functional diversity calculations

calc_f_div <- function(df, site_id, model_label, suffix, trait_cols) {
  
  empty_res <- data.frame(
    FD_FDis = NA_real_,
    FD_RaoQ = NA_real_,
    FD_Beta_Dispersion = NA_real_,
    FD_Beta_Turnover = NA_real_
  ) %>% rename_with(~paste0(., "_", model_label, "_", suffix))
  
  df_keyed <- df %>%
    left_join(phylo_to_neon, by = "taxonID") %>%
    mutate(community_taxonID = coalesce(phyloName, taxonID))
  
  # Build site-specific trait matrix: site mean where available, cross-site fallback
  site_traits_raw <- neon_traits_cross_keyed %>%
    select(community_taxonID, all_of(trait_cols)) %>%
    left_join(
      neon_traits_site_keyed %>%
        filter(siteID == site_id) %>%
        select(community_taxonID, all_of(trait_cols)),
      by = "community_taxonID",
      suffix = c("_global", "_site")
    ) %>%
    mutate(
      LMA_gm2         = coalesce(LMA_gm2_site,         LMA_gm2_global),
      N_pct           = coalesce(N_pct_site,           N_pct_global),
      C_pct           = coalesce(C_pct_site,           C_pct_global),
      chlorophyll_mgg = coalesce(chlorophyll_mgg_site, chlorophyll_mgg_global),
      lignin_pct      = coalesce(lignin_pct_site,      lignin_pct_global),
      cellulose_pct   = coalesce(cellulose_pct_site,   cellulose_pct_global),
      EWT_g_cm2       = coalesce(EWT_g_cm2_site,       EWT_g_cm2_global)
    ) %>%
    select(community_taxonID, all_of(trait_cols)) %>%
    distinct(community_taxonID, .keep_all = TRUE)
  
  # Apply model trait filters
  if (model_label == "m1") {
    traits_use <- site_traits_raw %>%
      filter(!is.na(LMA_gm2), !is.na(N_pct), !is.na(C_pct)) %>%
      select(community_taxonID, LMA = LMA_gm2, N = N_pct, C = C_pct)
  } else {
    traits_use <- site_traits_raw %>%
      drop_na(all_of(trait_cols)) %>%
      select(community_taxonID,
             LMA = LMA_gm2, N = N_pct, C = C_pct,
             Chl = chlorophyll_mgg, Lignin = lignin_pct,
             Cellulose = cellulose_pct, Water = EWT_g_cm2)
  }
  
  if(nrow(traits_use) < 2) return(empty_res)
  traits_matrix <- traits_use %>% column_to_rownames("community_taxonID")
  
  comm_raw <- df_keyed %>%
    group_by(plotID, community_taxonID) %>%
    summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
    filter(community_taxonID %in% rownames(traits_matrix)) %>%
    pivot_wider(names_from = community_taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  #Need at least 2 plots and 2 species per site
  comm_filtered <- comm_raw[rowSums(comm_raw > 0) >= 2, , drop = FALSE]
  if(nrow(comm_filtered) < 1) return(empty_res)
  
  species_present <- colnames(comm_filtered)[colSums(comm_filtered) > 0]
  if(length(species_present) < 2) return(empty_res)
  
  comm <- comm_filtered[, species_present, drop = FALSE]
  site_traits_mat <- traits_matrix[species_present, , drop = FALSE]
  site_traits_scaled <- scale(site_traits_mat)
  
  # Run dbFD
  fd_out <- FD::dbFD(site_traits_scaled,
                     as.matrix(comm),
                     corr = "cailliez",
                     calc.FRic = TRUE,
                     stand.x = FALSE,
                     messages = FALSE)
  
  fdis_val <- mean(fd_out$FDis, na.rm = TRUE)
  raoq_val <- mean(fd_out$RaoQ, na.rm = TRUE)
  
  # Functional Beta Diversity (Dispersion)
  beta_disp <- NA
  if(nrow(comm) >= 2) {
    comm_norm   <- decostand(comm, method = "total")
    plot_traits <- as.matrix(comm_norm) %*% site_traits_scaled
    centroid    <- colMeans(plot_traits)
    beta_disp   <- mean(sqrt(rowSums((sweep(plot_traits, 2, centroid))^2)), na.rm = TRUE)
  }
  
  # Functional Beta Diversity (Turnover/Distance)
  beta_turn <- NA
  if(nrow(comm) >= 2) {
    beta_turn <- tryCatch({
      comm_norm   <- decostand(comm, method = "total")
      plot_traits <- as.matrix(comm_norm) %*% site_traits_scaled
      mean(dist(plot_traits, method = "euclidean"), na.rm = TRUE)
    }, error = function(e) NA)
  }
  
  return(data.frame(
    FD_FDis = fdis_val,
    FD_RaoQ = raoq_val,
    FD_Beta_Dispersion = beta_disp,
    FD_Beta_Turnover = beta_turn
  ) %>% rename_with(~paste0(., "_", model_label, "_", suffix)))
}

# -------------------------------------------------------------------------------------
# Run functional div with site filters

trait_cols <- c("LMA_gm2", "N_pct", "C_pct", "chlorophyll_mgg",
                "lignin_pct", "cellulose_pct", "EWT_g_cm2")

# only the sites that met your 80% coverage criteria
layered_fd_m1 <- layered_ready %>%
  filter(siteID %in% sites_to_keep) %>%
  group_by(siteID) %>%
  group_modify(~ calc_f_div(.x, site_id = .y$siteID, "m1", "Layered", trait_cols))

layered_fd_m2 <- layered_ready %>%
  filter(siteID %in% sites_to_keep) %>%
  group_by(siteID) %>%
  group_modify(~ calc_f_div(.x, site_id = .y$siteID, "m2", "Layered", trait_cols))

topdown_fd_m1 <- topdown_ready %>%
  filter(siteID %in% sites_to_keep) %>%
  group_by(siteID) %>%
  group_modify(~ calc_f_div(.x, site_id = .y$siteID, "m1", "TopDown", trait_cols))

topdown_fd_m2 <- topdown_ready %>%
  filter(siteID %in% sites_to_keep) %>%
  group_by(siteID) %>%
  group_modify(~ calc_f_div(.x, site_id = .y$siteID, "m2", "TopDown", trait_cols))

# -------------------------------------------------------------------------------------
# Bring it together and save

layered_final_all <- layered_master_table %>%
  filter(siteID %in% sites_to_keep) %>%
  left_join(layered_fd_m1, by = "siteID") %>%
  left_join(layered_fd_m2, by = "siteID")

topdown_final_all <- final_topdown_results %>%
  filter(siteID %in% sites_to_keep) %>%
  left_join(topdown_fd_m1, by = "siteID") %>%
  left_join(topdown_fd_m2, by = "siteID")

write.csv(layered_final_all, "data_out/NEON_Layered_FullPlantDiv_2024.csv", row.names = FALSE)
write.csv(topdown_final_all, "data_out/NEON_TopDown_FullPlantDiv_2024.csv", row.names = FALSE)

save(layered_final_all, topdown_final_all,
     traits_m1, comm_m1, traits_m2, comm_m2,
     file = "data_out/plantDiv_finalAll.RData")

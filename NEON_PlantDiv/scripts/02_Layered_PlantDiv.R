## 02_Layered_PlantDiv: "Layered" diversity calculation taking footprints from various levels of canopy

library(vegan)
library(tidyverse)
library(V.PhyloMaker2)
library(picante)

# load community data if not
load("data_work/community_2024.RData")

# Diversity functions:

## Taxonomic Function
calc_t_div <- function(df) {
  comm <- df %>%
    group_by(plotID, taxonID) %>%
    summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  # filter empty plots
  comm <- comm[rowSums(comm) > 0, , drop = FALSE]
  
  site_totals <- colSums(comm)
  richness <- specnumber(site_totals)
  shannon_eff <- exp(diversity(site_totals, index = "shannon"))
  simpson_eff <- 1/diversity(site_totals, index = "simpson")
  
  beta_bray <- NA
  beta_dispersion <- NA
  
  if(nrow(comm) > 1) { # only if comm over 1 
    dist_mat <- vegdist(comm, method = "bray")
    beta_bray <- mean(dist_mat, na.rm = TRUE)
    
    disp_mod <- betadisper(dist_mat, group = rep("site", nrow(comm)), type = "centroid")
    beta_dispersion <- mean(disp_mod$distances)
  }
  
  return(data.frame(TD_Richness = richness, TD_Shannon_Eff = shannon_eff, 
                    TD_Simpson_Eff = simpson_eff, TD_Beta_Bray = beta_bray,
                    TD_Beta_Dispersion = beta_dispersion))
}

## Phylogenetic Function
calc_p_div <- function(df, tre_obj) {
  require(ape)
  require(picante)
  
  #matrix setup
  comm <- df %>%
    group_by(plotID, taxonID) %>%
    summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  # tax match
  present_in_tree <- intersect(colnames(comm), tre_obj$tip.label)
  
  if(length(present_in_tree) < 2) {
    return(data.frame(PD_Faith = NA, PD_Beta_Sorensen = NA, PD_Dispersion = NA))
  }
  
  # keep only plots w/ 1 species in tree
  comm_matched <- comm[, present_in_tree, drop = FALSE]
  comm_matched <- comm_matched[rowSums(comm_matched) > 0, , drop = FALSE]
  
  if(nrow(comm_matched) == 0) return(data.frame(PD_Faith = NA, PD_Beta_Sorensen = NA, PD_Dispersion = NA))
  
  comm_pa <- ifelse(comm_matched > 0, 1, 0)
  
  # grooming
  site_tree <- keep.tip(tre_obj, present_in_tree)
  site_tree <- multi2di(site_tree)
  
  if(!is.rooted(site_tree)) {
    site_tree <- root(site_tree, outgroup = site_tree$tip.label[1], resolve.root = TRUE)
  }
  
  # scale branches to rel time (w/ try catch for errors...)
  site_tree <- tryCatch({
    chronoMPL(site_tree)
  }, error = function(e) {
    return(site_tree) })
  
  # Diversity Calcs
  pd_val <- picante::pd(comm_pa, site_tree, include.root = FALSE)
  
  p_beta_sor <- NA
  p_dispersion <- NA
  
  if(nrow(comm_pa) > 1) { # take if comm is bigger than 1
    sor_sim <- picante::phylosor(as.matrix(comm_pa), site_tree)
    p_beta_sor <- mean(as.matrix(sor_sim), na.rm = TRUE)
    
    phylo_dist <- 1 - sor_sim
    disp_mod_p <- betadisper(as.dist(phylo_dist), group = rep("site", nrow(comm_pa)), 
                             type = "centroid", add = TRUE)
    p_dispersion <- mean(disp_mod_p$distances)
  }
  
  return(data.frame(
    PD_Faith = mean(pd_val$PD, na.rm = TRUE),
    PD_Beta_Sorensen = p_beta_sor,
    PD_Beta_Dispersion = p_dispersion
  ))
}

# Generate Phylogeny
sp_list <- taxon_lookup %>%
  filter(taxonRank == "species") %>%
  rename(species = phyloName) %>%
  mutate(genus = gsub("_.*", "", species)) %>%
  distinct(species, .keep_all = TRUE)

data("nodes.info.1.TPL", package = "V.PhyloMaker2")
sp_list_corrected <- sp_list %>%
  mutate(family_new = nodes.info.1.TPL$family[match(genus, nodes.info.1.TPL$genus)]) %>%
  mutate(family = coalesce(family_new, family)) %>%
  select(species, genus, family)

tre <- phylo.maker(sp.list = sp_list_corrected, tree = GBOTB.extended.TPL, 
                   nodes = nodes.info.1.TPL, scenarios = "S3")
dist_phylo <- as.matrix(cophenetic(tre$scenario.3))

# Get the names to match the phylogeny for calcs
layered_ready <- bind_rows(
  ground_df %>% mutate(area_m2 = (mean_pct / 100) * 400),
  shrubs_df,
  trees_df %>% select(siteID, plotID, taxonID, area_m2, is_abiotic)
) %>%
  filter(!is_abiotic) %>%
  left_join(taxon_lookup %>% select(taxonID, phyloName), by = "taxonID") %>%
  mutate(taxonID = coalesce(phyloName, taxonID)) %>%
  inner_join(plot_metadata %>% filter(plotType == "distributed") %>% select(plotID), by = "plotID") %>%
  group_by(siteID, plotID, taxonID) %>%
  summarize(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(area_m2 = pmin(area_m2, 400))

# ok now actually calculate div
layered_td <- layered_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_t_div(.x))

layered_pd <- layered_ready %>%
  group_by(siteID) %>%
  group_modify(~ calc_p_div(.x, tre$scenario.3))

# join and save it
layered_master_table <- layered_td %>%
  left_join(layered_pd, by = "siteID") %>%
  mutate(across(where(is.numeric), ~round(., 4)))

write.csv(layered_master_table, "data_out/NEON_Layered_PlantDiv_Results_2024.csv", row.names = FALSE)

# save phylogeny to use for the top down view :) 
save(dist_phylo, taxon_lookup, file = "data_work/phylo.RData")

# save community objects needed by 04_FuncDiv
save(layered_ready, tre, layered_master_table, file = "data_out/layered_results_2024.RData")

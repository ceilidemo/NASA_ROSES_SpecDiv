## 01_ProcessFiles: Process the downloaded data to get working files for diversity calculations

library(dplyr)
library(tidyr)
        
# load in in the veg data 
load("data_in/veg_data_2024.RData")

# retrieve taxa from % cover and struct data
taxon_lookup <- bind_rows(
  cover_all$div_1m2Data %>% select(taxonID, scientificName, family, taxonRank),
  structure_all$vst_mappingandtagging %>% select(taxonID, scientificName, family, taxonRank)
) %>% 
  distinct(taxonID, .keep_all = TRUE) %>%
  filter(!is.na(taxonID)) %>%
  mutate(phyloName = case_when(
    taxonRank == "species" ~ gsub(" ", "_", scientificName),
    TRUE ~ taxonID
  ))

# process the ground cover data (getting most recent visit)
ground_df <- cover_all$div_1m2Data %>%
  group_by(plotID, subplotID) %>%
  slice_max(endDate, n = 1, with_ties = FALSE) %>% # Latest per subplot
  ungroup() %>%
  mutate(final_id = coalesce(taxonID, otherVariables)) %>%
  filter(!is.na(final_id)) %>%
  group_by(siteID, plotID, final_id) %>%
  summarize(mean_pct = mean(percentCover, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    is_abiotic = grepl("litter|rock|soil|water|wood|other|scat|fungi|standing|biocrust|lichen|moss", 
                       tolower(final_id)),
    source = "ground"
  ) %>%
  rename(taxonID = final_id)

# do the same for shrubs (getting most recent visit)
shrubs_df <- structure_all$vst_shrubgroup %>%
  group_by(plotID, groupID) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(siteID, plotID, taxonID) %>%
  summarize(area_m2 = sum(canopyArea, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "shrub", is_abiotic = FALSE)

# and now the mapped trees (getting most recent visit)
taxon_map <- structure_all$vst_mappingandtagging %>%
  select(individualID, taxonID) %>%
  distinct(individualID, .keep_all = TRUE)

trees_df <- structure_all$vst_apparentindividual %>%
  group_by(individualID) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(taxon_map, by = "individualID") %>%
  filter(!is.na(maxCrownDiameter)) %>%
  mutate(area_m2 = pi * (maxCrownDiameter / 2)^2,
         source = "canopy", 
         is_abiotic = FALSE)

# also grab plot metadata so you can filter by distributed baseplots later on
plot_metadata <- structure_all$vst_perplotperyear %>% 
  select(plotID, plotType) %>% distinct()


# save this compilation of most recent data
save(ground_df, shrubs_df, trees_df, taxon_lookup, plot_metadata, structure_all,
     file = "data_work/community_2024.RData")

##############################
## Process the foliar traits##

load("data_in/foliar_data.RData")

# Field metadata
field_meta <- foliar_all$cfc_fieldData %>%
  select(sampleID, siteID, plotID, taxonID, scientificName,
         plantStatus, samplingImpractical) %>%
  filter(
    is.na(samplingImpractical) | samplingImpractical %in% c("OK", "")
  ) %>%
  distinct(sampleID, .keep_all = TRUE)

# LMA & EWT
lma_ewt <- foliar_all$cfc_LMA %>%
  filter(!is.na(freshMass), !is.na(dryMass), !is.na(leafArea),
         freshMass > dryMass, leafArea > 0) %>%
  mutate(
    LMA_gm2   = leafMassPerArea,
    EWT_g_cm2 = (freshMass - dryMass) / leafArea
  ) %>%
  group_by(sampleID) %>%
  summarize(
    LMA_gm2   = mean(LMA_gm2, na.rm = TRUE),
    EWT_g_cm2 = mean(EWT_g_cm2, na.rm = TRUE),
    .groups = "drop"
  )

# C and N (averaged per sampleID if dups)
cn <- foliar_all$cfc_carbonNitrogen %>%
  select(sampleID, N_pct = nitrogenPercent, C_pct = carbonPercent, d13C, d15N) %>%
  group_by(sampleID) %>%
  summarize(
    N_pct = mean(N_pct, na.rm = TRUE),
    C_pct = mean(C_pct, na.rm = TRUE),
    d13C  = mean(d13C, na.rm = TRUE),
    d15N  = mean(d15N, na.rm = TRUE),
    .groups = "drop"
  )

# Chlorophyll 
chl <- foliar_all$cfc_chlorophyll %>%
  select(sampleID, extractChlAConc, extractChlBConc, extractCarotConc, freshMass, solventVolume, dilutionFactor) %>%
  mutate(
    chl_a_ugg       = (extractChlAConc * solventVolume * dilutionFactor) / freshMass,
    chl_b_ugg       = (extractChlBConc * solventVolume * dilutionFactor) / freshMass,
    carot_ugg       = (extractCarotConc * solventVolume * dilutionFactor) / freshMass,
    chlorophyll_mgg = (chl_a_ugg + chl_b_ugg) / 1000,
    carotenoids_mgg = carot_ugg / 1000
  ) %>%
  group_by(sampleID) %>%
  summarize(
    chlorophyll_mgg = mean(chlorophyll_mgg, na.rm = TRUE),
    carotenoids_mgg = mean(carotenoids_mgg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(chlorophyll_mgg), chlorophyll_mgg > 0)

# Lignin 
lig <- foliar_all$cfc_lignin %>%
  select(sampleID, ligninPercent, cellulosePercent) %>%
  group_by(sampleID) %>%
  summarize(
    lignin_pct    = mean(ligninPercent, na.rm = TRUE),
    cellulose_pct = mean(cellulosePercent, na.rm = TRUE),
    .groups = "drop"
  )

# Join 
neon_traits_raw <- field_meta %>%
  left_join(lma_ewt, by = "sampleID") %>%
  left_join(cn,      by = "sampleID") %>%
  left_join(chl,     by = "sampleID") %>%
  left_join(lig,     by = "sampleID")

# Get site-level means
neon_traits_site <- neon_traits_raw %>%
  group_by(siteID, taxonID, scientificName) %>%
  summarize(across(c(LMA_gm2, N_pct, C_pct, d13C, d15N,
                     chlorophyll_mgg, carotenoids_mgg,
                     lignin_pct, cellulose_pct, EWT_g_cm2),
                   ~mean(.x, na.rm = TRUE)),
            n_samples = n(), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA_real_, .)))

# Get trait means across all sites together to fill gaps
neon_traits_fullCross <- neon_traits_raw %>%
  group_by(taxonID, scientificName) %>%
  summarize(across(c(LMA_gm2, N_pct, C_pct, d13C, d15N,
                     chlorophyll_mgg, carotenoids_mgg,
                     lignin_pct, cellulose_pct, EWT_g_cm2),
                   ~mean(.x, na.rm = TRUE)),
            n_sites = n_distinct(siteID), n_samples = n(), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA_real_, .)))

# Check coverage
cat("Species with site-level traits:", n_distinct(neon_traits_site$taxonID), "\n")
cat("Sites with foliar data:", n_distinct(neon_traits_site$siteID), "\n")

print(neon_traits_fullCross %>%
        summarize(across(c(LMA_gm2, N_pct, C_pct, chlorophyll_mgg, lignin_pct, d13C, d15N,
                           carotenoids_mgg, cellulose_pct, EWT_g_cm2),
                         ~sum(!is.na(.)))))

# Save
save(foliar_all, neon_traits_site, neon_traits_fullCross, neon_traits_raw,
     file = "data_in/foliar_data_2024.RData")
write.csv(neon_traits_site,   "data_in/neon_traits_site_2024.csv",   row.names = FALSE)
write.csv(neon_traits_fullCross, "data_in/neon_traits_fullCross_2024.csv", row.names = FALSE)



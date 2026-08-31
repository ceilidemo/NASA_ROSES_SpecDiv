## qDTM Phylogenetic and Functional Diversity Workflow across NEON Sites
# Using Shan Kothari's qDTM script with top-down sensor-visible canopy abundances
# Sampling plots serve as alpha communities; site-level metrics aggregate alpha and beta diversity.
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ape)
library(vegan)
library(V.PhyloMaker2)

# Load processed core datasets
load("NEON_PlantDiv/data_in/foliar_data_2024.RData")
load("NEON_PlantDiv/data_work/community_2024.RData")

# -------------------------------------------------------------------------------------
# Shan Kothari's qDTM Functions
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# 1. Shan Kothari's qDTM Functions (Fixed Class Checks)
# -------------------------------------------------------------------------------------
FTD <- function(tdmat, weights = NULL, q = 1) {
  if (length(tdmat) == 1 && tdmat == 0) {
    tdmat <- as.matrix(tdmat)
  }
  is_mat <- is.matrix(tdmat)
  is_dst <- inherits(tdmat, "dist")
  
  if (!is_mat && !is_dst) {
    stop("distances must be class dist or class matrix")
  } else if (is_mat && !isSymmetric(unname(tdmat))) {
    warning("trait matrix not symmetric")
  } else if (is_dst) {
    tdmat <- as.matrix(tdmat)
  }
  if (!isTRUE(all.equal(sum(diag(tdmat)), 0))) {
    warning("non-zero diagonal; species appear to have non-zero trait distance from themselves")
  }
  if (max(tdmat) > 1 || min(tdmat) < 0) {
    tdmat <- (tdmat - min(tdmat)) / (max(tdmat) - min(tdmat))
    warning("trait distances must be between 0 and 1; rescaling")
  }
  if (is.null(weights)) {
    nsp <- nrow(tdmat)
    weights <- rep(1 / nsp, nsp)
  } else {
    nsp <- sum(weights > 0)
  }
  if (!isTRUE(all.equal(sum(weights), 1))) {
    weights <- weights / sum(weights)
    warning("input proportional abundances do not sum to 1; summation to 1 forced")
  }
  tdmat.abund <- diag(weights) %*% tdmat %*% diag(weights)
  M <- sum(tdmat.abund)
  M.prime <- ifelse(nsp == 1, 0, M * nsp / (nsp - 1))
  fij <- tdmat.abund / M
  
  if (isTRUE(all.equal(M, 0))) {
    qHt <- 0
  } else if (q == 1) {
    fijlog <- ifelse(fij == 0, 0, fij * log(fij))
    qHt <- exp(-1 * sum(fijlog))
  } else if (q == 0) {
    qHt <- sum(fij > 0)
  } else {
    qHt <- sum(fij^q)^(1 / (1 - q))
  }
  
  qDT <- (1 + sqrt(1 + 4 * qHt)) / 2
  qDTM <- 1 + qDT * M
  qEt <- qDT / nsp
  list(nsp = nsp, q = q, M = M, M.prime = M.prime, qHt = qHt, qEt = qEt, qDT = qDT, qDTM = qDTM)
}

FTD.comm <- function(tdmat, spmat, q = 1, abund = FALSE, match.names = FALSE) {
  is_mat <- is.matrix(tdmat)
  is_dst <- inherits(tdmat, "dist")
  
  if (!is_mat && !is_dst) {
    stop("distances must be class dist or class matrix")
  } else if (is_mat && !isSymmetric(unname(tdmat))) {
    warning("trait matrix not symmetric")
  } else if (is_dst) {
    tdmat <- as.matrix(tdmat)
  }
  if (!isTRUE(all.equal(sum(diag(tdmat)), 0))) {
    warning("non-zero diagonal; species appear to have non-zero trait distance from themselves")
  }
  if (max(tdmat) > 1 || min(tdmat) < 0) {
    tdmat <- (tdmat - min(tdmat)) / (max(tdmat) - min(tdmat))
    warning("trait distances must be between 0 and 1; rescaling")
  }
  if (abund == FALSE) {
    spmat[spmat > 0] <- 1
    spmat <- spmat / rowSums(spmat)
  }
  n.comm <- nrow(spmat)
  if (match.names == TRUE) {
    sp.arr <- match(rownames(as.matrix(tdmat)), colnames(spmat))
    spmat <- spmat[, sp.arr]
  }
  out <- apply(spmat, 1, function(x) unlist(FTD(tdmat = tdmat, weights = x, q = q)))
  df.out <- data.frame(t(out))
  rownames(df.out) <- rownames(spmat)
  if (sum(df.out$nsp == 0) > 0) {
    warning("at least one community has no species")
  }
  nsp <- sum(colSums(spmat > 0))
  u.nsp <- mean(df.out$nsp)
  u.M <- sum(df.out$nsp * df.out$M) / sum(df.out$nsp)
  if (q == 1) {
    u.qDT <- prod(df.out$qDT)^(1 / n.comm)
  } else {
    u.qDT <- (sum(df.out$qDT^(1 - q)) / n.comm)^(1 / (1 - q))
  }
  u.M.prime <- u.M * u.nsp / (u.nsp - 1)
  u.qDTM <- 1 + u.qDT * u.M
  u.qEt <- u.qDT / u.nsp
  list(com.FTD = df.out, nsp = nsp, u.nsp = u.nsp, u.M = u.M, u.M.prime = u.M.prime, u.qEt = u.qEt, u.qDT = u.qDT, u.qDTM = u.qDTM)
}

comm.disp <- function(tdmat, com1, com2) {
  is_mat <- is.matrix(tdmat)
  is_dst <- inherits(tdmat, "dist")
  
  if (!is_mat && !is_dst) {
    stop("distances must be class dist or class matrix")
  } else if (is_mat && !isSymmetric(unname(tdmat))) {
    warning("trait matrix not symmetric")
  } else if (is_dst) {
    tdmat <- as.matrix(tdmat)
  }
  if (isTRUE(all.equal(sum(com1), 1)) == FALSE || isTRUE(all.equal(sum(com2), 1)) == FALSE) {
    com1 <- com1 / sum(com1)
    com2 <- com2 / sum(com2)
    warning("input proportional abundances do not sum to 1; summation to 1 forced")
  }
  mAB <- sum(diag(com1) %*% tdmat %*% diag(com2))
  mAA <- sum(diag(com1) %*% tdmat %*% diag(com1))
  mBB <- sum(diag(com2) %*% tdmat %*% diag(com2))
  dmAB <- mAB - (mAA + mBB) / 2
  return(dmAB)
}

comm.disp.mat <- function(tdmat, spmat, abund = FALSE, sp.weighted = FALSE) {
  n.comm <- nrow(spmat)
  if (abund == FALSE) {
    spmat[spmat > 0] <- 1
    spmat <- spmat / rowSums(spmat)
  }
  if (FALSE %in% sapply(rowSums(spmat), function(x) isTRUE(all.equal(x, 1)))) {
    spmat <- spmat / rowSums(spmat)
    warning("proportional abundances don't always sum to 1; summation to 1 forced")
  }
  disp.mat <- outer(1:n.comm, 1:n.comm, FUN = Vectorize(function(i, j) comm.disp(tdmat, com1 = spmat[i, ], com2 = spmat[j, ])))
  if (sp.weighted == TRUE) {
    nsp.comm <- rowSums(spmat > 0)
    disp.mat.weight <- diag(nsp.comm) %*% disp.mat %*% diag(nsp.comm)
    return(disp.mat.weight)
  } else {
    return(disp.mat)
  }
}

FTD.beta <- function(tdmat, spmat, abund = FALSE, q = 1) {
  nsp <- sum(colSums(spmat) > 0)
  St <- sum(spmat > 0)
  n.comm <- nrow(spmat)
  disp.mat.weight <- comm.disp.mat(tdmat, spmat, abund = abund, sp.weighted = TRUE)
  M.beta <- sum(disp.mat.weight) / St^2
  M.beta.prime <- M.beta * n.comm / (n.comm - 1)
  fAB <- disp.mat.weight / sum(disp.mat.weight)
  if (isTRUE(all.equal(sum(disp.mat.weight), 0))) {
    Ht.beta <- 0
  } else if (q == 1) {
    fABlog <- ifelse(fAB == 0, 0, fAB * log(fAB))
    Ht.beta <- exp(-1 * sum(fABlog))
  } else if (q == 0) {
    Ht.beta <- sum(fAB > 0)
  } else {
    Ht.beta <- sum(fAB^q)^(1 / (1 - q))
  }
  qDT.beta <- (1 + sqrt(1 + 4 * Ht.beta)) / 2
  qDTM.beta <- 1 + qDT.beta * M.beta
  Et.beta <- qDT.beta / n.comm
  list(nsp = nsp, St = St, n.comm = n.comm, q = q, M.beta = M.beta, M.beta.prime = M.beta.prime, Ht.beta = Ht.beta, qDT.beta = qDT.beta, qDTM.beta = qDTM.beta, disp.mat.weight = disp.mat.weight)
}

# -------------------------------------------------------------------------------------
# Build Phylogeny via V.PhyloMaker2
# -------------------------------------------------------------------------------------

taxon_lookup <- taxon_lookup %>%
  mutate(
    derived_binomial = sapply(scientificName, function(x) {
      words <- unlist(strsplit(trimws(x), "\\s+"))
      if (length(words) >= 2) paste(words[1], words[2], sep = "_") else NA_character_
    }),
    phyloName = case_when(
      nchar(phyloName) <= 6 ~ coalesce(derived_binomial, phyloName),
      is.na(phyloName) ~ derived_binomial,
      TRUE ~ phyloName
    )
  )

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

# -------------------------------------------------------------------------------------
# Construct Top-Down Sensor-Visible Community Data Frame
# -------------------------------------------------------------------------------------

phylo_to_neon <- taxon_lookup %>%
  select(taxonID, phyloName) %>%
  distinct() %>%
  filter(!is.na(phyloName))

# Combine trees, shrubs, and ground cover using area_m2 / mean_pct (top-down footprint view)
topdown_keyed <- bind_rows(
  trees_df %>% select(siteID, plotID, taxonID, area_m2, source),
  shrubs_df %>% select(siteID, plotID, taxonID, area_m2, source),
  ground_df %>% rename(area_m2 = mean_pct) %>% select(siteID, plotID, taxonID, area_m2, source)
) %>%
  left_join(phylo_to_neon, by = "taxonID") %>%
  mutate(community_taxonID = coalesce(phyloName, taxonID)) %>%
  inner_join(plot_metadata %>% filter(plotType == "distributed") %>% select(plotID), by = "plotID")

neon_traits_site_keyed <- neon_traits_site %>%
  left_join(phylo_to_neon, by = "taxonID") %>%
  mutate(community_taxonID = coalesce(phyloName, taxonID))

neon_traits_cross_keyed <- neon_traits_fullCross %>%
  left_join(phylo_to_neon, by = "taxonID") %>%
  mutate(community_taxonID = coalesce(phyloName, taxonID))

trait_cols <- c("LMA_gm2", "N_pct", "C_pct")

# Filter sites meeting the >= 80% sensor-visible trait coverage threshold
functional_sites <- topdown_keyed %>%
  group_by(siteID) %>%
  summarize(
    pct_area_covered = 100 * sum(area_m2[community_taxonID %in% neon_traits_cross_keyed$community_taxonID], na.rm = TRUE) / sum(area_m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(pct_area_covered >= 80.0) %>%
  pull(siteID)

# Keep all sites that have valid community data for phylogenetic analysis
sites_to_keep <- topdown_keyed %>%
  distinct(siteID) %>%
  pull(siteID)

cat("Processing", length(sites_to_keep), "total sites (Functional metrics will be NA for sites below 80% trait coverage).\n")

# -------------------------------------------------------------------------------------
# 4. Top-Down qDTM Site-Level Execution Function (Conditional Functional Evaluation)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# 4. Top-Down qDTM Site-Level Execution Function (Conditional Functional Evaluation)
# -------------------------------------------------------------------------------------

run_topdown_qdtm_site <- function(site_df, site_id, trait_cols, phylo_tree = tre$scenario.3, q = 1, functional_allowed_sites = functional_sites) {
  empty_res <- data.frame(
    siteID = site_id,
    Alpha_FD_qDTM = NA_real_,
    Beta_FD_qDTM = NA_real_,
    Alpha_PD_qDTM = NA_real_,
    Beta_PD_qDTM = NA_real_,
    Alpha_Richness_FD = NA_real_,
    Total_Richness_FD = NA_real_,
    Alpha_Dispersion_FD = NA_real_,
    Beta_Dispersion_FD = NA_real_,
    Alpha_Evenness_FD = NA_real_,
    Beta_Evenness_FD = NA_real_,
    Alpha_Richness_PD = NA_real_,
    Total_Richness_PD = NA_real_,
    Alpha_Dispersion_PD = NA_real_,
    Beta_Dispersion_PD = NA_real_,
    Alpha_Evenness_PD = NA_real_,
    Beta_Evenness_PD = NA_real_
  )
  
  # --- Base Community Preparation ---
  comm_raw_base <- site_df %>%
    group_by(plotID, community_taxonID) %>%
    summarize(total_area = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = community_taxonID, values_from = total_area, values_fill = 0) %>%
    column_to_rownames("plotID")
  
  comm_filtered_base <- comm_raw_base[rowSums(comm_raw_base > 0) >= 2, , drop = FALSE]
  if (nrow(comm_filtered_base) < 2) return(empty_res)
  
  # --- Functional qDTM & Components (Only if site meets >= 80% trait coverage) ---
  alpha_fd_qdtm <- beta_fd_qdtm <- alpha_rich_fd <- total_rich_fd <- alpha_disp_fd <- beta_disp_fd <- alpha_even_fd <- beta_even_fd <- NA_real_
  
  if (site_id %in% functional_allowed_sites) {
    site_traits_raw <- neon_traits_cross_keyed %>%
      select(community_taxonID, all_of(trait_cols)) %>%
      left_join(
        neon_traits_site_keyed %>%
          filter(siteID == site_id) %>%
          select(community_taxonID, all_of(trait_cols)),
        by = "community_taxonID",
        suffix = c("_global", "_site"),
        relationship = "many-to-many"
      ) %>%
      mutate(
        LMA_gm2         = coalesce(LMA_gm2_site,         LMA_gm2_global),
        N_pct           = coalesce(N_pct_site,           N_pct_global),
        C_pct           = coalesce(C_pct_site,           C_pct_global)
      ) %>%
      select(community_taxonID, all_of(trait_cols)) %>%
      drop_na() %>%
      distinct(community_taxonID, .keep_all = TRUE)
    
    if (nrow(site_traits_raw) >= 2) {
      traits_mat <- site_traits_raw %>% column_to_rownames("community_taxonID")
      comm_func_raw <- comm_raw_base[, colnames(comm_raw_base) %in% rownames(traits_mat), drop = FALSE]
      comm_func_filtered <- comm_func_raw[rowSums(comm_func_raw > 0) >= 2, , drop = FALSE]
      
      if (nrow(comm_func_filtered) >= 2) {
        species_present_func <- colnames(comm_func_filtered)[colSums(comm_func_filtered) > 0]
        if (length(species_present_func) >= 2) {
          comm_func_mat <- comm_func_filtered[, species_present_func, drop = FALSE]
          comm_func_mat <- comm_func_mat / rowSums(comm_func_mat)
          
          site_traits_mat <- traits_mat[species_present_func, , drop = FALSE]
          site_traits_scaled <- scale(site_traits_mat)
          trait_dist <- as.matrix(dist(site_traits_scaled, method = "euclidean"))
          trait_dist_norm <- (trait_dist - min(trait_dist)) / (max(trait_dist) - min(trait_dist))
          
          func_comm_res <- tryCatch({ FTD.comm(trait_dist_norm, comm_func_mat, q = q, abund = TRUE, match.names = TRUE) }, error = function(e) NULL)
          func_beta_res <- tryCatch({ FTD.beta(trait_dist_norm, comm_func_mat, abund = TRUE, q = q) }, error = function(e) NULL)
          
          alpha_fd_qdtm  <- if (!is.null(func_comm_res)) func_comm_res$u.qDTM else NA_real_
          beta_fd_qdtm   <- if (!is.null(func_beta_res)) func_beta_res$qDTM.beta else NA_real_
          alpha_rich_fd  <- if (!is.null(func_comm_res)) func_comm_res$u.nsp else NA_real_
          total_rich_fd  <- if (!is.null(func_beta_res)) func_beta_res$St else NA_real_
          alpha_disp_fd  <- if (!is.null(func_comm_res)) func_comm_res$u.M else NA_real_
          beta_disp_fd   <- if (!is.null(func_beta_res)) func_beta_res$M.beta else NA_real_
          alpha_even_fd  <- if (!is.null(func_comm_res)) func_comm_res$u.qEt else NA_real_
          beta_even_fd   <- if (!is.null(func_beta_res)) func_beta_res$qDT.beta / func_beta_res$n.comm else NA_real_
        }
      }
    }
  }
  
  # --- Phylogenetic qDTM & Components (Evaluated for all sites independently) ---
  alpha_pd_qdtm <- beta_pd_qdtm <- alpha_rich_pd <- total_rich_pd <- alpha_disp_pd <- beta_disp_pd <- alpha_even_pd <- beta_even_pd <- NA_real_
  
  if (!is.null(phylo_tree)) {
    species_present_comm <- colnames(comm_filtered_base)[colSums(comm_filtered_base) > 0]
    comm_comm_base <- comm_filtered_base[, species_present_comm, drop = FALSE]
    
    matched_phylo <- match.phylo.comm(phylo_tree, comm_comm_base)
    if (!is.null(matched_phylo) && length(matched_phylo$phy$tip.label) >= 2) {
      valid_plots <- which(rowSums(matched_phylo$comm > 0) > 0)
      if (length(valid_plots) >= 2) {
        clean_phylo_comm <- matched_phylo$comm[valid_plots, , drop = FALSE]
        
        phylo_dist <- cophenetic(matched_phylo$phy)
        phylo_dist_norm <- (phylo_dist - min(phylo_dist)) / (max(phylo_dist) - min(phylo_dist))
        
        comm_phylo <- clean_phylo_comm / rowSums(clean_phylo_comm)
        
        phylo_comm_res <- tryCatch({ FTD.comm(phylo_dist_norm, comm_phylo, q = q, abund = TRUE, match.names = TRUE) }, error = function(e) NULL)
        phylo_beta_res <- tryCatch({ FTD.beta(phylo_dist_norm, comm_phylo, abund = TRUE, q = q) }, error = function(e) NULL)
        
        alpha_pd_qdtm  <- if (!is.null(phylo_comm_res)) phylo_comm_res$u.qDTM else NA_real_
        beta_pd_qdtm   <- if (!is.null(phylo_beta_res)) phylo_beta_res$qDTM.beta else NA_real_
        alpha_rich_pd  <- if (!is.null(phylo_comm_res)) phylo_comm_res$u.nsp else NA_real_
        total_rich_pd  <- if (!is.null(phylo_beta_res)) phylo_beta_res$St else NA_real_
        alpha_disp_pd  <- if (!is.null(phylo_comm_res)) phylo_comm_res$u.M else NA_real_
        beta_disp_pd   <- if (!is.null(phylo_beta_res)) phylo_beta_res$M.beta else NA_real_
        alpha_even_pd  <- if (!is.null(phylo_comm_res)) phylo_comm_res$u.qEt else NA_real_
        beta_even_pd   <- if (!is.null(phylo_beta_res)) phylo_beta_res$qDT.beta / phylo_beta_res$n.comm else NA_real_
      }
    }
  }
  
  data.frame(
    siteID = site_id,
    Alpha_FD_qDTM = alpha_fd_qdtm,
    Beta_FD_qDTM = beta_fd_qdtm,
    Alpha_PD_qDTM = alpha_pd_qdtm,
    Beta_PD_qDTM = beta_pd_qdtm,
    Alpha_Richness_FD = alpha_rich_fd,
    Total_Richness_FD = total_rich_fd,
    Alpha_Dispersion_FD = alpha_disp_fd,
    Beta_Dispersion_FD = beta_disp_fd,
    Alpha_Evenness_FD = alpha_even_fd,
    Beta_Evenness_FD = beta_even_fd,
    Alpha_Richness_PD = alpha_rich_pd,
    Total_Richness_PD = total_rich_pd,
    Alpha_Dispersion_PD = alpha_disp_pd,
    Beta_Dispersion_PD = beta_disp_pd,
    Alpha_Evenness_PD = alpha_even_pd,
    Beta_Evenness_PD = beta_even_pd
  )
}

# -------------------------------------------------------------------------------------
# Save Top-Down qDTM Results
# -------------------------------------------------------------------------------------

topdown_qdtm_master <- map_dfr(sites_to_keep, function(s) {
  site_df <- topdown_keyed %>% filter(siteID == s)
  run_topdown_qdtm_site(site_df, site_id = s, trait_cols = trait_cols, phylo_tree = tre$scenario.3, q = 1)
}) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

write.csv(topdown_qdtm_master, "NEON_PlantDiv/new_data_out/NEON_TopDown_qDTM_Results_2024.csv", row.names = FALSE)
save(topdown_qdtm_master, tre, file = "NEON_PlantDiv/new_data_out/topdown_qdtm_results_2024.RData")



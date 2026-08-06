# 01_EMIT_Run: Script for processing tiles and calculating specdiv for EMIT sites

library(tidyverse)
library(ncdf4)
library(adespatial)
library(terra)

# set wd
setwd("~/EMIT_specdiv")

# pull in the specdiv functions
source("scripts/00_Specdiv_Func.R")

# Create a local, safe directory for temp files and register it with terra
temp_dir <- file.path(getwd(), "data_work/r_temp")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = temp_dir, memfrac = 0.8)

# Base for year of choice
base_path <- "data_in/2024_samplingBoundary"

# Set up output paths (and ensure results folder exists)
results_csv <- "data_out/results/EMIT_SpecDiv_ssmplingBoundary_2023.csv"
dir.create(dirname(results_csv), recursive = TRUE, showWarnings = FALSE)
dir.create("data_out/results_2024", recursive = TRUE, showWarnings = FALSE)
dir.create("data_out/maps_2024", recursive = TRUE, showWarnings = FALSE)
dir.create("data_work/new_cubes_2024", recursive = TRUE, showWarnings = FALSE)

# Get all site directories
site_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
results <- list()

# run through all: 
for (site_dir in site_dirs) {
  site <- basename(site_dir)
  
  nc_files <- list.files(site_dir, pattern = "\\.nc$", full.names = TRUE)
  num_files <- length(nc_files)
  
  if (num_files == 0) {
    message("Nothin there for: ", site)
    next
  }
  
  for (f_idx in 1:num_files) {
    nc_file <- nc_files[f_idx]
    file_name <- basename(nc_file)
    site_suffix <- if (num_files > 1) paste0(site, "_", f_idx) else site
    
    message("Processing NetCDF and calculating div for: ", site_suffix)
    
    tryCatch({
      gc()
      dat <- nc_open(nc_file)
      wavelengths <- ncvar_get(dat, "wavelengths")
      hyperspectral_cube <- ncvar_get(dat, "reflectance") # dimensions: [lon, lat, bands]
      nc_close(dat)
      
      # Masking out atmospheric water vapor absorption bands
      wvl_df <- data.frame(index = 1:length(wavelengths), wvl = wavelengths)
      good_indices <- wvl_df %>%
        filter(wvl > 400 & wvl < 2450) %>%
        filter(!(wvl > 1320 & wvl < 1440)) %>%
        filter(!(wvl > 1770 & wvl < 1970)) %>%
        pull(index)
      
      sub_cube <- hyperspectral_cube[, , good_indices]
      sub_cube[sub_cube < 0] <- 0
      
      # Permute dimensions from [lon, lat, bands] to terra's expected [lat, lon, bands] layout
      arr_permuted <- aperm(sub_cube, c(2, 1, 3))
      
      # Build SpatRaster directly from the 3D array in memory
      cube <- rast(arr_permuted)
      
      # Run diversity pipeline directly on the in-memory SpatRaster
      cube_norm <- bright_norm(cube)
      
      disp_res <- beta_diversity_dispersion(cube_norm)
      PW_res <- beta_diversity_pairwise(cube_norm) 
      agg_res <- beta_diversity_dispersion_aggregated(cube_norm, agg_factor = 10)
      agg_PW_res <- beta_diversity_pairwise_aggregated(cube_norm, agg_factor = 10)
      
      site_row <- data.frame(
        site = site_suffix,
        original_file = file_name,
        sum_squares = disp_res$sum_squares,
        Beta_dispersion = disp_res$beta_dispersion,
        Beta_avg_pairwise = PW_res$avg_beta_pairwise,
        agg_factor_used = agg_res$agg_factor_used,
        gamma_sum_squares = agg_res$gamma_sum_squares,
        gamma_dispersion = agg_res$gamma_dispersion,
        beta_agg_sum_squares = agg_res$beta_agg_sum_squares,
        beta_agg_dispersion = agg_res$beta_agg_dispersion,
        Beta_agg_avg_pairwise = agg_PW_res$avg_beta_pairwise,
        n_aggregate_pixels = agg_res$n_aggregate_pixels,
        timestamp = Sys.time()
      )
      
      write.table(site_row, results_csv, 
                  append = file.exists(results_csv), 
                  sep = ",", row.names = FALSE, 
                  col.names = !file.exists(results_csv))
      
      # Save maps out safely
      if (!is.null(PW_res$lcbd_map)) {
        terra::writeRaster(PW_res$lcbd_map, file.path("data_out/maps_2023", paste0(site_suffix, "_uniqueness.tif")), overwrite = TRUE)
      }
      if (!is.null(agg_res$lcbd_agg_map)) {
        terra::writeRaster(agg_res$lcbd_agg_map, file.path("data_out/maps_2023", paste0(site_suffix, "_agg_uniqueness.tif")), overwrite = TRUE)
      }
      if (!is.null(agg_PW_res$lcbd_map)) {
        terra::writeRaster(agg_PW_res$lcbd_map, file.path("data_out/maps_2023", paste0(site_suffix, "_agg_pw_uniqueness.tif")), overwrite = TRUE)
      }
      
      results[[site_suffix]] <- list(disp_res = disp_res, PW_res = PW_res, agg_res = agg_res, agg_PW_res = agg_PW_res)
      message("SUCCESS for: ", site_suffix)
      
      rm(cube, cube_norm, sub_cube, arr_permuted, disp_res, PW_res, agg_res, agg_PW_res, site_row)
      gc() 
      
    }, error = function(e) {
      message("ERROR processing diversity metrics for ", site_suffix, ": ", e$message)
    })
  }
}

# STEP 3: Save final combined table
if (length(results) > 0) {
  results_df <- bind_rows(lapply(results, function(x) {
    data.frame(
      sum_squares = x$disp_res$sum_squares,
      Beta_dispersion = x$disp_res$beta_dispersion,
      Beta_avg_pairwise = x$PW_res$avg_beta_pairwise,
      gamma_sum_squares = x$agg_res$gamma_sum_squares,
      gamma_dispersion = x$agg_res$gamma_dispersion,
      beta_agg_sum_squares = x$agg_res$beta_agg_sum_squares,
      beta_agg_dispersion = x$agg_res$beta_agg_dispersion,
      Beta_agg_avg_pairwise = x$agg_PW_res$avg_beta_pairwise
    )
  }), .id = "site")
  
  write_csv(results_df, "data_out/results_2024/EMIT_SpecDiv_2024.csv")
}


# 01_EMIT_Run: Script for processing tiles and calculating specdiv for EMIT sites

library(sp)
library(raster)
library(tidyverse)
library(ncdf4)
library(adespatial)

# Download all EMIT tiles through AppEEARS
  # Download the sampling boundaries through NEON's data portal and upload to AppEEARS
  # Set dates etc. 

# set wd
setwd(dir = "EMIT_SpecDiv/")

# pull in the specdiv functions
source("scripts/00_Specdiv_Func.R")


# Base for year of choice
base_path <- "data_in/2023_samplingBoundary"

# Set up output paths
results_csv <- "data_out/results/EMIT_SpecDiv_ssmplingBoundary_2023.csv"
dir.create("data_out/results_2023", recursive = TRUE, showWarnings = FALSE)
dir.create("data_out/maps_2023", recursive = TRUE, showWarnings = FALSE)
dir.create("data_work/new_cubes_2023", recursive = TRUE, showWarnings = FALSE)

# Get all site directories
site_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
results <- list()

# run through all: 
for (site_dir in site_dirs) {
  site <- basename(site_dir)
  
  # list the nc files
  nc_files <- list.files(site_dir, pattern = "\\.nc$", full.names = TRUE)
  num_files <- length(nc_files)
  
  # catch error
  if (num_files == 0) {
    message("Nothin there for: ", site)
    next
  }
  
  # Loop through each file for the site 
  for (f_idx in 1:num_files) {
    nc_file <- nc_files[f_idx]
    file_name <- basename(nc_file)
    
    # If there's more than one file, append _1, _2, etc. 
    site_suffix <- if (num_files > 1) paste0(site, "_", f_idx) else site
    
    out_tif <- file.path("data_work/new_cubes_2023", paste0(site_suffix, "_hyperspectral_cube.tif"))
    
    # STEP 1! Process NetCDF to GeoTIFF
    if (!file.exists(out_tif)) {
      message("Processing NetCDF for: ", site_suffix)
      
      tryCatch({
        dat <- nc_open(nc_file)
        wavelengths <- ncvar_get(dat, "wavelengths")
        hyperspectral_cube <- ncvar_get(dat, "reflectance")
        
        # Masking out atmospheric water vapor absorption bands and ends
        wvl_df <- data.frame(index = 1:length(wavelengths), wvl = wavelengths)
        good_indices <- wvl_df %>%
          filter(wvl > 400 & wvl < 2450) %>%
          filter(!(wvl > 1320 & wvl < 1440)) %>%
          filter(!(wvl > 1770 & wvl < 1970)) %>%
          pull(index)
        
        # Build stack with good bands
        hyperspectral_stack <- stack(lapply(good_indices, function(b) {
          r <- raster(hyperspectral_cube[, , b])
          r[r < 0] <- 0   # Force negative artifacts to 0
          return(r)
        }))
        
        writeRaster(hyperspectral_stack, out_tif, format = "GTiff", overwrite = TRUE)
        nc_close(dat)
        message("Saved tif with ", length(good_indices), " bands.")
        
      }, error = function(e) { 
        message("Skipping conversion for ", site_suffix, ": ", e$message) 
      })
    }
    
    # STEP 2! Calculate Spectral Turnover & Dispersion
    if (file.exists(out_tif)) {
      message("Calculating div for ", site_suffix)
      
      tryCatch({
        cube <- readAll(brick(out_tif)) 
        cube_norm <- bright_norm(cube)
        
        # Calculate dispersion, pairwise beta diversity, aggregated dispersion, and aggregated pairwise metrics
        disp_res <- beta_diversity_dispersion(cube_norm)
        PW_res <- beta_diversity_pairwise(cube_norm) 
        agg_res <- beta_diversity_aggregated(cube_norm, agg_factor = 10)
        agg_PW_res <- beta_diversity_pairwise_aggregated(cube_norm, agg_factor = 10)
        
        # Create a row for this specific tile
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
        
        # Append rows to the primary output CSV as it goes
        write.table(site_row, results_csv, 
                    append = file.exists(results_csv), 
                    sep = ",", row.names = FALSE, 
                    col.names = !file.exists(results_csv))
        
        # Save LCBD map if available
        if (!is.null(PW_res$lcbd_map)) {
          writeRaster(PW_res$lcbd_map, 
                      file.path("data_out/maps_2023", paste0(site_suffix, "_uniqueness.tif")), 
                      overwrite = TRUE)
          message("-> Saved uniqueness map for: ", site_suffix)
        }
        
        # Save aggregated LCBD map if available
        if (!is.null(agg_res$lcbd_agg_map)) {
          writeRaster(agg_res$lcbd_agg_map, 
                      file.path("data_out/maps_2023", paste0(site_suffix, "_agg_uniqueness.tif")), 
                      overwrite = TRUE)
          message("-> Saved aggregated uniqueness map for: ", site_suffix)
        }
        
        # Save aggregated pairwise LCBD map if available
        if (!is.null(agg_PW_res$lcbd_map)) {
          writeRaster(agg_PW_res$lcbd_map, 
                      file.path("data_out/maps_2023", paste0(site_suffix, "_agg_pw_uniqueness.tif")), 
                      overwrite = TRUE)
          message("-> Saved aggregated pairwise uniqueness map for: ", site_suffix)
        }
        
        # Store structured list for final backup df export
        results[[site_suffix]] <- list(disp_res = disp_res, PW_res = PW_res, agg_res = agg_res, agg_PW_res = agg_PW_res)
        message("SUCCESS for: ", site_suffix)
        
        # Clean up memory explicitly per iteration
        rm(cube, cube_norm, disp_res, PW_res, agg_res, agg_PW_res, site_row)
        gc() 
        
      }, error = function(e) {
        message("ERROR processing diversity metrics for ", site_suffix, ": ", e$message)
      })
    }
  }
}

# STEP 3! Save it bb
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
  
  write_csv(results_df, "data_out/results/EMIT_SpecDiv_2023.csv")
}
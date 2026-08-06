## 00_Specdiv_Func: terra-optimized functions for EMIT work

# --- Brightness normalization for raster cubes ---
bright_norm <- function(cube) {
  library(terra)
  
  mat <- values(cube, mat = TRUE)
  mat[is.na(mat)] <- 0
  
  pixel_brightness <- sqrt(rowSums(mat^2, na.rm = TRUE))
  mat_norm <- sweep(mat, 1, pixel_brightness, "/")
  mat_norm[is.nan(mat_norm) | is.infinite(mat_norm)] <- NA
  
  cube_norm <- rast(cube)
  values(cube_norm) <- mat_norm
  
  return(cube_norm)
}

# --- Sum of Squares ---
sum_squares <- function(Y) {
  Y <- as.matrix(Y)   
  n <- nrow(Y)
  if (n < 2) return(list(ss = NA, sdiv = NA, fcsd = rep(NA, ncol(Y))))
  
  Y.cent <- scale(Y, center = TRUE, scale = FALSE)
  sij <- Y.cent^2
  SS.total <- sum(sij, na.rm = TRUE)
  SS.col <- colSums(sij, na.rm = TRUE)
  fcsd <- SS.col / SS.total
  sdiv <- SS.total / (n - 1)
  list(ss = SS.total, sdiv = sdiv, fcsd = fcsd)
}

# --- Dispersion Beta Diversity Function ---
beta_diversity_dispersion <- function(cube) {
  library(terra)
  library(tidyverse)
  
  mat <- values(cube, mat = TRUE)
  mat <- na.omit(mat)
  
  sdiv_gamma <- sum_squares(mat)
  
  output <- list(
    sum_squares = sdiv_gamma$ss,
    beta_dispersion = sdiv_gamma$sdiv,
    beta_disp_fcsd = sdiv_gamma$fcsd
  )
  return(output)
}

# --- Pairwise Beta Diversity Function ---
beta_diversity_pairwise <- function(cube) {
  library(terra)
  library(vegan)
  library(adespatial)
  
  xy_full <- terra::xyFromCell(cube, 1:ncell(cube))
  mat_full <- values(cube, mat = TRUE)
  
  valid_idx <- complete.cases(mat_full)
  xy_full <- xy_full[valid_idx, ]
  mat_full <- mat_full[valid_idx, ]
  mat_full[mat_full < 0] <- 0
  n_pixels <- nrow(mat_full)
  
  col_means <- colMeans(mat_full)
  sq_dist <- rowSums(sweep(mat_full, 2, col_means)^2)
  lcbd_vals <- sq_dist / sum(sq_dist)
  
  lcbd_df <- data.frame(x = xy_full[,1], y = xy_full[,2], lcbd = lcbd_vals)
  lcbd_raster <- rast(lcbd_df, type = "xyz", crs = crs(cube))
  
  if (n_pixels > 30000) {
    message("huge site.. will take longer bc averaging")
    skip <- ceiling(sqrt(n_pixels / 30000))
    rows_to_keep <- seq(1, n_pixels, by = skip)
    mat_sub <- mat_full[rows_to_keep, ]
  } else {
    mat_sub <- mat_full
  }
  
  dist_matrix <- vegan::vegdist(mat_sub, method = "euclidean")
  avg_beta <- mean(dist_matrix)
  
  return(list(avg_beta_pairwise = avg_beta, lcbd_map = lcbd_raster))
}

# --- Aggregated Spectral Diversity & Variance Partitioning ---
beta_diversity_dispersion_aggregated <- function(cube, agg_factor = 10) {
  library(terra)
  library(tidyverse)
  
  ncols <- ncol(cube)
  nrows <- nrow(cube)
  
  if (ncols < agg_factor || nrows < agg_factor) {
    message("site too small for agg_factor = ", agg_factor, ".. scaling down")
    agg_factor <- min(floor(ncols / 2), floor(nrows / 2))
    if (agg_factor < 2) {
      stop("site is wayyyy to small... cannot aggregate")
    }
  }
  
  mat_gamma <- na.omit(values(cube, mat = TRUE))
  gamma_res <- sum_squares(mat_gamma)
  
  cube_agg <- terra::aggregate(cube, fact = agg_factor, fun = mean, na.rm = TRUE)
  mat_beta <- na.omit(values(cube_agg, mat = TRUE))
  
  beta_res <- sum_squares(mat_beta)
  
  xy_agg <- terra::xyFromCell(cube_agg, 1:ncell(cube_agg))
  valid_agg <- complete.cases(values(cube_agg, mat = TRUE))
  col_means_agg <- colMeans(mat_beta)
  sq_dist_agg <- rowSums(sweep(mat_beta, 2, col_means_agg)^2)
  lcbd_agg_vals <- sq_dist_agg / sum(sq_dist_agg)
  
  lcbd_agg_df <- data.frame(x = xy_agg[valid_agg, 1], y = xy_agg[valid_agg, 2], lcbd = lcbd_agg_vals)
  lcbd_agg_raster <- rast(lcbd_agg_df, type = "xyz", crs = crs(cube))
  
  output <- list(
    agg_factor_used      = agg_factor,
    gamma_sum_squares    = gamma_res$ss,
    gamma_dispersion     = gamma_res$sdiv,
    beta_agg_sum_squares = beta_res$ss,
    beta_agg_dispersion  = beta_res$sdiv,
    beta_agg_fcsd        = beta_res$fcsd,
    n_aggregate_pixels   = nrow(mat_beta),
    lcbd_agg_map         = lcbd_agg_raster
  )
  
  return(output)
}

# --- Pairwise Beta Diversity Function for Aggregates ---
beta_diversity_pairwise_aggregated <- function(cube, agg_factor = 10) {
  library(terra)
  library(vegan)
  library(adespatial)
  
  ncols <- ncol(cube)
  nrows <- nrow(cube)
  
  if (ncols < agg_factor || nrows < agg_factor) {
    message("site too small for agg_factor = ", agg_factor, ".. scaling down")
    agg_factor <- min(floor(ncols / 2), floor(nrows / 2))
    if (agg_factor < 2) {
      stop("site is wayyyy to small... cannot aggregate")
    }
  }
  
  cube_agg <- terra::aggregate(cube, fact = agg_factor, fun = mean, na.rm = TRUE)
  xy_full <- terra::xyFromCell(cube_agg, 1:ncell(cube_agg))
  mat_full <- values(cube_agg, mat = TRUE)
  
  valid_idx <- complete.cases(mat_full)
  xy_full <- xy_full[valid_idx, ]
  mat_full <- mat_full[valid_idx, ]
  mat_full[mat_full < 0] <- 0
  n_pixels <- nrow(mat_full)
  
  col_means <- colMeans(mat_full)
  sq_dist <- rowSums(sweep(mat_full, 2, col_means)^2)
  lcbd_vals <- sq_dist / sum(sq_dist)
  
  lcbd_df <- data.frame(x = xy_full[,1], y = xy_full[,2], lcbd = lcbd_vals)
  lcbd_raster <- rast(lcbd_df, type = "xyz", crs = crs(cube))
  
  if (n_pixels > 30000) {
    message("huge aggregate site.. will take longer bc averaging")
    skip <- ceiling(sqrt(n_pixels / 30000))
    rows_to_keep <- seq(1, n_pixels, by = skip)
    mat_sub <- mat_full[rows_to_keep, ]
  } else {
    mat_sub <- mat_full
  }
  
  dist_matrix <- vegan::vegdist(mat_sub, method = "euclidean")
  avg_beta <- mean(dist_matrix)
  
  return(list(
    agg_factor_used    = agg_factor,
    avg_beta_pairwise  = avg_beta, 
    lcbd_map           = lcbd_raster,
    n_aggregate_pixels = n_pixels
  ))
}
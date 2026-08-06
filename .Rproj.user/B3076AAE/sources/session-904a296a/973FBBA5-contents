## 00_Specdiv_Func: functions for EMIT work
# adapted from Ettienne's specdiv functions
# source this code in scriot 01_EMIT_Run

# --- Brightness normalization for raster cubes ---
bright_norm <- function(cube) {
  require(raster)
  
  mat <- as.matrix(cube)
  
  pixel_brightness <- sqrt(rowSums(mat^2, na.rm = TRUE))
  mat_norm <- sweep(mat, 1, pixel_brightness, "/")
  mat_norm[is.nan(mat_norm) | is.infinite(mat_norm)] <- NA
  cube_norm <- brick(cube)
  values(cube_norm) <- mat_norm
  
  return(cube_norm)
}

# --- Sum of Squares ---
sum_squares <- function(Y) {
  Y <- as.matrix(Y)   # force numeric matrix
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

# --- Dispersin Beta Diversity Function ---
beta_diversity_dispersion <- function(cube) {
  require(raster)
  require(tidyverse)
  
  cube_points <- rasterToPoints(cube, spatial = FALSE)
  
  if(is.list(cube_points)) { cube_points <- do.call(rbind, cube_points) }
  cube_points <- as_tibble(cube_points)
  
  value_columns <- colnames(cube_points)[3:ncol(cube_points)]
  cube_points_sel <- cube_points %>% dplyr::select(all_of(value_columns))
  
  cube_mat <- as.matrix(cube_points_sel)
  sdiv_gamma <- sum_squares(cube_mat)
  
  output <- list(
    sum_squares = sdiv_gamma$ss,
    beta_dispersion = sdiv_gamma$sdiv,
    beta_disp_fcsd = sdiv_gamma$fcsd
  )
  return(output)
}

# --- Pairwise Beta Diversity Function ---
beta_diversity_pairwise <- function(cube) {
  require(raster)
  require(vegan)
  require(adespatial)
  
  cube_points_full <- rasterToPoints(cube, spatial = FALSE)
  
  # force list back for large sites
  if(is.list(cube_points_full)) { cube_points_full <- do.call(rbind, cube_points_full) }
  
  xy_full <- cube_points_full[, 1:2]
  mat_full <- cube_points_full[, 3:ncol(cube_points_full)]
  mat_full[mat_full < 0] <- 0
  n_pixels <- nrow(mat_full)
  
  # calc map
  col_means <- colMeans(mat_full)
  sq_dist <- rowSums(sweep(mat_full, 2, col_means)^2)
  lcbd_vals <- sq_dist / sum(sq_dist)
  
  lcbd_df <- data.frame(x = xy_full[,1], y = xy_full[,2], lcbd = lcbd_vals)
  lcbd_raster <- rasterFromXYZ(lcbd_df, crs = proj4string(cube))
  
  # calc beta
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
  require(raster)
  require(tidyverse)
  
  ncols <- ncol(cube)
  nrows <- nrow(cube)
  
  if (ncols < agg_factor || nrows < agg_factor) {
    message("site too small for agg_factor = ", agg_factor, ".. scaling down")
    agg_factor <- min(floor(ncols / 2), floor(nrows / 2))
    if (agg_factor < 2) {
      stop("site is wayyyy to small... cannot aggregate")
    }
  }
  
  cube_points <- rasterToPoints(cube, spatial = FALSE)
  if(is.list(cube_points)) { cube_points <- do.call(rbind, cube_points) }
  cube_points <- as_tibble(cube_points)
  
  value_columns <- colnames(cube_points)[3:ncol(cube_points)]
  mat_gamma <- as.matrix(cube_points %>% dplyr::select(all_of(value_columns)))
  
  gamma_res <- sum_squares(mat_gamma)
  
  cube_agg <- raster::aggregate(cube, fact = agg_factor, fun = mean, na.rm = TRUE)
  
  agg_points <- rasterToPoints(cube_agg, spatial = FALSE)
  if(is.list(agg_points)) { agg_points <- do.call(rbind, agg_points) }
  agg_points <- as_tibble(agg_points)
  
  mat_beta <- as.matrix(agg_points %>% dplyr::select(all_of(value_columns)))
  
  beta_res <- sum_squares(mat_beta)
  
  output <- list(
    agg_factor_used      = agg_factor,
    gamma_sum_squares    = gamma_res$ss,
    gamma_dispersion     = gamma_res$sdiv,
    beta_agg_sum_squares = beta_res$ss,
    beta_agg_dispersion  = beta_res$sdiv,
    beta_agg_fcsd        = beta_res$fcsd,
    n_aggregate_pixels   = nrow(mat_beta),
    lcbd_agg_map         = rasterFromXYZ(agg_points[, 1:2], crs = proj4string(cube))
  )
  
  return(output)
}

# --- Pairwise Beta Diversity Function for Aggregates ---
beta_diversity_pairwise_aggregated <- function(cube, agg_factor = 10) {
  require(raster)
  require(vegan)
  require(adespatial)
  
  ncols <- ncol(cube)
  nrows <- nrow(cube)
  
  if (ncols < agg_factor || nrows < agg_factor) {
    message("site too small for agg_factor = ", agg_factor, ".. scaling down")
    agg_factor <- min(floor(ncols / 2), floor(nrows / 2))
    if (agg_factor < 2) {
      stop("site is wayyyy to small... cannot aggregate")
    }
  }
  
  cube_agg <- raster::aggregate(cube, fact = agg_factor, fun = mean, na.rm = TRUE)
  cube_points_full <- rasterToPoints(cube_agg, spatial = FALSE)
  
  if(is.list(cube_points_full)) { cube_points_full <- do.call(rbind, cube_points_full) }
  
  xy_full <- cube_points_full[, 1:2]
  mat_full <- cube_points_full[, 3:ncol(cube_points_full)]
  mat_full[mat_full < 0] <- 0
  n_pixels <- nrow(mat_full)
  
  col_means <- colMeans(mat_full)
  sq_dist <- rowSums(sweep(mat_full, 2, col_means)^2)
  lcbd_vals <- sq_dist / sum(sq_dist)
  
  lcbd_df <- data.frame(x = xy_full[,1], y = xy_full[,2], lcbd = lcbd_vals)
  lcbd_raster <- rasterFromXYZ(lcbd_df, crs = proj4string(cube))
  
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
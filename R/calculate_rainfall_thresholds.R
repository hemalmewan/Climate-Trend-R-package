#' Calculate Pixel-wise Rainfall Thresholds (Equal Cumulative Volume)
#'
#' @description
#' Calculates the "Low" and "High" daily rainfall thresholds that divide the 
#' total cumulative rainfall volume into three equal parts over the entire time series.
#'
#' @details
#' This function dynamically calculates intensity thresholds for each pixel. 
#' It filters for wet days (> 1 mm), sorts the daily rainfall amounts, and calculates 
#' the cumulative sum. 
#' 
#' \strong{Mathematical Definition:}
#' Let \eqn{P} be the time series of daily precipitation for a single pixel. 
#' Strictly positive wet days are extracted and sorted into an ascending sequence \eqn{S}:
#' \deqn{S = \{s_1, s_2, \dots, s_n\} \quad \text{where } s_i > \text{threshold}}
#' 
#' The cumulative volume \eqn{C_k} at index \eqn{k} is the sum of all sorted rainfall up to that point:
#' \deqn{C_k = \sum_{i=1}^{k} s_i}
#' 
#' Let the total cumulative volume \eqn{V = C_n}. The Low Threshold (\eqn{T_{low}}) and 
#' High Threshold (\eqn{T_{high}}) are defined as the specific daily rainfall amounts 
#' \eqn{s_k} where the cumulative volume first crosses \eqn{\frac{V}{3}} and \eqn{\frac{2V}{3}}:
#' \deqn{T_{low} = s_a \quad \text{where } a = \min \left\{ k \mid C_k \ge \frac{V}{3} \right\}}
#' \deqn{T_{high} = s_b \quad \text{where } b = \min \left\{ k \mid C_k \ge \frac{2V}{3} \right\}}
#'
#' @param nc_files Character vector. Full file paths to the year-wise NetCDF files.
#' @param threshold Numeric. Daily rainfall threshold to define a wet day. Default is 1.0 mm.
#' @param cores Numeric. Number of CPU cores to use for parallel processing. Default is 1.
#' @param save_path Character. (Optional) File path to save the resulting 2-layer raster (.tif or .nc).
#'
#' @return A SpatRaster with two layers: "Low_Threshold" and "High_Threshold".
#' @export
#' @importFrom terra rast app writeRaster
calculate_rainfall_thresholds <- function(nc_files, 
                                          threshold = 1.0, 
                                          cores = 1, 
                                          save_path = NULL) {
  
  message("Virtually stacking NetCDF files (Lazy Loading)...")
  r_stack <- terra::rast(nc_files)
  
  message("Applying cumulative volume threshold math across all pixels...")
  
  # --- THE CORE MATH FUNCTION ---
  core_threshold_math <- function(x) {
    
    # Extract strictly positive wet days
    x_wet <- x[!is.na(x) & x > threshold]
    
    # Return NA if there is no valid rain data
    if (length(x_wet) == 0) return(c(Low_Threshold = NA, High_Threshold = NA))
    
    # Sort ascending and calculate cumulative sum
    x_sorted <- base::sort(x_wet)
    cum_vol <- base::cumsum(x_sorted)
    
    # Total volume
    total_vol <- cum_vol[length(cum_vol)]
    
    # Define volume targets
    target1 <- total_vol / 3
    target2 <- (total_vol * 2) / 3
    
    # Find the indices that cross the targets
    idx_low <- which(cum_vol >= target1)[1]
    idx_high <- which(cum_vol >= target2)[1]
    
    low_thresh <- x_sorted[idx_low]
    high_thresh <- x_sorted[idx_high]
    
    return(c(Low_Threshold = low_thresh, High_Threshold = high_thresh))
  }
  
  # --- PARALLEL SPATIAL EXECUTION ---
  threshold_maps <- terra::app(r_stack, fun = core_threshold_math, cores = cores)
  names(threshold_maps) <- c("Low_Threshold", "High_Threshold")
  
  # --- SAVE OUTPUT ---
  if (!is.null(save_path)) {
    message(paste("Saving threshold map to:", save_path))
    terra::writeRaster(threshold_maps, filename = save_path, overwrite = TRUE)
  }
  
  message("✅ Threshold calculation complete!")
  return(threshold_maps)
}
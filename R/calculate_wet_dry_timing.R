#' Calculate Spatial Timing for Wet and Dry Seasons
#'
#' @description
#' Calculates the exact Day of the Year (DOY) for the onset and retreat of both 
#' the Wet Season and the Dry Season using a rolling 10-day rainfall sum.
#'
#' @details
#' This function evaluates the occurrence of distinct seasons using a rolling sum 
#' approach to filter out erratic "false start" rainfall events.
#' 
#' Mathematically, the function computes a rolling sum \eqn{S_t} of daily precipitation \eqn{P} 
#' over a window of \eqn{w} days (default \eqn{w = 10}):
#' \deqn{S_t = \sum_{i = t - w + 1}^{t} P_i}
#' 
#' The Day of the Year (DOY) for each seasonal event is defined by finding the earliest 
#' day \eqn{t} that satisfies specific threshold conditions, adjusted back to the start 
#' of the \eqn{w}-day window:
#' 
#' \strong{Wet Season Onset:}
#' \deqn{Onset_{wet} = \min(t) - (w - 1) \quad \text{where } S_t \ge W_{th} \text{ and } t > d_{spring}}
#' 
#' \strong{Wet Season Retreat:}
#' \deqn{Retreat_{wet} = \min(t) - (w - 1) \quad \text{where } S_t < D_{th} \text{ and } t > d_{autumn}}
#' 
#' \strong{Dry Season Onset & Retreat:}
#' Uses symmetrical logic, isolating periods where \eqn{S_t < D_{th}} during winter months 
#' and \eqn{S_t \ge W_{th}} during the spring transition.
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01").
#' @param region_name (Optional) Character. Admin region name (e.g., "Amhara") to filter.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH").
#' @param wet_threshold Numeric. The 10-day mm sum required to trigger the Wet Season. Default is 20.
#' @param dry_threshold Numeric. The 10-day mm sum required to trigger the Dry Season. Default is 5.
#' @param window_size Numeric. The number of days in the rolling window. Default is 10.
#' @param save_path (Optional) Character. File path to save the output raster.
#' @param plot_result Logical. If TRUE, plots the four spatial maps for the first year.
#'
#' @return A SpatRaster stack containing four layers per year: 
#'   `Wet_Onset`, `Wet_Retreat`, `Dry_Onset`, `Dry_Retreat`.
#' @export
#' @importFrom terra rast vect crop mask time nlyr app subset project crs writeRaster
calculate_wet_dry_timing <- function(nc_input, 
                                     start_date = NULL, 
                                     region_name = NULL,
                                     country = NULL,
                                     wet_threshold = 20,
                                     dry_threshold = 5,
                                     window_size = 10,
                                     save_path = NULL,
                                     plot_result = TRUE) {
  
  # --- 1. Load Data ---
  if (is.character(nc_input)) {
    if (!file.exists(nc_input)) stop("NetCDF file not found.")
    r <- terra::rast(nc_input)
  } else {
    r <- nc_input
  }
  
  # --- 2. Handle Boundary ---
  region_vect <- NULL
  if (!is.null(country)) {
    level_needed <- ifelse(is.null(region_name), 0, 1) 
    region_vect <- geodata::gadm(country, level = level_needed, path = tempdir())
    if (!is.null(region_name)) {
      region_vect <- region_vect[region_vect$NAME_1 == region_name, ]
    }
    r <- terra::crop(r, terra::project(region_vect, terra::crs(r)), mask = TRUE)
  }
  
  # --- 3. Date Handling ---
  raw_dates <- terra::time(r)
  if ((is.null(raw_dates) || all(is.na(raw_dates))) && !is.null(start_date)) {
    terra::time(r) <- seq(from = as.Date(start_date), by = "day", length.out = terra::nlyr(r))
    dates <- terra::time(r)
  } else {
    dates <- as.Date(raw_dates)
  }
  
  years <- as.numeric(format(dates, "%Y"))
  unique_years <- unique(years)
  
  # --- 4. The Core Math: 10-Day Rolling Sums ---
  core_timing <- function(x) {
    if (all(is.na(x))) return(c(Wet_Onset = NA, Wet_Retreat = NA, Dry_Onset = NA, Dry_Retreat = NA))
    x[is.na(x)] <- 0 
    
    # Calculate the rolling sum of rainfall
    roll_sum <- as.numeric(stats::filter(x, rep(1, window_size), sides = 1))
    days_seq <- seq_along(x)
    
    # -- WET SEASON LOGIC --
    wet_onset_cands <- which(roll_sum >= wet_threshold & days_seq > 90)
    wet_onset <- if(length(wet_onset_cands) > 0) min(wet_onset_cands) - (window_size - 1) else NA
    
    wet_retreat_cands <- which(roll_sum < dry_threshold & days_seq > 240)
    wet_retreat <- if(length(wet_retreat_cands) > 0) min(wet_retreat_cands) - (window_size - 1) else NA
    
    # -- DRY SEASON LOGIC --
    dry_onset_cands <- which(roll_sum < dry_threshold & days_seq > 270)
    dry_onset <- if(length(dry_onset_cands) > 0) min(dry_onset_cands) - (window_size - 1) else NA
    
    dry_retreat_cands <- which(roll_sum >= wet_threshold & days_seq < 90)
    dry_retreat <- if(length(dry_retreat_cands) > 0) min(dry_retreat_cands) - (window_size - 1) else NA
    
    return(c(Wet_Onset = wet_onset, Wet_Retreat = wet_retreat, Dry_Onset = dry_onset, Dry_Retreat = dry_retreat))
  }
  
  # --- 5. Process Year by Year ---
  message("Calculating spatial Julian dates for Wet and Dry seasons...")
  timing_layers <- list()
  
  for (yr in unique_years) {
    idx <- which(years == yr)
    if (length(idx) == 0) next 
    
    r_year <- terra::subset(r, idx)
    year_timing <- terra::app(r_year, fun = core_timing)
    names(year_timing) <- c(paste0("Wet_Onset_", yr), paste0("Wet_Retreat_", yr),
                            paste0("Dry_Onset_", yr), paste0("Dry_Retreat_", yr))
    
    timing_layers[[length(timing_layers) + 1]] <- year_timing
  }
  
  # --- 6. Stack and Format ---
  timing_stack <- terra::rast(timing_layers)
  if (!is.null(save_path)) terra::writeRaster(timing_stack, filename = save_path, overwrite = TRUE)
  
  if (plot_result) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(2, 2))
    
    pal_wet <- grDevices::hcl.colors(50, "YlGnBu", rev = TRUE)
    pal_dry <- grDevices::hcl.colors(50, "YlOrRd", rev = TRUE)
    
    terra::plot(timing_stack[[1]], main = "Wet Season Onset (DOY)", col = pal_wet)
    terra::plot(timing_stack[[2]], main = "Wet Season Retreat (DOY)", col = pal_wet)
    terra::plot(timing_stack[[3]], main = "Dry Season Onset (DOY)", col = pal_dry)
    terra::plot(timing_stack[[4]], main = "Dry Season Retreat (DOY)", col = pal_dry)
  }
  
  message("✅ Wet and Dry Season timing calculation complete!")
  return(timing_stack)
}
#' Calculate Spatial Wet-Season Dry Spells (WSDS)
#'
#' @description
#' Ingests a daily precipitation NetCDF file, filters the data to a specific 
#' wet season, and spatially calculates the maximum number of consecutive 
#' dry days for every pixel across the map.
#'
#' @details
#' Unlike standard CDD which evaluates an entire calendar year, WSDS isolates 
#' mid-season agricultural droughts. 
#' 
#' Mathematically, WSDS evaluates an indicator sequence \eqn{I_t} for daily precipitation \eqn{P_t}:
#' \deqn{I_t = 1 \text{ if } P_t < \text{threshold}}
#' \deqn{I_t = 0 \text{ if } P_t \ge \text{threshold}}
#' 
#' The index calculates the maximum length of consecutive days \eqn{[a, b]} within 
#' the wet season where \eqn{I_t = 1}:
#' \deqn{WSDS = \max \left( \sum_{t=a}^{b} I_t \right)}
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01") if the NetCDF lacks built-in time.
#' @param region (Optional) Shapefile path or SpatVector to crop the map.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Amhara") to filter.
#' @param target_season Numeric vector. The specific months of your wet season (e.g., 6:9 for Kiremt).
#' @param season_name Character. Name of the season for labeling the output (e.g., "Kiremt").
#' @param threshold Numeric. Daily precipitation threshold to define a "dry" day. Default is 1.0 mm.
#' @param save_path (Optional) Character. File path to save the output raster (e.g., "output/Kiremt_WSDS_1981.tif"). 
#'   Supports .tif or .nc extensions.
#' @param plot_result Logical. If TRUE, plots the spatial map for the first year.
#'
#' @return A SpatRaster stack where each layer represents the WSDS for a specific year.
#' @export
#' @importFrom terra rast vect crop mask time nlyr app subset project crs writeRaster
calculate_spatial_wsds <- function(nc_input, 
                                   start_date = NULL, 
                                   region = NULL,
                                   country = NULL,
                                   region_name = NULL,
                                   target_season = 6:9,  
                                   season_name = "Kiremt",
                                   threshold = 1.0,
                                   save_path = NULL,
                                   plot_result = TRUE) {
  
  # --- 1. Load Data ---
  if (is.character(nc_input)) {
    if (!file.exists(nc_input)) stop("NetCDF file not found.")
    r <- terra::rast(nc_input)
  } else {
    r <- nc_input
  }
  
  # --- 2. Handle Auto-Boundary & Cropping ---
  region_vect <- NULL
  if (!is.null(region)) {
    region_vect <- if (is.character(region)) terra::vect(region) else region
  } else if (!is.null(country)) {
    if (!requireNamespace("geodata", quietly = TRUE)) {
      stop("The 'geodata' package is required. Please install it using install.packages('geodata').")
    }
    level_needed <- ifelse(is.null(region_name), 0, 1) 
    message(paste("Downloading boundary data for ISO3:", country, "..."))
    region_vect <- tryCatch(
      geodata::gadm(country, level = level_needed, path = tempdir()), 
      error = function(e) stop("Boundary download failed. Check country code or internet connection.")
    )
    if (!is.null(region_name)) {
      if (!(region_name %in% region_vect$NAME_1)) {
        stop(paste("Region not found. Available regions are:", paste(unique(region_vect$NAME_1), collapse=", ")))
      }
      region_vect <- region_vect[region_vect$NAME_1 == region_name, ]
    }
  }
  
  if (!is.null(region_vect)) {
    if (terra::crs(region_vect) != terra::crs(r)) {
      region_vect <- terra::project(region_vect, terra::crs(r))
    }
    message("Clipping NetCDF data to the boundary...")
    r <- terra::crop(r, region_vect, mask = TRUE)
  }
  
  # --- 3. Date Handling ---
  raw_dates <- terra::time(r)
  if ((is.null(raw_dates) || all(is.na(raw_dates))) && !is.null(start_date)) {
    terra::time(r) <- seq(from = as.Date(start_date), by = "day", length.out = terra::nlyr(r))
    dates <- terra::time(r)
  } else if (!is.null(raw_dates) && !all(is.na(raw_dates))) {
    dates <- as.Date(raw_dates)
  } else {
    stop("Could not read dates. Please provide 'start_date'.")
  }
  
  years <- as.numeric(format(dates, "%Y"))
  months <- as.numeric(format(dates, "%m"))
  unique_years <- unique(years)
  
  # --- 4. Define the Core Math for terra::app ---
  core_wsds <- function(x) {
    x_valid <- stats::na.omit(x)
    if (length(x_valid) == 0) return(NA)
    
    streaks <- base::rle(x_valid < threshold)
    if (any(streaks$values)) {
      return(max(streaks$lengths[streaks$values == TRUE]))
    } else {
      return(0)
    }
  }
  
  # --- 5. Process Spatial WSDS Year by Year ---
  message(paste("Calculating Spatial WSDS for season:", season_name))
  wsds_layers <- list()
  valid_years <- c()
  
  for (yr in unique_years) {
    idx <- which(years == yr & months %in% target_season)
    if (length(idx) == 0) next 
    
    r_season <- terra::subset(r, idx)
    year_wsds <- terra::app(r_season, fun = core_wsds)
    
    wsds_layers[[length(wsds_layers) + 1]] <- year_wsds
    valid_years <- c(valid_years, yr)
  }
  
  # --- 6. Stack, Save, and Format Output ---
  if (length(wsds_layers) > 0) {
    wsds_stack <- terra::rast(wsds_layers)
    names(wsds_stack) <- paste0("WSDS_", season_name, "_", valid_years)
    
    # Save the raster to disk if a path is provided
    if (!is.null(save_path)) {
      message(paste("Saving output to:", save_path))
      terra::writeRaster(wsds_stack, filename = save_path, overwrite = TRUE)
    }
    
    # --- 7. Plotting ---
    if (plot_result) {
      drought_pal <- grDevices::hcl.colors(50, "YlOrRd", rev = TRUE)
      map_to_plot <- wsds_stack[[1]]
      
      terra::plot(map_to_plot, 
                  main = paste("Longest Dry Spell (Days) -", season_name, valid_years[1]),
                  col = drought_pal)
      
      if (!is.null(region_vect)) {
        terra::plot(region_vect, add = TRUE, lwd = 2, border = "black")
      }
    }
    
    message("✅ WSDS spatial calculation complete!")
    return(wsds_stack)
  } else {
    stop("No data found for the specified target season.")
  }
}
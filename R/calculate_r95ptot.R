#' Calculate R95pTOT (Contribution of Very Wet Days to Total Rainfall)
#'
#' Calculates the percentage (%) of total annual precipitation that comes from 
#' days exceeding the 95th percentile.
#' 
#' @details
#' The R95pTOT index is calculated using the following formula:
#' \deqn{R95pTOT = \left( \frac{R95p}{PRCPTOT} \right) \times 100}
#' 
#' Where:
#' \itemize{
#'   \item \strong{R95p}: Annual sum of precipitation from days > 95th percentile.
#'   \item \strong{PRCPTOT}: Annual total precipitation from wet days (> 1mm).
#' }
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01").
#' @param baseline_range Numeric Vector. Baseline years (e.g., c(1981, 2010)).
#' @param region (Optional) Shapefile path or SpatVector to mask data.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia") to filter.
#' @param wet_threshold Numeric. Minimum rain (mm) for wet days. Default 1.
#' @param plot_result Logical. If TRUE, plots the results.
#' @param save_result Logical. If TRUE, saves the output as a .tif file.
#' @param output_dir Character. Directory to save the file.
#' @param output_name Character. Prefix for the saved file name.
#'
#' @return A SpatRaster object containing the R95pTOT values (%) for each year.
#' @export
#' @importFrom terra rast vect crop mask time nlyr app time<- subset quantile project crs writeRaster
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot
calculate_r95ptot <- function(nc_input, start_date = NULL, 
                              baseline_range = c(1981, 2010),
                              region = NULL,
                              country = NULL,       
                              region_name = NULL,   
                              wet_threshold = 1, 
                              plot_result = TRUE,
                              save_result = FALSE,              # <--- NEW
                              output_dir = "Output_Maps",       # <--- NEW
                              output_name = "R95pTOT_Index") {  # <--- NEW
  
  # --- 1. Load Data ---
  if (is.character(nc_input)) {
    if (!file.exists(nc_input)) stop("NetCDF file not found.")
    r <- terra::rast(nc_input)
  } else {
    r <- nc_input
  }
  
  # --- 2. Handle Region & Country (THE NEW LOGIC) ---
  region_vect <- NULL
  
  # Case A: User provided a custom shapefile path or object
  if (!is.null(region)) {
    if (is.character(region)) {
      region_vect <- terra::vect(region)
    } else {
      region_vect <- region
    }
  } 
  # Case B: User provided Country Code (Auto-Download)
  else if (!is.null(country)) {
    # Check if geodata is installed
    if (!requireNamespace("geodata", quietly = TRUE)) {
      stop("Package 'geodata' is needed for this feature. Please run: install.packages('geodata')")
    }
    
    message(paste("Downloading boundary for:", country, "..."))
    
    # Download Level 1 (Regions) if region_name is asked, otherwise Level 0 (Whole Country)
    level_needed <- ifelse(is.null(region_name), 0, 1)
    
    # Download GADM data (saved to tempdir so it doesn't clutter user's PC)
    region_vect <- tryCatch({
      geodata::gadm(country = country, level = level_needed, path = tempdir())
    }, error = function(e) {
      stop("Could not download country data. Check internet connection or ISO3 code.")
    })
    
    # Filter for specific region name if requested
    if (!is.null(region_name)) {
      # Check if the requested region exists in the downloaded map
      if (!(region_name %in% region_vect$NAME_1)) {
        stop(paste0("Region '", region_name, "' not found in ", country, ".\n",
                    "Available regions: ", paste(head(unique(region_vect$NAME_1)), collapse=", ")))
      }
      message(paste("Selecting Region:", region_name))
      region_vect <- region_vect[region_vect$NAME_1 == region_name, ]
    }
  }
  
  # --- Perform Clipping if a region was found ---
  if (!is.null(region_vect)) {
    # Match Coordinate Systems
    if (terra::crs(region_vect) != terra::crs(r)) {
      region_vect <- terra::project(region_vect, terra::crs(r))
    }
    message("Clipping data to region...")
    r <- terra::crop(r, region_vect, mask = TRUE)
  }
  
  # --- 3. Date Handling ---
  raw_dates <- terra::time(r)
  if ((is.null(raw_dates) || all(is.na(raw_dates))) && !is.null(start_date)) {
    message("Generating daily sequence from: ", start_date)
    terra::time(r) <- seq(from = as.Date(start_date), by = "day", length.out = terra::nlyr(r))
    dates <- terra::time(r)
  } else if (!is.null(raw_dates) && !all(is.na(raw_dates))) {
    dates <- as.Date(raw_dates)
  } else {
    stop("Could not read dates. Please provide 'start_date'.")
  }
  
  # --- 4. STEP 1: Calculate Baseline Threshold ---
  message(paste0("Calculating 95th percentile baseline (", baseline_range[1], "-", baseline_range[2], ")..."))
  
  years <- as.numeric(format(dates, "%Y"))
  baseline_idx <- which(years >= baseline_range[1] & years <= baseline_range[2])
  
  if (length(baseline_idx) == 0) {
    warning("Baseline years not found. Using all data as baseline.")
    baseline_idx <- 1:terra::nlyr(r)
  }
  
  baseline_stack <- terra::subset(r, baseline_idx)
  
  calc_quantile <- function(x) {
    if (all(is.na(x))) return(NA)
    wet_values <- x[x >= wet_threshold]
    if (length(wet_values) == 0) return(0)
    return(stats::quantile(wet_values, probs = 0.95, na.rm = TRUE))
  }
  
  # 95th Percentile Map
  percentile_map <- terra::app(baseline_stack, fun = calc_quantile)
  
  # --- 5. STEP 2: Annual R95pTOT Calculation ---
  message("Calculating Annual R95pTOT (%) ...")
  
  unique_years <- unique(years)
  annual_layers <- list()
  valid_years <- c()
  
  for (yr in unique_years) {
    yr_idx <- which(years == yr)
    if (length(yr_idx) == 0) next
    
    yr_stack <- terra::subset(r, yr_idx)
    
    # A. Calculate Total Rain (PRCPTOT)
    # Mask out values < 1mm (optional, but standard for wet-day totals)
    wet_mask <- yr_stack >= wet_threshold
    prcptot <- sum(yr_stack * wet_mask, na.rm = TRUE)
    
    # B. Calculate Extreme Rain (R95p)
    extreme_mask <- yr_stack > percentile_map
    r95p <- sum(yr_stack * extreme_mask, na.rm = TRUE)
    
    # C. Calculate Percentage: (R95p / PRCPTOT) * 100
    # Handle division by zero (if total rain is 0, result is 0)
    r95ptot <- (r95p / prcptot) * 100
    r95ptot[prcptot == 0] <- 0 
    
    annual_layers[[length(annual_layers) + 1]] <- r95ptot
    valid_years <- c(valid_years, yr)
  }
  
  result_stack <- terra::rast(annual_layers)
  names(result_stack) <- paste0("R95pTOT_", valid_years)
  
  # --- 6. Plotting ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Use "Spectral" or "YlGnBu" to show percentage intensity
    # Blue/Purple = High contribution of extremes
    pct_pal <- hcl.colors(50, "Spectral", rev = TRUE)
    
    par(mfrow = c(1, 1), mar = c(3, 3, 3, 4))
    
    terra::plot(result_stack[[1]], 
                main = paste("R95pTOT (Extreme Rain Contribution) - Year", valid_years[1]), 
                col = pct_pal,
                plg = list(title = "Percentage (%)", title.cex = 0.8))
    
    # Add boundary if it exists
    if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "black")
  }
  
  # --- 7. Save Results (NEW BLOCK) ---
  if (save_result) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    if (inherits(result_stack, "SpatRaster")) {
      out_path <- file.path(output_dir, paste0(output_name, ".tif"))
      terra::writeRaster(result_stack, out_path, overwrite = TRUE)
      message("Result saved to: ", out_path)
    } else if (is.list(result_stack)) {
      for (item_name in names(result_stack)) {
        out_path <- file.path(output_dir, paste0(output_name, "_", item_name, ".tif"))
        terra::writeRaster(result_stack[[item_name]], out_path, overwrite = TRUE)
        message("Result saved to: ", out_path)
      }
    }
  }
  
  return(result_stack)
}
#' Calculate R99pTOT (Contribution of Extremely Wet Days)
#'
#' Calculates the percentage (%) of total annual precipitation that comes from 
#' days exceeding the 99th percentile (Extremely Wet Days).
#'
#' @details
#' The R99pTOT index represents the contribution of extreme events to the total water budget:
#' \deqn{R99pTOT = \left( \frac{R99p}{PRCPTOT} \right) \times 100}
#' 
#' Where:
#' \itemize{
#'   \item \strong{R99p}: Annual sum of precipitation from days > 99th percentile.
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
#' @return A SpatRaster object containing the R99pTOT values (%) for each year.
#' @export
#' @importFrom terra rast vect crop mask time nlyr app time<- subset quantile global project crs writeRaster
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot
calculate_r99ptot <- function(nc_input, start_date = NULL, 
                              baseline_range = c(1981, 2010),
                              region = NULL,
                              country = NULL,       
                              region_name = NULL,   
                              wet_threshold = 1, 
                              plot_result = TRUE,
                              save_result = FALSE,             
                              output_dir = "Output_Maps",       
                              output_name = "R99pTOT_Index") {  
  
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
  
  # --- 4. STEP 1: Calculate 99th Percentile Baseline ---
  message(paste0("Calculating 99th percentile baseline (", baseline_range[1], "-", baseline_range[2], ")..."))
  
  years <- as.numeric(format(dates, "%Y"))
  baseline_idx <- which(years >= baseline_range[1] & years <= baseline_range[2])
  
  if (length(baseline_idx) == 0) {
    warning("Baseline years not found. Using all data as baseline.")
    baseline_idx <- 1:terra::nlyr(r)
  }
  
  baseline_stack <- terra::subset(r, baseline_idx)
  
  calc_quantile_99 <- function(x) {
    if (all(is.na(x))) return(NA)
    wet_values <- x[x >= wet_threshold]
    if (length(wet_values) == 0) return(0)
    return(stats::quantile(wet_values, probs = 0.99, na.rm = TRUE))
  }
  
  percentile_map <- terra::app(baseline_stack, fun = calc_quantile_99)
  
  # --- 5. STEP 2: Annual R99pTOT Calculation ---
  message("Calculating Annual R99pTOT (%) ...")
  
  unique_years <- unique(years)
  annual_layers <- list()
  valid_years <- c()
  
  for (yr in unique_years) {
    yr_idx <- which(years == yr)
    if (length(yr_idx) == 0) next
    
    yr_stack <- terra::subset(r, yr_idx)
    
    # A. Total Rain (PRCPTOT)
    wet_mask <- yr_stack >= wet_threshold
    prcptot <- sum(yr_stack * wet_mask, na.rm = TRUE)
    
    # B. Extreme Rain (R99p)
    extreme_mask <- yr_stack > percentile_map
    r99p <- sum(yr_stack * extreme_mask, na.rm = TRUE)
    
    # C. Percentage
    r99ptot <- (r99p / prcptot) * 100
    
    # Fix NaNs (0/0) or Infs (X/0)
    r99ptot[prcptot == 0] <- 0 
    r99ptot[is.na(r99ptot)] <- 0
    
    annual_layers[[length(annual_layers) + 1]] <- r99ptot
    valid_years <- c(valid_years, yr)
  }
  
  result_stack <- terra::rast(annual_layers)
  names(result_stack) <- paste0("R99pTOT_", valid_years)
  
  # --- 6. Plotting (Wrapped in TryCatch) ---
  if (plot_result) {
    tryCatch({
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par))
      
      # Use "Spectral" (Red-Yellow-Blue) for Percentage intensity
      pct_pal <- hcl.colors(50, "Spectral", rev = TRUE)
      
      # Check data range safely
      val_range <- terra::global(result_stack[[1]], "range", na.rm = TRUE)
      min_val <- val_range[1, 1]
      max_val <- val_range[1, 2]
      
      par(mfrow = c(1, 1), mar = c(3, 3, 3, 4))
      
      # If flat map (e.g., 0% contribution everywhere)
      if (isTRUE(min_val == max_val)) {
        message(paste0("Map contains constant value (", min_val, "%). Plotting without legend."))
        terra::plot(result_stack[[1]], 
                    main = paste("Annual R99pTOT - Year", valid_years[1]),
                    col = "lightgrey", 
                    legend = FALSE)
      } else {
        terra::plot(result_stack[[1]], 
                    main = paste("R99pTOT (Extreme Contribution) - Year", valid_years[1]), 
                    col = pct_pal,
                    plg = list(title = "Percentage (%)", title.cex = 0.8))
      }
      
      # Add boundary if it exists
      if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "black")
      
    }, error = function(e) {
      warning(paste("Plotting failed, but calculations are correct. Error:", e$message))
    })
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
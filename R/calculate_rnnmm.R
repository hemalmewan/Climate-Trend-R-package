#' Calculate Rnnmm (Number of Days > Threshold)
#'
#' Counts the number of days where precipitation exceeds a specific threshold.
#' 
#' @details
#' The Rnnmm index (e.g., R10mm, R20mm) is defined as:
#' \deqn{Rnnmm = \sum_{i=1}^{N} I(RR_{i} > \text{threshold})}
#' where \eqn{I} is an indicator function (1 if true, 0 if false).
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01") if file has no dates.
#' @param region (Optional) Shapefile path or SpatVector to mask the data.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia") to filter.
#' @param threshold Numeric. The precipitation threshold in mm. 
#'    Default is 10 (calculates R10mm). Change to 20 for R20mm.
#' @param timescale Character vector. "monthly", "annual", or both.
#' @param plot_result Logical. If TRUE, plots the results.
#' @param save_result Logical. If TRUE, saves the output as a .tif file.
#' @param output_dir Character. Directory to save the file.
#' @param output_name Character. Prefix for the saved file name.
#'
#' @return A named list containing 'monthly' (SpatRaster) and/or 'annual' (SpatRaster).
#' @export
#' @importFrom terra rast vect crop mask time nlyr time<- app project crs writeRaster
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot mtext
calculate_rnnmm <- function(nc_input, start_date = NULL, 
                            region = NULL,
                            country = NULL,       
                            region_name = NULL,   
                            threshold = 10, 
                            timescale = c("monthly", "annual"), 
                            plot_result = TRUE,
                            save_result = FALSE,              # <--- NEW
                            output_dir = "Output_Maps",       # <--- NEW
                            output_name = "Rnnmm_Index") {    # <--- NEW
  
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
  
  # --- 4. Calculation ---
  results <- list()
  
  # A. Monthly Calculation
  if ("monthly" %in% timescale) {
    message(paste0("Calculating Monthly R", threshold, "mm..."))
    monthly_layers <- list()
    valid_months   <- c() 
    
    for (m in 1:12) {
      idx <- which(format(dates, "%m") == sprintf("%02d", m))
      if (length(idx) == 0) next 
      
      month_stack <- terra::subset(r, idx)
      
      # Logic: Count days > threshold
      count_days <- sum(month_stack > threshold, na.rm = TRUE)
      
      monthly_layers[[length(monthly_layers) + 1]] <- count_days
      valid_months <- c(valid_months, month.abb[m])
    }
    
    if (length(monthly_layers) > 0) {
      monthly_stack <- terra::rast(monthly_layers)
      names(monthly_stack) <- valid_months
      results$monthly <- monthly_stack
    }
  }
  
  # B. Annual Calculation
  if ("annual" %in% timescale) {
    message(paste0("Calculating Annual R", threshold, "mm..."))
    # Summing across the entire year for > threshold
    results$annual <- sum(r > threshold, na.rm = TRUE)
    names(results$annual) <- paste0("Annual_R", threshold, "mm")
  }
  
  # --- 5. Plotting (Viridis Style) ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Use "Viridis" palette (Standard for count data)
    count_pal <- hcl.colors(50, "Viridis", rev = FALSE) 
    
    # Title prefix based on threshold
    idx_name <- paste0("R", threshold, "mm")
    
    if ("monthly" %in% timescale && !is.null(results$monthly)) {
      n <- length(valid_months)
      par(mfrow = c(ceiling(n/3), 3), mar = c(2, 2, 2, 1))
      for (i in 1:n) {
        terra::plot(results$monthly[[i]], 
                    main = paste(idx_name, "-", valid_months[i]), 
                    col = count_pal)
        # Add boundary if it exists
        if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "white")
      }
    }
    
    if ("annual" %in% timescale) {
      if ("monthly" %in% timescale && dev.interactive()) dev.new()
      par(mfrow = c(1, 1), mar = c(3, 3, 3, 4))
      terra::plot(results$annual, 
                  main = paste("Annual Frequency:", idx_name, "(Days >", threshold, "mm)"), 
                  col = count_pal)
      
      # Add boundary if it exists
      if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "white")
    }
  }
  
  # --- 6. Save Results (NEW BLOCK) ---
  if (save_result) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    if (is.list(results)) {
      for (item_name in names(results)) {
        out_path <- file.path(output_dir, paste0(output_name, "_", item_name, ".tif"))
        terra::writeRaster(results[[item_name]], out_path, overwrite = TRUE)
        message("Result saved to: ", out_path)
      }
    } else if (inherits(results, "SpatRaster")) {
      out_path <- file.path(output_dir, paste0(output_name, ".tif"))
      terra::writeRaster(results, out_path, overwrite = TRUE)
      message("Result saved to: ", out_path)
    }
  }
  
  return(results)
}
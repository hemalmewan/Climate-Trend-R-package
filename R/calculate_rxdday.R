#' Calculate Rxdday (Maximum x-day Precipitation)
#'
#' Calculates the maximum precipitation amount accumulated over a moving window of x days.
#' Common uses are Rx1day (window=1) and Rx5day (window=5).
#'
#' @details
#' The Rxdday index represents the highest rainfall accumulation over a specified duration:
#' \deqn{Rx\text{day} = \max(\text{rolling\_sum}(\text{Precip}, \text{window}))}
#' 
#' Where:
#' \itemize{
#'   \item \strong{window}: The number of consecutive days (e.g., 1 or 5).
#'   \item \strong{rolling_sum}: The sum of precipitation for each overlapping window.
#' }
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01").
#' @param window_size Integer. The number of days for the rolling sum. Default 1 (Rx1day). Use 5 for Rx5day.
#' @param timescale Character vector. "monthly", "annual", or both.
#' @param region (Optional) Shapefile path or SpatVector to mask data.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia") to filter.
#' @param plot_result Logical. If TRUE, plots the results.
#' @param save_result Logical. If TRUE, saves the output as a .tif file.
#' @param output_dir Character. Directory to save the file.
#' @param output_name Character. Prefix for the saved file name.
#'
#' @return A named list containing 'monthly' and/or 'annual' SpatRasters.
#' @export
#' @importFrom terra rast vect crop mask time nlyr app time<- roll project crs writeRaster
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot
calculate_rxdday <- function(nc_input, start_date = NULL, 
                             window_size = 1,
                             timescale = c("monthly", "annual"),
                             region = NULL,
                             country = NULL,       
                             region_name = NULL,   
                             plot_result = TRUE,
                             save_result = FALSE,              # <--- NEW
                             output_dir = "Output_Maps",       # <--- NEW
                             output_name = "Rxdday_Index") {   # <--- NEW
  
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
  
  # --- 4. Rolling Sum Calculation (Crucial Step) ---
  if (window_size > 1) {
    message(paste("Calculating rolling sum for", window_size, "days..."))
    # terra::roll calculates moving sum across time (3rd dimension)
    # type="to" means the value at day X includes day X and previous (window-1) days
    r_rolled <- terra::roll(r, n = window_size, fun = sum, type = "to", na.rm = TRUE)
  } else {
    r_rolled <- r # Rx1day is just the raw daily data
  }
  
  # --- 5. Calculation Loop ---
  results <- list()
  
  # A. Monthly Calculation
  if ("monthly" %in% timescale) {
    message(paste0("Calculating Monthly Rx", window_size, "day..."))
    monthly_layers <- list()
    valid_months <- c()
    
    for (m in 1:12) {
      idx <- which(format(dates, "%m") == sprintf("%02d", m))
      if (length(idx) == 0) next 
      
      month_stack <- terra::subset(r_rolled, idx)
      
      # Max value in that month
      max_val <- max(month_stack, na.rm = TRUE)
      
      monthly_layers[[length(monthly_layers) + 1]] <- max_val
      valid_months <- c(valid_months, month.abb[m])
    }
    
    monthly_stack <- terra::rast(monthly_layers)
    names(monthly_stack) <- valid_months
    results$monthly <- monthly_stack
  }
  
  # B. Annual Calculation
  if ("annual" %in% timescale) {
    message(paste0("Calculating Annual Rx", window_size, "day..."))
    # Max value in the entire year
    results$annual <- max(r_rolled, na.rm = TRUE)
    names(results$annual) <- paste0("Annual_Rx", window_size, "day")
  }
  
  # --- 6. Plotting ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Use "Lajolla" or "Blues" for Intensity
    int_pal <- hcl.colors(50, "Lajolla", rev = TRUE)
    
    idx_name <- paste0("Rx", window_size, "day")
    
    if ("monthly" %in% timescale && !is.null(results$monthly)) {
      n <- length(valid_months)
      par(mfrow = c(ceiling(n/3), 3), mar = c(2, 2, 2, 1))
      for (i in 1:n) {
        terra::plot(results$monthly[[i]], 
                    main = paste(idx_name, "-", valid_months[i]), 
                    col = int_pal)
        
        # Add boundary if it exists
        if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "black")
      }
    }
    
    if ("annual" %in% timescale) {
      if ("monthly" %in% timescale && dev.interactive()) dev.new()
      par(mfrow = c(1, 1), mar = c(3, 3, 3, 4))
      terra::plot(results$annual, 
                  main = paste("Annual Maximum", idx_name, "(mm)"), 
                  col = int_pal,
                  plg = list(title = "Intensity (mm)", title.cex = 0.8))
      
      # Add boundary if it exists
      if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "black")
    }
  }
  
  # --- 7. Save Results (NEW BLOCK) ---
  if (save_result) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Check if the result is a List of SpatRasters
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
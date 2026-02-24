#' Calculate PRCPTOT (Total Wet-Day Precipitation)
#'
#' Calculates the total annual precipitation from wet days.
#' 
#' @details
#' The PRCPTOT index is defined as:
#' \deqn{PRCPTOT = \sum_{i=1}^{N} RR_{i} \quad \text{where } RR_{i} \ge 1mm}
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01").
#' @param region (Optional) Shapefile path or SpatVector.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia") to filter.
#' @param threshold Numeric. Minimum rain (mm). Default 1.
#' @param timescale Character. "monthly" or "annual".
#' @param plot_result Logical. If TRUE, plots the results.
#' @param save_result Logical. If TRUE, saves the output as a .tif file.
#' @param output_dir Character. Directory to save the file.
#' @param output_name Character. Prefix for the saved file name.
#'
#' @return A named list containing 'monthly' and/or 'annual' SpatRasters.
#' @export
#' @importFrom terra rast vect crop mask time nlyr time<- project crs writeRaster
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot
calculate_prcptot <- function(nc_input, start_date = NULL, 
                              region = NULL, 
                              country = NULL,       
                              region_name = NULL,   
                              threshold = 1, 
                              timescale = c("monthly", "annual"), 
                              plot_result = TRUE,
                              save_result = FALSE,              
                              output_dir = "Output_Maps",      
                              output_name = "PRCPTOT_Index") { 
  
  # --- 1. Load Data ---
  if (is.character(nc_input)) {
    if (!file.exists(nc_input)) stop("NetCDF file not found.")
    r <- terra::rast(nc_input)
  } else {
    r <- nc_input
  }
  
  # --- 2. Handle Region & Country ---
  region_vect <- NULL
  
  # Case A: User provided a custom shapefile path or object
  if (!is.null(region)) {
    if (is.character(region)) {
      region_vect <- terra::vect(region)
    } else {
      region_vect <- region
    }
  } 
  # Case B : User provided Country Code (Auto-Download)
  else if (!is.null(country)) {
    if (!requireNamespace("geodata", quietly = TRUE)) {
      stop("Package 'geodata' is needed for this feature. Please run: install.packages('geodata')")
    }
    
    message(paste("Downloading boundary for:", country, "..."))
    level_needed <- ifelse(is.null(region_name), 0, 1)
    
    region_vect <- tryCatch({
      geodata::gadm(country = country, level = level_needed, path = tempdir())
    }, error = function(e) {
      stop("Could not download country data. Check internet connection or ISO3 code.")
    })
    
    if (!is.null(region_name)) {
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
  monthly_layers <- list()
  valid_months   <- c() 
  
  for (m in 1:12) {
    idx <- which(format(dates, "%m") == sprintf("%02d", m))
    if (length(idx) == 0) next 
    
    month_stack <- terra::subset(r, idx)
    
    # PRCPTOT logic
    wet_mask <- month_stack >= threshold
    prcptot_month <- sum(month_stack * wet_mask, na.rm = TRUE)
    
    monthly_layers[[length(monthly_layers) + 1]] <- prcptot_month
    valid_months <- c(valid_months, month.abb[m])
  }
  
  if (length(monthly_layers) == 0) stop("No data found for any month.")
  
  monthly_stack <- terra::rast(monthly_layers)
  names(monthly_stack) <- valid_months
  
  # --- 5. Prepare Output ---
  results <- list()
  if ("monthly" %in% timescale) results$monthly <- monthly_stack
  if ("annual" %in% timescale) {
    results$annual <- sum(monthly_stack, na.rm = TRUE)
    names(results$annual) <- "Annual_PRCPTOT"
  }
  
  # --- 6. Plotting ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    blue_pal <- hcl.colors(50, "Blues 3", rev = FALSE)
    
    if ("annual" %in% timescale) {
      par(mfrow = c(1, 1))
      terra::plot(results$annual, main = "Annual PRCPTOT", col = blue_pal)
      
      if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "black")
    }
  }
  
  # --- 7. Save Results (NEW BLOCK) ---
  if (save_result) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Iterate through the list and save each component
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
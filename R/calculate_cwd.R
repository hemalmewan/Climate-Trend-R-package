#' Calculate CWD (Consecutive Wet Days)
#'
#' Calculates the maximum length of a wet spell (consecutive days >= threshold).
#' 
#' @details
#' The CWD index represents the maximum number of consecutive days with significant rainfall:
#' \deqn{CWD = \max(\text{count}(\text{days where } RR \ge 1mm))}
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01").
#' @param region (Optional) Shapefile path or SpatVector.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia") to filter.
#' @param threshold Numeric. Days with rain >= this value are considered "wet". Default 1.
#' @param plot_result Logical. If TRUE, plots the annual map. Default TRUE.
#' @param save_result Logical. If TRUE, saves the output as a .tif file.
#' @param output_dir Character. Directory to save the file.
#' @param output_name Character. Prefix for the saved file name.
#'
#' @return A named list containing 'annual' (SpatRaster).
#' @export
#' @importFrom terra rast vect crop mask time nlyr app time<- project crs writeRaster
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot
calculate_cwd <- function(nc_input, start_date = NULL, 
                          region = NULL,
                          country = NULL,       
                          region_name = NULL,   
                          threshold = 1, 
                          plot_result = TRUE,
                          save_result = FALSE,              # <--- NEW
                          output_dir = "Output_Maps",       # <--- NEW
                          output_name = "CWD_Index") {      # <--- NEW
  
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
    n_days <- terra::nlyr(r)
    terra::time(r) <- seq(from = as.Date(start_date), by = "day", length.out = n_days)
  } else if (is.null(raw_dates) && is.null(start_date)) {
    stop("Could not read dates. Please provide 'start_date'.")
  }
  
  # --- 4. Define Calculation Function ---
  calc_max_wet_spell <- function(x) {
    if (all(is.na(x))) return(NA)
    
    # Identify WET days (Greater than or Equal to threshold)
    is_wet <- x >= threshold
    is_wet[is.na(is_wet)] <- FALSE 
    
    if (!any(is_wet)) return(0) # No wet days
    
    # Run Length Encoding
    rle_res <- rle(is_wet)
    
    # Filter for streaks that are TRUE (wet)
    wet_streaks <- rle_res$lengths[rle_res$values == TRUE]
    
    if (length(wet_streaks) == 0) return(0)
    return(max(wet_streaks))
  }
  
  # --- 5. Run Calculation ---
  message("Calculating Annual CWD... (This processes pixel-by-pixel, please wait)")
  
  annual_cwd <- terra::app(r, fun = calc_max_wet_spell)
  names(annual_cwd) <- "Annual_CWD"
  
  results <- list(annual = annual_cwd)
  
  # --- 6. Plotting (Wet Theme) ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    wet_pal <- hcl.colors(50, "GnBu", rev = FALSE)
    
    par(mfrow = c(1, 1), mar = c(3, 3, 3, 4))
    terra::plot(results$annual, 
                main = paste("Annual Max Consecutive Wet Days (CWD) >=", threshold, "mm"), 
                col = wet_pal)
    
    if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "black")
  }
  
  # --- 7. Save Results (NEW BLOCK) ---
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
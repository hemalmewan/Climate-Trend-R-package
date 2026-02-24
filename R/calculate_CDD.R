#' Calculate CDD (Consecutive Dry Days)
#'
#' Calculates the maximum length of a dry spell (consecutive days < threshold).
#' 
#' @details
#' The CDD index represents the maximum number of consecutive days with low rainfall:
#' \deqn{CDD = \max(\text{count}(\text{days where } RR < 1mm))}
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01") if file has no dates.
#' @param region (Optional) Shapefile path or SpatVector.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia") to filter.
#' @param threshold Numeric. Days with rain < this value are considered "dry". Default 1.
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
calculate_cdd <- function(nc_input, start_date = NULL, 
                          region = NULL,
                          country = NULL,       
                          region_name = NULL,   
                          threshold = 1, 
                          plot_result = TRUE,
                          save_result = FALSE,              # <--- NEW
                          output_dir = "Output_Maps",       # <--- NEW
                          output_name = "CDD_Index") {      # <--- NEW
  
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
    n_days <- terra::nlyr(r)
    terra::time(r) <- seq(from = as.Date(start_date), by = "day", length.out = n_days)
  } else if (is.null(raw_dates) && is.null(start_date)) {
    stop("Could not read dates. Please provide 'start_date' to ensure data validity.")
  }
  
  # --- 4. Define Calculation Function ---
  calc_max_dry_spell <- function(x) {
    if (all(is.na(x))) return(NA)
    
    # Identify dry days
    is_dry <- x < threshold
    is_dry[is.na(is_dry)] <- FALSE 
    
    if (!any(is_dry)) return(0) # No dry days at all
    
    # Run Length Encoding to find streak lengths
    rle_res <- rle(is_dry)
    
    # Filter for streaks that are TRUE (dry)
    dry_streaks <- rle_res$lengths[rle_res$values == TRUE]
    
    if (length(dry_streaks) == 0) return(0)
    return(max(dry_streaks))
  }
  
  # --- 5. Run Calculation ---
  message("Calculating Annual CDD... (This processes pixel-by-pixel, please wait)")
  
  annual_cdd <- terra::app(r, fun = calc_max_dry_spell)
  names(annual_cdd) <- "Annual_CDD"
  
  results <- list(annual = annual_cdd)
  
  # --- 6. Plotting ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Use Drought colors: Yellow (short dry spells) -> Red (long dry spells)
    dry_pal <- hcl.colors(50, "YlOrRd", rev = TRUE)
    
    par(mfrow = c(1, 1), mar = c(3, 3, 3, 4))
    terra::plot(results$annual, 
                main = "Annual Maximum Consecutive Dry Days (CDD)", 
                col = dry_pal)
    
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
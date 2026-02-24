#' Internal Helper: Process Region Arguments
#'
#' Handles the logic for loading shapefiles OR downloading GADM boundaries.
#' This function is used internally by all climate indices.
#'
#' @param region SpatVector or character path.
#' @param country Character. ISO3 code (e.g., "ETH").
#' @param region_name Character. Admin region name (e.g., "Oromia").
#' @param ref_raster SpatRaster. Used to match Coordinate Reference Systems (CRS).
#'
#' @return A SpatVector object (projected to match ref_raster) or NULL.
#' @noRd
process_region_input <- function(region, country, region_name, ref_raster) {
  
  final_region <- NULL
  
  # --- Case A: User provided a Shapefile or Path ---
  if (!is.null(region)) {
    if (is.character(region)) {
      if (!file.exists(region)) stop("Region shapefile not found.")
      final_region <- terra::vect(region)
    } else {
      final_region <- region
    }
    
    # --- Case B: User provided Country Code (Auto-Download) ---
  } else if (!is.null(country)) {
    
    if (!requireNamespace("geodata", quietly = TRUE)) {
      stop("Package 'geodata' is required. Please run install.packages('geodata')")
    }
    
    message(paste("Downloading boundary for:", country, "..."))
    
    # Level 0 = Country Boundary, Level 1 = Province/State Boundaries
    level_needed <- ifelse(is.null(region_name), 0, 1)
    
    # Download GADM data (saved to tempdir so it doesn't clutter user's PC)
    shp <- tryCatch({
      geodata::gadm(country = country, level = level_needed, path = tempdir())
    }, error = function(e) {
      stop("Could not download country data. Check internet connection or ISO3 code.")
    })
    
    # Filter for specific region name if requested
    if (!is.null(region_name)) {
      # "NAME_1" is the standard GADM column for Level 1 names
      avail_names <- shp$NAME_1
      
      if (!(region_name %in% avail_names)) {
        # Intelligent Error Message: Show user the available options
        stop(paste0("Region '", region_name, "' not found in ", country, ".\n",
                    "Did you mean one of these? ", 
                    paste(head(avail_names, 5), collapse=", "), "..."))
      }
      
      message(paste("Selecting Region:", region_name))
      final_region <- shp[shp$NAME_1 == region_name, ]
    } else {
      final_region <- shp
    }
  }
  
  # --- Case C: Reproject if needed ---
  if (!is.null(final_region)) {
    if (terra::crs(final_region) != terra::crs(ref_raster)) {
      # message("Reprojecting region to match NetCDF coordinates...")
      final_region <- terra::project(final_region, terra::crs(ref_raster))
    }
  }
  
  return(final_region)
}
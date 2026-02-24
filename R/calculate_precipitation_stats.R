#' Calculate Annual and Seasonal Precipitation Totals
#'
#' Calculates the total accumulated precipitation for the entire year 
#' and/or for specific seasons defined by the user.
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param start_date Character/Date. Start date (e.g. "1981-01-01").
#' @param region (Optional) Shapefile path or SpatVector.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH") to auto-download boundary.
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia") to filter.
#' @param timescale Character vector. Options: "annual", "seasonal", or c("annual", "seasonal").
#' @param seasons Named List. Defines the seasons and their months. 
#'    Example for Ethiopia: \code{list(Belg = 2:5, Kiremt = 6:9)}.
#'    Default is standard meteorological seasons (MAM, JJA, SON, DJF).
#' @param plot_result Logical. If TRUE, plots the results.
#'
#' @return A named list containing 'annual' and/or 'seasonal' SpatRasters.
#' @export
#' @importFrom terra rast vect crop mask time nlyr app time<- subset project crs
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot
calculate_precipitation_stats <- function(nc_input, start_date = NULL, 
                                          region = NULL,
                                          country = NULL, 
                                          region_name = NULL,
                                          timescale = c("annual", "seasonal"),
                                          seasons = list(
                                            MAM = 3:5, 
                                            JJA = 6:8, 
                                            SON = 9:11, 
                                            DJF = c(12, 1, 2)
                                          ),
                                          plot_result = TRUE) {
  
  # --- 1. Load Data ---
  if (is.character(nc_input)) {
    if (!file.exists(nc_input)) stop("NetCDF file not found.")
    r <- terra::rast(nc_input)
  } else {
    r <- nc_input
  }
  
  # --- 2. Handle Region & Country ---
  region_vect <- NULL
  if (!is.null(region)) {
    if (is.character(region)) region_vect <- terra::vect(region) else region_vect <- region
  } else if (!is.null(country)) {
    if (!requireNamespace("geodata", quietly = TRUE)) stop("Package 'geodata' required.")
    level_needed <- ifelse(is.null(region_name), 0, 1)
    region_vect <- tryCatch(geodata::gadm(country, level=level_needed, path=tempdir()), error=function(e) stop("Download failed."))
    if (!is.null(region_name)) {
      if (!(region_name %in% region_vect$NAME_1)) stop("Region name not found.")
      region_vect <- region_vect[region_vect$NAME_1 == region_name, ]
    }
  }
  
  if (!is.null(region_vect)) {
    if (terra::crs(region_vect) != terra::crs(r)) region_vect <- terra::project(region_vect, terra::crs(r))
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
  
  years <- as.numeric(format(dates, "%Y"))
  unique_years <- unique(years)
  results <- list()
  
  # --- 4. ANNUAL Calculation ---
  if ("annual" %in% timescale) {
    message("Calculating Annual Total Precipitation...")
    annual_layers <- list()
    valid_years <- c()
    
    for (yr in unique_years) {
      idx <- which(years == yr)
      if (length(idx) == 0) next
      
      # Sum all rainfall in the year
      yr_sum <- sum(terra::subset(r, idx), na.rm = TRUE)
      annual_layers[[length(annual_layers) + 1]] <- yr_sum
      valid_years <- c(valid_years, yr)
    }
    
    if (length(annual_layers) > 0) {
      results$annual <- terra::rast(annual_layers)
      names(results$annual) <- paste0("Annual_Total_", valid_years)
    }
  }
  
  # --- 5. SEASONAL Calculation ---
  if ("seasonal" %in% timescale) {
    message("Calculating Seasonal Precipitation...")
    # Loop through each season defined in the 'seasons' list
    for (season_name in names(seasons)) {
      target_months <- seasons[[season_name]]
      
      season_layers <- list()
      season_years <- c()
      
      message(paste("Processing Season:", season_name, "(Months:", paste(target_months, collapse=","), ")"))
      
      for (yr in unique_years) {
        # Find indices that match BOTH the Year AND the target Months
        month_nums <- as.numeric(format(dates, "%m"))
        idx <- which(years == yr & month_nums %in% target_months)
        
        if (length(idx) == 0) next
        
        # Calculate Seasonal Sum
        season_sum <- sum(terra::subset(r, idx), na.rm = TRUE)
        season_layers[[length(season_layers) + 1]] <- season_sum
        season_years <- c(season_years, yr)
      }
      
      if (length(season_layers) > 0) {
        s_stack <- terra::rast(season_layers)
        names(s_stack) <- paste0(season_name, "_", season_years)
        # Store in results list under the season name
        results[[season_name]] <- s_stack
      }
    }
  }
  
  # --- 6. Plotting ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Use "Blues 3" or "YlGnBu" for rainfall totals
    precip_pal <- hcl.colors(50, "Blues 3", rev = TRUE)
    
    # Count plots needed
    n_plots <- 0
    if (!is.null(results$annual)) n_plots <- n_plots + 1
    # Add number of seasons calculated
    seasonal_keys <- setdiff(names(results), "annual")
    n_plots <- n_plots + length(seasonal_keys)
    
    if (n_plots > 0) {
      # Dynamic layout (rows/cols) based on number of plots
      par(mfrow = c(ceiling(n_plots/2), 2), mar = c(3, 3, 3, 4))
      
      # Plot Annual
      if (!is.null(results$annual)) {
        terra::plot(results$annual[[1]], 
                    main = paste("Annual Total Precip - Year", unique_years[1]), 
                    col = precip_pal,
                    plg = list(title = "Total (mm)", title.cex = 0.8))
        if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE)
      }
      
      # Plot Seasons
      for (s_name in seasonal_keys) {
        # Format month names for the legend title
        m_names <- month.abb[seasons[[s_name]]]
        title_months <- paste(m_names[1], "-", tail(m_names, 1))
        
        terra::plot(results[[s_name]][[1]], 
                    main = paste(s_name, "Season Precip -", unique_years[1]), 
                    col = precip_pal,
                    plg = list(title = "Total (mm)", title.cex = 0.8))
        if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE)
      }
    }
  }
  
  return(results)
}
#' Calculate Mann-Kendall Trend Statistics (Annual & Seasonal)
#'
#' Applies the Mann-Kendall test and Sen's Slope estimator to detect long-term trends.
#' Utilizes parallel processing for large datasets.
#'
#' @param data_input A SpatRaster OR a List of SpatRasters. 
#' @param year_range Numeric Vector (Optional). Start and end year to analyze (e.g., c(1981, 2010)). 
#' @param significance_level Numeric. The threshold for statistical significance (Default 0.05).
#' @param cores Integer. Number of CPU cores to use. Default is 1. Set higher (e.g., 4 or 8) for large maps.
#' @param plot_result Logical. If TRUE, plots the significant trends.
#' @param region (Optional) Shapefile path or SpatVector to draw the boundary on the plot.
#' @param save_result Logical. If TRUE, saves the output as a .tif file.
#' @param output_dir Character. Directory to save the file.
#' @param output_name Character. Prefix for the saved file name.
#'
#' @return A List of SpatRasters.
#' @export
#' @importFrom terra app rast plot mask nlyr subset names vect project crs writeRaster
#' @importFrom trend mk.test sens.slope
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot mtext
calculate_trend_stats <- function(data_input, 
                                  year_range = NULL,
                                  significance_level = 0.05, 
                                  cores = 1,          
                                  plot_result = TRUE,
                                  region = NULL,                    # <--- NEW
                                  save_result = FALSE,              # <--- NEW
                                  output_dir = "Output_Maps",       # <--- NEW
                                  output_name = "Trend_Stats") {    # <--- NEW
  
  if (!requireNamespace("trend", quietly = TRUE)) stop("Package 'trend' is required.")
  
  if (inherits(data_input, "SpatRaster")) {
    data_list <- list(Dataset = data_input)
  } else if (is.list(data_input)) {
    data_list <- data_input
  } else {
    stop("Input must be a SpatRaster or a List of SpatRasters.")
  }
  
  # --- Handle Region for Plotting ---
  region_vect <- NULL
  if (!is.null(region)) {
    if (is.character(region)) {
      region_vect <- terra::vect(region)
    } else {
      region_vect <- region
    }
  }
  
  # --- Pixel-Wise Function ---
  calc_mk_pixel <- function(x) {
    if (any(is.na(x)) || length(unique(x)) < 3) return(c(NA, NA, NA)) 
    tryCatch({
      # Use :: to ensure parallel worker nodes find the package
      mk_res <- trend::mk.test(x)
      ss_res <- trend::sens.slope(x)
      return(c(mk_res$statistic, mk_res$p.value, ss_res$estimates))
    }, error = function(e) return(c(NA, NA, NA)))
  }
  
  output_list <- list()
  
  for (name in names(data_list)) {
    r_stack <- data_list[[name]]
    
    if (!is.null(year_range)) {
      layer_years <- as.numeric(gsub("\\D", "", names(r_stack)))
      target_idx <- which(layer_years >= year_range[1] & layer_years <= year_range[2])
      if (length(target_idx) == 0) next
      r_stack <- terra::subset(r_stack, target_idx)
    }
    
    message(paste("\nCalculating Trend for:", name))
    message(paste("Processing on", cores, "core(s). This may take a few minutes for large regions..."))
    
    # --- Parallel Execution ---
    # The 'cores' argument drastically reduces computation time
    trend_map <- terra::app(r_stack, fun = calc_mk_pixel, cores = cores)
    names(trend_map) <- c("Z_Score", "P_Value", "Sen_Slope")
    
    output_list[[name]] <- trend_map
  }
  
  # --- Plotting ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    trend_pal <- hcl.colors(50, "RdBu", rev = FALSE) 
    
    for (name in names(output_list)) {
      stats <- output_list[[name]]
      slope_layer <- stats[["Sen_Slope"]]
      
      if (dev.interactive()) dev.new()
      
      # Adjusted layout: 1 row, 1 column
      par(mfrow = c(1, 1), mar = c(3, 3, 3, 4), oma = c(0, 0, 2, 0))
      
      # Plot only the Sen's Slope Magnitude
      terra::plot(slope_layer, main = "Magnitude (Sen's Slope)", col = trend_pal)
      
      # --- Draw the Basin/Region Boundary ---
      if (!is.null(region_vect)) {
        if (terra::crs(region_vect) != terra::crs(slope_layer)) {
          region_vect <- terra::project(region_vect, terra::crs(slope_layer))
        }
        terra::plot(region_vect, add = TRUE, border = "black", lwd = 1.5)
      }
      
      mtext(paste("Trend Analysis:", name), outer = TRUE, cex = 1.2, font = 2)
    }
  }
  
  # --- Save Results (NEW BLOCK) ---
  if (save_result) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    for (name in names(output_list)) {
      # Creates names like "Trend_Stats_Dataset.tif" or "Trend_Stats_Annual.tif"
      out_path <- file.path(output_dir, paste0(output_name, "_", name, ".tif"))
      terra::writeRaster(output_list[[name]], out_path, overwrite = TRUE)
      message("Result saved to: ", out_path)
    }
  }
  
  return(output_list)
}
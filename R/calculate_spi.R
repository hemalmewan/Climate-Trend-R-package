#' Calculate Standardized Precipitation Index (SPI)
#'
#' Calculates the SPI for a given timescale (e.g., 1-month, 3-month, 12-month) 
#' using the SPEI package. Automatically converts daily input data to monthly totals.
#'
#' @details
#' SPI values are interpreted as:
#' \itemize{
#'   \item \strong{> 2.0}: Extremely wet
#'   \item \strong{1.5 to 1.99}: Very wet
#'   \item \strong{-0.99 to 0.99}: Near normal
#'   \item \strong{-1.0 to -1.49}: Moderately dry
#'   \item \strong{-1.5 to -1.99}: Severely dry
#'   \item \strong{< -2.0}: Extremely dry
#' }
#'
#' @param nc_input Character path to NetCDF file OR a SpatRaster object.
#' @param scale Integer. The time scale in months (e.g., 3 for SPI-3, 12 for SPI-12).
#' @param start_date Character/Date. Start date of the dataset (e.g., "1981-01-01").
#' @param region (Optional) Shapefile path or SpatVector.
#' @param country (Optional) Character. ISO3 code (e.g., "ETH").
#' @param region_name (Optional) Character. Admin region name (e.g., "Oromia").
#' @param cores Integer. Number of CPU cores for parallel processing. Default 1.
#' @param plot_result Logical. If TRUE, plots the final month of the calculated SPI.
#' @param save_result Logical. If TRUE, saves the output as a .tif file.
#' @param output_dir Character. Directory to save the file.
#' @param output_name Character. Prefix for the saved file name.
#'
#' @return A SpatRaster containing the monthly SPI values.
#' @export
#' @importFrom terra rast vect crop mask time nlyr app time<- tapp project crs writeRaster
#' @importFrom grDevices hcl.colors
#' @importFrom graphics par plot
calculate_spi <- function(nc_input, scale = 3, start_date = NULL, 
                          region = NULL, country = NULL, region_name = NULL, 
                          cores = 1, plot_result = TRUE,
                          save_result = FALSE,              # <--- NEW
                          output_dir = "Output_Maps",       # <--- NEW
                          output_name = "SPI_Index") {      # <--- NEW
  
  # --- 1. Check Dependencies ---
  if (!requireNamespace("SPEI", quietly = TRUE)) {
    stop("Package 'SPEI' is required. Please run: install.packages('SPEI')")
  }
  
  # --- 2. Load Data ---
  if (is.character(nc_input)) {
    if (!file.exists(nc_input)) stop("NetCDF file not found.")
    r <- terra::rast(nc_input)
  } else {
    r <- nc_input
  }
  
  # --- 3. Handle Region & Country ---
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
  
  # --- 4. Date Handling ---
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
  
  # --- 5. Aggregate Daily to Monthly (BULLETPROOF FIX) ---
  message("Aggregating daily precipitation to monthly totals...")
  
  ym_index <- format(dates, "%Y-%m")
  unique_months <- unique(ym_index)
  
  # Calculate in RAM first
  r_monthly_raw <- terra::tapp(r, index = ym_index, fun = sum, na.rm = TRUE)
  
  # FIX 1: Isolate from memory by writing explicitly to disk and reloading.
  # This guarantees all parallel worker cores can read the file safely.
  temp_monthly_file <- tempfile(fileext = ".tif")
  terra::writeRaster(r_monthly_raw, temp_monthly_file, overwrite = TRUE)
  r_monthly <- terra::rast(temp_monthly_file)
  
  # Ensure we have enough data (Optional safety check)
  if (terra::nlyr(r_monthly) < 36) {
    warning("SPI generally requires 30+ years of data for accurate distributions. You provided less than 3 years. Results may be NA or unreliable.")
  }
  
  # --- 6. Define the SPI Pixel Function ---
  message(paste("Calculating SPI-", scale, " (Scale: ", scale, " months)...", sep=""))
  
  calc_spi_pixel <- function(x, s) {
    if (all(is.na(x))) return(rep(NA, length(x)))
    
    # FIX 2: Prevent Mathematical Crashes
    # If a pixel has no variation (e.g., constantly 0), the Gamma fit crashes 
    # the underlying C code, instantly killing the parallel worker node.
    if (stats::var(x, na.rm = TRUE) == 0) return(rep(NA, length(x)))
    
    tryCatch({
      res <- SPEI::spi(x, scale = s, na.rm = TRUE)
      return(as.numeric(res$fitted))
    }, error = function(e) {
      return(rep(NA, length(x)))
    })
  }
  
  # --- 7. Apply to Raster ---
  spi_raster <- terra::app(r_monthly, fun = calc_spi_pixel, cores = cores, s = scale)
  names(spi_raster) <- paste0("SPI_", scale, "_", unique_months)
  
  # Cleanup temp file
  if (file.exists(temp_monthly_file)) unlink(temp_monthly_file)
  
  # --- 8. Plotting ---
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    spi_pal <- hcl.colors(50, "RdBu", rev = FALSE)
    par(mfrow = c(1, 1), mar = c(3, 3, 3, 4)) # Standardized margins
    
    last_layer_idx <- terra::nlyr(spi_raster)
    last_layer_name <- names(spi_raster)[last_layer_idx]
    
    terra::plot(spi_raster[[last_layer_idx]], 
                main = paste("Standardized Precipitation Index (", last_layer_name, ")", sep=""),
                col = spi_pal,
                range = c(-3, 3), 
                plg = list(title = "SPI Value", title.cex = 0.8))
    
    if (!is.null(region_vect)) terra::plot(region_vect, add = TRUE, border = "black")
  }
  
  # --- 9. Save Results (NEW BLOCK) ---
  if (save_result) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    if (inherits(spi_raster, "SpatRaster")) {
      out_path <- file.path(output_dir, paste0(output_name, ".tif"))
      terra::writeRaster(spi_raster, out_path, overwrite = TRUE)
      message("Result saved to: ", out_path)
    }
  }
  
  return(spi_raster)
}
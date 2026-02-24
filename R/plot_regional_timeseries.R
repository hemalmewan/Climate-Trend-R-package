#' Plot Regional Time Series of Climate Indices
#'
#' Averages the spatial pixel values of a raster stack over time and creates a 
#' time-series line or bar chart. Automatically adapts to monthly or annual stacks.
#'
#' @param raster_obj A SpatRaster containing multiple time layers.
#' @param stat Character. "mean", "max", "min", or "sum". How to aggregate the pixels. Default "mean".
#' @param title Character. The main title of the plot.
#' @param y_label Character. The Y-axis label.
#' @param plot_type Character. "line" for a line graph, "bar" for a bar chart. Default "line".
#' @param color Character. The color of the line or bars. Default "firebrick".
#'
#' @return A data.frame containing the extracted time series values (invisibly).
#' @export
#' @importFrom terra global nlyr
#' @importFrom graphics par plot axis grid barplot
plot_regional_timeseries <- function(raster_obj, 
                                     stat = "mean", 
                                     title = "Regional Time Series", 
                                     y_label = "Index Value",
                                     plot_type = "line",
                                     color = "firebrick") {
  
  if (!inherits(raster_obj, "SpatRaster")) stop("Input must be a SpatRaster object.")
  
  n_layers <- terra::nlyr(raster_obj)
  if (n_layers < 2) {
    warning("Raster only has 1 layer. Time series plots require 2 or more time steps.")
  }
  
  # --- 1. Extract Spatial Statistic ---
  message(paste("Calculating spatial", stat, "across all pixels for each layer..."))
  
  # terra::global calculates the stat (e.g., mean) for the whole map layer
  ts_data <- terra::global(raster_obj, fun = stat, na.rm = TRUE)
  
  # --- 2. Prepare Data ---
  y_values <- ts_data[[1]]
  x_names <- rownames(ts_data)
  
  # --- 3. Plotting ---
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Increase bottom margin slightly in case names are long
  par(mar = c(6, 4, 4, 2) + 0.1) 
  
  if (plot_type == "bar") {
    # Draw Bar Chart (Great for PRCPTOT or Count indices)
    bp <- barplot(y_values, names.arg = x_names, col = color, 
                  main = title, ylab = y_label, las = 2, border = NA)
    grid(nx = NA, ny = NULL)
    
  } else {
    # Draw Line Chart (Great for SPI or Percentiles)
    x_idx <- 1:n_layers
    plot(x_idx, y_values, type = "b", col = color, lwd = 2, pch = 16,
         main = title, ylab = y_label, xlab = "", xaxt = "n")
    
    # las = 2 turns the x-axis labels vertically so they don't overlap
    axis(1, at = x_idx, labels = x_names, las = 2) 
    grid()
  }
  
  # --- 4. Return Data ---
  # Silently return the raw dataframe in case you want to save it to a CSV later
  df_out <- data.frame(Time = x_names, Value = y_values)
  invisible(df_out)
}
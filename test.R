#### ===================================================
## Import ClimateTrend 
library(ClimateTrend)
library(parallel)
##=====================================================
## =======================================
## Help Function
##=======================================
help("calculate_prcptot")
help("calculate_cdd")
help("calculate_cwd")
help("calculate_r95p")
help("calculate_r95ptot")
help("calculate_r99p")
help("calculate_r99ptot")
help("calculate_rnnmm")
help("calculate_rxdday")
help("calculate_precipitation_stats")
help("calculate_spi")
help("calculate_trend_stats")
help("plot_regional_timeseries")
help("pre_whiten")
help("calculate_spatial_wsds")

##======================================
##======================================
##======================================


##============================================================
## Man Kendell Trend Detection For each Climate Indices
##=============================================================

# 1. Your list of files
nc_files <- c(
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1981.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1982.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1983.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1984.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1985.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1986.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1987.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1988.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1989.nc",
  "C:/Hydrology-Project/Rainfall Trend/NCDF/ethiopia/Daily_nc_1990.nc"
  
)
my_cores <- parallel::detectCores() - 1 #computer corses (default 23)


##======================================================
## Consecutive Dry Days(CDD)
##======================================================

cdd_results <- list()

message("Starting CDD processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  cdd_results[[file_index]] <- calculate_cdd(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    country = "ETH",
    region_name = "Amhara",
    threshold = 1,
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(cdd_results)[file_index] <- paste0("CDD_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
cdd_10yr_stack <- terra::rast(lapply(cdd_results, function(x) x$annual))

## Apply priwheting for a specific location
lon_lat <- cbind(37.5, 11.5) 
pixel_10yr_cdd <- terra::extract(cdd_10yr_stack, lon_lat)
my_vector <- as.numeric(pixel_10yr_cdd[-1]) 


 


trend_results <- calculate_trend_stats(
   data_input = cdd_10yr_stack,
   significance_level = 0.05,
   cores = my_cores,
   plot_result = TRUE,
 )


# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))                                                                                      


##=======================================
## Consecutive Wet Days(CWD)
##=======================================

cwd_results <- list()

message("Starting CWD processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  cwd_results[[file_index]] <- calculate_cwd(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    country = "ETH",
    region_name = "Amhara",
    threshold = 1,
    plot_result = FALSE,
    save_result = TRUE,
    output_dir ="C:/Hydrology-Project/",
    output_name = paste0("CWD_",year_string)
  )
  
  # Name the list element so it's easy to find later
  names(cwd_results)[file_index] <- paste0("CWD_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
cwd_10yr_stack <- terra::rast(lapply(cwd_results, function(x) x$annual))



trend_results <- calculate_trend_stats(
  data_input = cwd_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))                                                                                      

##==========================================
## PRCPTOT
##==========================================

prcptot_results <- list()

message("Starting PRCPTOT processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  prcptot_results[[file_index]] <- calculate_prcptot(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    country = "ETH",
    region_name = "Amhara",
    threshold = 1,
    plot_result = FALSE,
    timescale ="annual",
    save_result = TRUE,
    output_dir ="C:/Hydrology-Project/",
    output_name =paste0("PRCPTOT_",year_string)
  )
  
  # Name the list element so it's easy to find later
  names(prcptot_results)[file_index] <- paste0("PRCPTOT_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
prcptot_10yr_stack <- terra::rast(lapply(prcptot_results, function(x) x$annual))


trend_results <- calculate_trend_stats(
  data_input = prcptot_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration)) 

##=============================================
## Rxdday
##=============================================
rxdday_results <- list()

message("Starting Rxdday processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  rxdday_results[[file_index]] <- calculate_rxdday(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    country = "ETH",
    region_name = "Amhara",
    window_size =5,
    timescale ="annual",
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(rxdday_results)[file_index] <- paste0("Rxdday_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
rxdday_10yr_stack <- terra::rast(lapply(rxdday_results, function(x) x$annual))


trend_results <- calculate_trend_stats(
  data_input = rxdday_10yr_stack ,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration)) 

##=================================================
## Rnnmm
##=================================================
rnnmm_results <- list()

message("Starting Rnnmm processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
   rnnmm_results[[file_index]] <- calculate_rnnmm(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    country = "ETH",
    region_name = "Amhara",
    threshold = 1,
    timescale = "annual",
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(rnnmm_results)[file_index] <- paste0("R1mm_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
rnnmm_10yr_stack <- terra::rast(lapply(rnnmm_results, function(x) x$annual))


trend_results <- calculate_trend_stats(
  data_input = rnnmm_10yr_stack ,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration)) 

##=======================================
## R95
##======================================

r95p_results <- list()

message("Starting R95p processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  r95p_results[[file_index]] <- calculate_r95p(
    nc_input = current_file,
    baseline_range = c(1981,2010),
    start_date = dynamic_start_date,
    percentile = 0.95,  
    country = "ETH",
    region_name = "Amhara",
    wet_threshold =1,
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(r95p_results)[file_index] <- paste0("R95p_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
r95p_10yr_stack <- terra::rast(r95p_results)



trend_results <- calculate_trend_stats(
  data_input = r95p_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))

##=========================================
## R95pTOT
##=====================================
r95ptot_results <- list()

message("Starting R95pTOT processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  r95ptot_results[[file_index]] <- calculate_r95ptot(
    nc_input = current_file,
    baseline_range = c(1981,2010),
    start_date = dynamic_start_date,
    country = "ETH",
    region_name = "Amhara",
    wet_threshold =1,
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(r95ptot_results)[file_index] <- paste0("R95pTOT_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
r95ptot_10yr_stack <- terra::rast(r95ptot_results)


trend_results <- calculate_trend_stats(
  data_input = r95ptot_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))

##==========================
## R99p
##=========================
r99p_results <- list()

message("Starting R99p processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  r99p_results[[file_index]] <- calculate_r99p(
    nc_input = current_file,
    baseline_range = c(1981,2010),
    start_date = dynamic_start_date,
    country = "ETH",
    region_name = "Amhara",
    wet_threshold =1,
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(r99p_results)[file_index] <- paste0("R99p_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
r99p_10yr_stack <- terra::rast(r99p_results)


trend_results <- calculate_trend_stats(
  data_input = r99p_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()



# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))


trend_results <- calculate_trend_stats(
  data_input = r99ptot_10yr_stack ,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))



##==========================================
## Standard Precipitation Index(SPI)
##==========================================

message("Merging 10 years of daily data............")
r_list <- lapply(nc_files, terra::rast)
full_data <- terra::rast(r_list)

# CORRECTED: Use nlyr() to get the total number of layers (days) as an integer
total_days <- terra::nlyr(full_data)

# Now seq() gets a proper number for length.out
terra::time(full_data) <- seq(from = as.Date("1981-01-01"), by = "day", length.out = total_days)

## calculate SPI 
message("Calculating SPI-3...")

spi3_results <- calculate_spi(
  nc_input = full_data,
  scale = 3,
  start_date = "1981-01-01",
  country = "ETH",
  cores = my_cores,
  plot_result = TRUE,
  save_result = TRUE,
  output_dir = "C:/Hydrology-Project/",
  output_name = "SPI_INDEX"
)

# Let's plot August 1984 to see the drought maps!
dec_1985_layer <- which(names(spi3_results) == "SPI_3_1985-12")

if(length(dec_1985_layer) > 0) {
  # Added the margin fix here so your legend doesn't get cut off
  par(mar = c(3, 3, 3, 6)) 
  
  # EXPLICIT FIX: Use terra::plot instead of just plot
  terra::plot(spi3_results[[dec_1985_layer]], main = "SPI-3 for December 1985")
  
} else {
  message("Layer not found. Run 'names(spi3_results)' to see available dates.")
}

##=======================================
## Wet Season Dry Spell
##======================================
WSDS_results <- list()

message("Starting WSDS processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  WSDS_results[[file_index]] <- calculate_spatial_wsds(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    country = "ETH",
    region_name = "Amhara",
    season_name = "Kiremt",
    threshold = 1,
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(WSDS_results)[file_index] <- paste0("WSDS_", year_string)
}

message("\nStacking the 10-year WSDS maps...")
WSDS_10yr_stack <- terra::rast(WSDS_results)

message("Running Mann-Kendall Trend Analysis...")
trend_results <- calculate_trend_stats(
  data_input = WSDS_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)


# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))                                                                                      


##============================================================
## Time Series Plot 
##==========================================

##==================================
## CDD
##=================================
# 2. Plot the monthly time series as a Bar Chart
plot_regional_timeseries(
  raster_obj =cdd_10yr_stack , 
  stat = "max",                   
  title = "Consecutive Dry Days over the 1981-1990",
  y_label = "Consecutive Dry Days",
  plot_type = "line",               
  color = "steelblue"
)


##===============================
## cwd
##===============================
plot_regional_timeseries(
  raster_obj =cwd_10yr_stack, 
  stat = "max",                  
  title = "Consecutive Wet Days over the 1981-1990",
  y_label = "Consecutive Wet Days",
  plot_type = "line",              
  color = "steelblue"
)


##============================
## PRCPTOT
##===========================
plot_regional_timeseries(
  raster_obj =prcptot_10yr_stack, 
  stat = "max",                  
  title = "PRCPTOT over the 1981-1990",
  y_label = "PRCPTOT(mm)",
  plot_type = "line",              
  color = "steelblue"
)

##=========================
## Rxdday
##========================
plot_regional_timeseries(
  raster_obj =rxdday_10yr_stack, 
  stat = "max",                  
  title = "Rxdday over the 1981-1990",
  y_label = "Rxdday(mm)",
  plot_type = "line",              
  color = "steelblue"
)

##========================
## Rnnmm
##=======================
plot_regional_timeseries(
  raster_obj =rnnmm_10yr_stack, 
  stat = "max",                  
  title = "Rnnmm over the 1981-1990",
  y_label = "Rnnmm(Days)",
  plot_type = "line",              
  color = "steelblue"
)

##===================
## R95pToT
##===================
plot_regional_timeseries(
  raster_obj =r95ptot_10yr_stack, 
  stat = "max",                  
  title = "R95pTOT over the 1981-1990",
  y_label = "R95pTOT(%)",
  plot_type = "line",              
  color = "steelblue"
)

##======================
## R99p
##=====================
plot_regional_timeseries(
  raster_obj =r99p_10yr_stack, 
  stat = "max",                  
  title = "R99p over the 1981-1990",
  y_label = "R99p(mm)",
  plot_type = "line",              
  color = "steelblue"
)


##===============================
## WSDS
##===============================
plot_regional_timeseries(
  raster_obj =WSDS_10yr_stack, 
  stat = "max",                  
  title = "Consecutive Wet Season Dry Spell(WSDS) the 1981-1990",
  y_label = "Consecutive Dry Days in Wet Season",
  plot_type = "line",              
  color = "steelblue"
)





##====================================
## Basins
##====================================
# Load your new basin shapefile
all_basins <- terra::vect("C:/Hydrology-Project/Shapefiles/Ethiopia_Major_Basins.shp")



# Select a specific basin (You can inspect 'all_basins' to find the specific ID or area you want)
# For example, selecting the 5th basin in the list:
my_target_basin <- all_basins[26, ] 

##==================
## CDD
##==================
basin_cdd_results <- list()

message("Starting CDD processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  basin_cdd_results[[file_index]] <- calculate_cdd(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    region = my_target_basin,
    threshold = 1,
    plot_result = FALSE
  )
  
  # Name the list element so it's easy to find later
  names(basin_cdd_results)[file_index] <- paste0("CDD_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
basin_cdd_10yr_stack <- terra::rast(lapply(basin_cdd_results, function(x) x$annual))

 
trend_results <- calculate_trend_stats(
  data_input = basin_cdd_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  region = my_target_basin,
  plot_result = TRUE
)

##Time series plot for 10 year period in selected basin
plot_regional_timeseries(
  raster_obj =basin_cdd_10yr_stack, 
  stat = "max",                  
  title = "CDD of Basin over the 1981-1990",
  y_label = "CDD(Days)",
  plot_type = "line",              
  color = "steelblue"
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))


##==================================================
## Seasonal and Annual Climate Trend Analysis
##==================================================

##===========================
## 1. Annual Trend Result
##===========================
annual_results <- list()

message("Starting Annual Precipitation processing...")
start_time <- Sys.time() ## start time

for (file_index in 1:length(nc_files)) {
  
  # Get the current file path
  current_file <- nc_files[file_index]
  
  # Extract the 4-digit year from the file path using Regular Expressions
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  
  # Build the dynamic start date (e.g., "1981-01-01")
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # 3. Calculate and store in the list (Notice the double brackets [[ ]])
  annual_results[[file_index]] <- calculate_precipitation_stats(
    nc_input = current_file,
    start_date = dynamic_start_date, # Passes the dynamic date here
    country = "ETH",
    region_name = "Amhara",
    plot_result = FALSE,
    timescale = "annual"

  )
  
  # Name the list element so it's easy to find later
  names(annual_results)[file_index] <- paste0("Annual_", year_string)
}

# Extract just the 'annual' raster from each list item and combine them into one stack
annual_10yr_stack <- terra::rast(lapply(annual_results, function(x) x$annual))


trend_results <- calculate_trend_stats(
  data_input = annual_10yr_stack,
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

# --- STOP THE CLOCK ---
end_time <- Sys.time()


# Calculate and print the total duration
total_duration <- round(end_time - start_time, 2)
message("\n--- Processing Complete ---")
message("Total time taken: ", format(total_duration))

##===============================
## 2.Seasonal Trend Result
##==============================

## Define seasons
# 1. Define the Ethiopian seasons as a NAMED list
# Your function strictly requires this list to have names!
eth_seasons <- list(
  Kiremt = 6:9,
  Belg = 2:5
)

# 2. Initialize empty lists to store the rasters for each season
Kiremt_list <- list()
Belg_list <- list()

message("Starting Seasonal Precipitation processing...")
start_time <- Sys.time() 

# 3. Loop through your 10 years of NetCDF files
for (file_index in 1:length(nc_files)) {
  
  current_file <- nc_files[file_index]
  year_string <- regmatches(current_file, regexpr("[0-9]{4}", current_file))
  dynamic_start_date <- paste0(year_string, "-01-01")
  
  message("\nProcessing Year: ", year_string)
  
  # Run your function and pass the ENTIRE eth_seasons list
  yearly_result <- calculate_precipitation_stats(
    nc_input = current_file,
    start_date = dynamic_start_date,
    country = "ETH",
    region_name = "Amhara",
    timescale = "seasonal",    # Tells it to only do seasons
    seasons = eth_seasons,     # Feeds both Kiremt and Belg into the function
    plot_result = FALSE
  )
  
  
  # Extract them and put them into our master lists
  Kiremt_list[[file_index]] <- yearly_result$Kiremt
  Belg_list[[file_index]] <- yearly_result$Belg
}

# --- 4. Stack the Results ---
# Now we use terra::rast to stack the bare rasters directly!
message("\nStacking Kiremt and Belg 10-Year Maps...")
Kiremt_10yr_stack <- terra::rast(Kiremt_list)
Belg_10yr_stack <- terra::rast(Belg_list)

# --- 5. Run the Trend Analysis ---
message("Running Mann-Kendall Trend Analysis for Kiremt Season...")
kiremt_trend_results <- calculate_trend_stats(
  data_input = Kiremt_10yr_stack, # Pass the Kiremt stack
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

end_time <- Sys.time()
message("Total time taken: ", format(round(end_time - start_time, 2)))


# --- 5. Run the Trend Analysis ---
message("Running Mann-Kendall Trend Analysis for Belg Season...")
kiremt_trend_results <- calculate_trend_stats(
  data_input = Belg_10yr_stack, # Pass the Kiremt stack
  significance_level = 0.05,
  cores = my_cores,
  plot_result = TRUE
)

## 📖 User Manual: Climate Trend Analysis & Index Toolkit
Welcome to the Climate Trend Analysis Toolkit!This R package is designed to process massive NetCDF climate datasets,calculate standard ETCCDI precipitation and drought indices,and perform robudst Mann-Kendall trend analysis.

Wheather you are analyzing a single country or a specific hydrological river basin,this toolkit handles the heavy lifting,parrallel processing,and spatial masking for you.

## 📑 Table of Contents

1.Core Concepts & Golden Rules
2.Data Preparation & Masking
3.Workflow A: Counting Indices (Year-by-Year)
4.Workflow B: Percentile & Multi-Month Indices (Merged Data)
5.Trend Analysis (Mann-Kendall & Sen's Slope)
6.Visualizing Time Series
7.Function Glossary

## 1.Core Concepts & Golden Rules
Before you start writing code,you must understand the two types of climate indices in this package.Treating them the same way will result in scientifically incorrent data.

### Type 1:Counting/Absolute Indices (PRCPTOT,CDD,CWD,Rnnmm,Rxdday)
    How they work:They look at a single year in isolation and count something(e.g.,"How many dry days were in 1985?").

    How to process:You can process these inside a for loop,feeding the function one yearly .nc file at a time.

## Type 2:Percentile & Distribution Indices(R95p,R99p,R95pTOT,SPI)
    How they work:They compare a single day or month against massive historical baseline(e.g., "Was today's rain extreme compared to the last 30 years?")

    The Golden Rule("Merge First"):You MUST merge your NetCDF file into a single,multi-year dataset BEFORE running these functions.Never run these inside a year-by-year loop.

## 2.Data Preparation & Masking
This package is optimized to use the **terra** package for spatial data.To save memory(RAM) and speed up caclulations,always clip your data to your specific study area before doing the math.

### Auto-Donwloading Country Boundaries
If you want to analyze a country  or state,the toolkit can download the boundary for you automatically using the **country** and **region_name** arguments.


library(terra)
data <- rast("Daily_nc_1981.nc")

#### Calculates total rainfall just for the Amhara region of Ethiopia
result <- calculate_prcptot(
  nc_input = data,
  country = "ETH",
  region_name = "Amhara" 
)

### Using Custom Shapefiles (River Basins)
If you are studying a specific watershed,load your shapefile via **terra::vect()** and pass it to the **region** argument.

my_basin <- terra::vect("Path/To/Awash_Basin.shp")

#### The function will automatically crop the data to your exact basin shape
result <- calculate_prcptot(
  nc_input = data,
  region = my_basin 
)

## 3.Workflow A:Counting Indices(Year-by-Year)
For basic indices,you can loop through your NetCDF files,calculate the index,and let the package automatically save the **.tif** maps to your hard drive.

nc_files <- list.files("Data/NCDF/", pattern = "\\.nc$", full.names = TRUE)
my_basin <- terra::vect("Shapefiles/Awash_Basin.shp")

#### Process CDD year by year
for (file in nc_files) {
  #### Extract year for naming (e.g., "1981")
  year_str <- regmatches(file, regexpr("[0-9]{4}", file))
  start_date <- paste0(year_str, "-01-01")
  
  calculate_cdd(
    nc_input = file,
    start_date = start_date,
    region = my_basin,
    plot_result = FALSE, 
    save_result = TRUE,                      # Automatically saves the map!
    output_dir = "Results/CDD",
    output_name = paste0("CDD_", year_str)   # Saves as CDD_1981.tif
  )
}

## 4.Workflow B:Percentile & Multi-Month Indices(Merged Data)
For complex indices like **R95p** or **SPI-3**,you must merge the dataset first so the function can calculate the historical baseline.

#### 1.Merge files
nc_files <- list.files("Data/NCDF/", pattern = "\\.nc$", full.names = TRUE)
full_data <- terra::rast(nc_files)

#### Fix the timeline
total_days <- terra::nlyr(full_data)
terra::time(full_data) <- seq(from = as.Date("1981-01-01"), by = "day", length.out = total_days)

#### 2.Crop early to save memory
my_basin <- terra::vect("Shapefiles/Awash_Basin.shp")
full_data <- terra::crop(full_data, terra::project(my_basin, terra::crs(full_data)), mask = TRUE)

#### 3.Calculate R95p for the whole 30-year stack at once
r95p_stack <- calculate_r95p(
  nc_input = full_data,
  baseline_range = c(1981, 2010),
  region = my_basin,
  save_result = TRUE,
  output_dir = "Results/R95p",
  output_name = "R95p_30yr_Stack"
)

## 5.Trend Analysis (Mann-Kendall & Sen's Slope)
Once you have a multi-year stack of an index(e.g.,30 layers of CDD),you can run a pixel-by-pixel trend analysis.This uses parallel processing to speed up the heavy math.

#### Pass your multi-year stack into the trend function
trends <- calculate_trend_stats(
  data_input = r95p_stack,
  significance_level = 0.05,
  cores = 4,                   # Utilize 4 CPU cores for speed
  region = my_basin,           # Draws the basin border on the final plot
  save_result = TRUE,
  output_dir = "Results/Trends",
  output_name = "R95p_Trend"
)

#### The result is a list containing Z_Score, P_Value, and Sen_Slope layers

## 6.Visualizing Time Series
To create regional averages(converting a 2D map into a 1D line or bar chart over time),use the build-in time series tool.

#### Plots the spatial average of your R95p stack over 30 years
ts_data <- plot_regional_timeseries(
  raster_obj = r95p_stack,
  stat = "max",                   # Aggregates pixels using the Mean
  title = "Maximum Extreme Rain (R95p) for the Awash Basin",
  y_label = "Precipitation (mm)",
  plot_type = "line",              # Use "line" or "bar"
  color = "firebrick"
)

#### You can even save the raw numerical data to a CSV!
write.csv(ts_data, "Awash_R95p_Timeseries.csv", row.names = FALSE)

## 7.Function Glossary
### Precipiptation & Extremes
    calculate_prcptot():Total annual/monthly precipitation from wet days(>=1mm).
    calculate_rnnmm():Count of days excedding a threshold(e.g.,R10mm,R20mm).
    calculate_rxdday():Max rainfall over a rolloing window(e.g.,Rx1day,Rx5day).
    calculate_r95p()/calculate_r99p():Total rain from extreme(96th/99th percentile) days.
    calculate_r95ptot()/calculate_r99ptot():Percentage contribution of extreme days to the annual total.

### Drought & Spells
    calculate_cdd():Consecutive Dry Days(longest streak of days<1mm).
    calculate_cwd():Consective Wet Days(longest streak of days>=1mm).
    calculate_spi():Standardized Precipitation Index(meteorological drought at 1,3,6, or 12-month scales).

### Analysis & Plotting
    calculate_trend_stats():Pixel-by-pixel Mann-Kendall significance and Sen's Slope magnitude.
    plot_regional_timeseries():Spatial averaging to genertes regional time-series line/bar graphs.


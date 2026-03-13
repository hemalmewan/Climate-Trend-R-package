#' Detect and Remove Serial Correlation via Pre-Whitening
#'
#' @description
#' Evaluates a time series vector for significant Lag-1 autocorrelation. If the 
#' autocorrelation exceeds the 95% confidence interval, it applies a standard 
#' Pre-Whitening (PW) transformation to scrub the memory from the dataset prior 
#' to trend analysis.
#'
#' @details
#' Hydrological and climate time series often exhibit positive serial correlation, 
#' which artificially inflates the significance of the Mann-Kendall trend test 
#' (Type I Error). This function automatically detects and removes it.
#' 
#' The significance bound is calculated as:
#' \deqn{B = 1.96 / \sqrt{n}}
#' 
#' If the Lag-1 autocorrelation (\eqn{r_1}) exceeds \eqn{B}, the series is pre-whitened:
#' \deqn{X'_t = X_t - r_1 X_{t-1}}
#' Note that pre-whitening reduces the length of the time series by exactly one time step.
#'
#' @param x A numeric vector representing a sequential time series (e.g., annual CDD values).
#' @param plot Logical. If `TRUE`, generates Autocorrelation Function (ACF) plots. Default is `FALSE`.
#' @param verbose Logical. If `TRUE`, prints messages detailing the detection and scrubbing process. Default is `FALSE`.
#'
#' @return A numeric vector. If significant autocorrelation was detected, returns 
#' the pre-whitened series (length `n - 1`). If no significant autocorrelation 
#' was found, returns the original series (length `n`).
#'
#' @export
#'
#' @examples
#' # Generate dummy climate data with artificial serial correlation
#' set.seed(42)
#' simulated_cdd <- as.numeric(arima.sim(model = list(ar = 0.7), n = 30) + 20)
#' 
#' # Clean the data with logging and plotting turned ON
#' clean_cdd <- pre_whiten(simulated_cdd, plot = TRUE, verbose = TRUE)
#' 
pre_whiten <- function(x, plot = FALSE, verbose = FALSE) {
  # Ensure input is numeric and handle missing values
  if (!is.numeric(x)) stop("Input data must be a numeric vector.")
  n <- length(x)
  if (n < 4) stop("Time series is too short to calculate meaningful autocorrelation.")
  
  if (verbose) message("--- Pre-Whitening Check Initiated ---")
  
  # Calculate the Lag-1 Autocorrelation coefficient (r1)
  # Only plot here if the user requested it, with a custom title
  if (plot) {
    acf_result <- stats::acf(x, plot = TRUE, na.action = stats::na.pass, main = "Original Series ACF")
  } else {
    acf_result <- stats::acf(x, plot = FALSE, na.action = stats::na.pass)
  }
  
  r1 <- acf_result$acf[2] 
  
  # Calculate the 95% confidence bound
  conf_bound <- 1.96 / sqrt(n)
  
  if (verbose) {
    message(sprintf("Lag-1 Autocorrelation (r1): %.3f", r1))
    message(sprintf("95%% Confidence Bound:     +/- %.3f", conf_bound))
  }
  
  # Check if the Lag-1 correlation exceeds the confidence bound
  if (abs(r1) > conf_bound) {
    
    if (verbose) message("⚠️ Significant memory detected. Applying Pre-Whitening (PW) formula...")
    
    # Initialize an empty vector for the scrubbed data (n - 1)
    x_clean <- numeric(n - 1)
    
    # Apply the Pre-Whitening formula
    for (t in 2:n) {
      x_clean[t - 1] <- x[t] - (r1 * x[t - 1])
    }
    
    if (plot) {
      # Plot the newly scrubbed data to visually prove the memory is gone
      stats::acf(x_clean, plot = TRUE, na.action = stats::na.pass, main = "Pre-Whitened Series ACF (Memory Removed)")
    }
    
    if (verbose) message("✅ Pre-Whitening complete. Returning scrubbed series (length n-1).")
    return(x_clean)
    
  } else {
    if (verbose) message("✅ No significant memory detected. Returning original series (length n).")
    return(x)
  }
}
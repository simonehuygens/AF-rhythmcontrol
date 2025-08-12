rm(list = ls())
options(scipen = 999)
df <- read.csv(here::here("input", "bedragen_per_periode_v  v2 20250429.csv"), sep = ",", header = TRUE)

# Rename variables
df$c_SAF          <- df$v_p0
df$c_before_CA    <- df$v_p1
df$c_after_CA1    <- df$v_p2 
df$c_after_CA2    <- df$v_p3
df$c_SFAF         <- df$v_p4

# Calculate differences
df$c_d_before_CA  <- df$c_before_CA - df$c_SFAF 
df$c_d_after_CA1  <- df$c_after_CA1 - df$c_SFAF
df$c_d_after_CA2  <- df$c_after_CA2 - df$c_SFAF
df$c_d_SAF        <- df$c_SAF - df$c_SFAF

# Subset data of patients without events
df_no_ev_24 <- df[df$ev_24 == "n" , ] 
# Only include 2021 because there seems to be impact of COVID 
df_no_ev_24 <- df[df$ev_24 == "n" & df$jaar_KA_n == "2021", ]
# Exclude patient with more than one DBC of catheter ablation
df_no_ev_24 <- df_no_ev_24[df_no_ev_24$v_pr_khm_p1 < 20000, ]

# Function to find the best fitting distribution for a variable
find_best_distribution <- function(data, distributions = c("norm", "lnorm", "exp", "gamma", "weibull"), 
                                   plot_results = TRUE, return_details = FALSE) {
  # Load required packages
  if (!require(fitdistrplus)) install.packages("fitdistrplus")
  library(fitdistrplus)
  if (!require(ggplot2)) install.packages("ggplot2")
  library(ggplot2)
  
  # Ensure data is numeric and remove NA values
  data <- as.numeric(data)
  data <- data[!is.na(data)]
  
  if (length(data) == 0) {
    stop("No valid data after cleaning")
  }
  
  # Store fitting results
  fit_results <- list()
  aic_values <- numeric(length(distributions))
  bic_values <- numeric(length(distributions))
  
  # Fit each distribution and calculate AIC/BIC
  for (i in seq_along(distributions)) {
    dist_name <- distributions[i]
    tryCatch({
      fit <- fitdist(data, dist_name)
      fit_results[[dist_name]] <- fit
      aic_values[i] <- fit$aic
      bic_values[i] <- fit$bic
      cat(sprintf("Fitted %s distribution: AIC = %.4f, BIC = %.4f\n", dist_name, fit$aic, fit$bic))
    }, error = function(e) {
      cat(sprintf("Failed to fit %s distribution: %s\n", dist_name, e$message))
      aic_values[i] <- Inf
      bic_values[i] <- Inf
    })
  }
  
  # Find best distribution based on AIC
  valid_indices <- which(is.finite(aic_values))
  if (length(valid_indices) == 0) {
    stop("Could not fit any of the specified distributions")
  }
  
  best_aic_index <- valid_indices[which.min(aic_values[valid_indices])]
  best_dist_aic <- distributions[best_aic_index]
  
  # Find best distribution based on BIC
  best_bic_index <- valid_indices[which.min(bic_values[valid_indices])]
  best_dist_bic <- distributions[best_bic_index]
  
  # Get the best fitting distribution (using AIC as default)
  best_fit <- fit_results[[best_dist_aic]]
  
  # Plot results if requested
  if (plot_results && length(valid_indices) > 0) {
    
    # Create a density plot comparing all fitted distributions
    densities <- data.frame(x = numeric(0), y = numeric(0), Distribution = character(0))
    x_range <- seq(min(data) - 0.1 * sd(data), max(data) + 0.1 * sd(data), length.out = 500)
    
    for (dist_name in names(fit_results)) {
      fit <- fit_results[[dist_name]]
      params <- fit$estimate
      
      # Calculate density based on distribution type
      if (dist_name == "norm") {
        dens <- dnorm(x_range, mean = params["mean"], sd = params["sd"])
      } else if (dist_name == "lnorm") {
        dens <- dlnorm(x_range, meanlog = params["meanlog"], sdlog = params["sdlog"])
      } else if (dist_name == "exp") {
        dens <- dexp(x_range, rate = params["rate"])
      } else if (dist_name == "gamma") {
        dens <- dgamma(x_range, shape = params["shape"], rate = params["rate"])
      } else if (dist_name == "weibull") {
        dens <- dweibull(x_range, shape = params["shape"], scale = params["scale"])
      } else {
        next
      }
      
      densities <- rbind(densities, data.frame(
        x = x_range,
        y = dens,
        Distribution = rep(dist_name, length(x_range))
      ))
    }
    
    # Reset plot layout
    par(mfrow = c(1, 1))
    
    # Create a data histogram with density curves
    hist(data, freq = FALSE, breaks = 30, 
         xlab = "Value", border = "darkgray", col = "lightgray")
    
    colors <- c("red", "blue", "green", "purple", "orange", "brown", "black")
    line_types <- c(1, 2, 3, 4, 5, 6, 1)
    
    legend_text <- character(0)
    legend_colors <- c()
    legend_lines <- c()
    
    for (i in seq_along(names(fit_results))) {
      dist_name <- names(fit_results)[i]
      dist_data <- subset(densities, Distribution == dist_name)
      
      if (nrow(dist_data) > 0) {
        lines(dist_data$x, dist_data$y, col = colors[i], lty = line_types[i], lwd = 2)
        
        # Format AIC for legend
        aic_value <- fit_results[[dist_name]]$aic
        legend_entry <- sprintf("%s (AIC: %.2f)", dist_name, aic_value)
        
        # Mark the best fit
        if (dist_name == best_dist_aic) {
          legend_entry <- paste(legend_entry, "- BEST FIT")
        }
        
        legend_text <- c(legend_text, legend_entry)
        legend_colors <- c(legend_colors, colors[i])
        legend_lines <- c(legend_lines, line_types[i])
      }
    }
    
    legend("topright", legend = legend_text, col = legend_colors, 
           lty = legend_lines, cex = 0.8, bty = "n")
  }
  
  # Summary of results
  cat("\n=========== SUMMARY ===========\n")
  cat(sprintf("Best distribution by AIC: %s (AIC = %.4f)\n", best_dist_aic, min(aic_values[valid_indices])))
  cat(sprintf("Best distribution by BIC: %s (BIC = %.4f)\n", best_dist_bic, min(bic_values[valid_indices])))
  cat("Parameters of best fitting distribution (AIC):\n")
  print(best_fit$estimate)
  
  # Return results
  if (return_details) {
    return(list(
      best_distribution_aic = best_dist_aic,
      best_distribution_bic = best_dist_bic,
      best_fit = best_fit,
      all_fits = fit_results,
      aic_values = aic_values,
      bic_values = bic_values,
      distributions = distributions
    ))
  } else {
    return(best_dist_aic)
  }
}

# Determine best fitting distribution (weibull)
best_dist <- find_best_distribution(df_no_ev_24$v_pr_khm_p1, 
                                    distributions = c("norm", "lnorm", "exp", "gamma", "weibull"),
                                    plot_results = TRUE,
                                    return_details = TRUE)

c_CA_DBC_params  <- fitdist(df_no_ev_24$v_pr_khm_p1, "weibull")

# Determine best fitting distribution (lnorm)
best_dist <- find_best_distribution(df_no_ev_24$v_farm_e_aad_p0, 
                                    distributions = c("norm", "lnorm", "exp", "gamma", "weibull"),
                                    plot_results = TRUE,
                                    return_details = TRUE)

c_AAD_params    <- fitdist(df_no_ev_24$v_farm_e_aad_p0, "lnorm")

# Export input data for model
# Means
c_d_before_CA   <- mean(df_no_ev_24$c_d_before_CA)
c_d_after_CA1   <- mean(df_no_ev_24$c_d_after_CA1 )
c_d_after_CA2   <- mean(df_no_ev_24$c_d_after_CA2 )
c_d_SAF         <- mean(df_no_ev_24$c_d_SAF      )
c_CA_DBC        <- mean(df_no_ev_24$v_pr_khm_p1)
c_AAD           <- mean(df_no_ev_24$v_farm_e_aad_p0)

# Standard deviations
c_sd_d_before_CA   <- sd(df_no_ev_24$c_d_before_CA)
c_sd_d_after_CA1   <- sd(df_no_ev_24$c_d_after_CA1)
c_sd_d_after_CA2   <- sd(df_no_ev_24$c_d_after_CA2)
c_sd_d_SAF         <- sd(df_no_ev_24$c_d_SAF      )

# Standard errors
c_se_d_before_CA   <- sd(df_no_ev_24$c_d_before_CA) / sqrt(length(df_no_ev_24$c_d_before_CA))
c_se_d_after_CA1   <- sd(df_no_ev_24$c_d_after_CA1) / sqrt(length(df_no_ev_24$c_d_after_CA1))
c_se_d_after_CA2   <- sd(df_no_ev_24$c_d_after_CA2) / sqrt(length(df_no_ev_24$c_d_after_CA2))
c_se_d_SAF         <- sd(df_no_ev_24$c_d_SAF      ) / sqrt(length(df_no_ev_24$c_d_SAF))

# Number of patients
n_c_Vektis         <- nrow(df_no_ev_24)

# Save in .RData
save(c_d_SAF, c_d_before_CA, c_d_after_CA1, c_d_after_CA2, c_CA_DBC, c_AAD,
    c_sd_d_SAF, c_sd_d_before_CA, c_sd_d_after_CA1, c_sd_d_after_CA2,
    c_CA_DBC_params, c_AAD_params, n_c_Vektis, file = here::here("input", "c_Vektis_update.RData"))






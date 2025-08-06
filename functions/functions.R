#### Rounding ####
# Function to round age upwards when >= .5 and downwards when < .5
round_age <- function(x) {
  floor(x + 0.5)
}

#### Conversion of proportion to and from logit scale ####
# Convert proportions and CIs to logit scale
logit <- function(p) log(p/(1-p))
# Convert back to original scale using inverse logit (expit)
expit <- function(x) exp(x)/(1+exp(x))

#### Convert probability of event from one time period to cycle length ####
convert_probability_to_cl <- function(prob, time, cl) {
  return(1 - (1 - prob)^(cl/time))
}

#### Adjustment for inflation using consumer price index ####
# Cost year = 2024, except for drug costs (2025)
# Correct for inflation
adjust_inflation <- function(amount, source_year, target_year = 2024) {
  # Define yearly inflation rates (as decimals)
  infl_rates <- c(
    "2017" = 0.014,
    "2018" = 0.016,
    "2019" = 0.012,
    "2020" = 0.025,
    "2021" = 0.118,
    "2022" = 0.030,
    "2023" = 0.024
  )
  
  # Get years to include in calculation
  years <- as.character(source_year:(target_year-1))
  
  # Calculate inflation factor and apply
  factor <- prod(1 + infl_rates[years])
  return(amount * factor)
}

#### Plot health state trace ####
plot_m_TR <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_s,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_s,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}

# Plot line trace ####
plot_trace_line <- function(m_L) 
{
  m_TR <- t(apply(m_L, 2, function(x) table(factor(x, levels = c(1,2,3,4,5,6, 555, 999), 
                                                   ordered = TRUE))))
  m_TR <- m_TR/n_i
  colnames(m_TR) <- v_names_lines
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")
  plot(0:n_cycles, m_TR[, 1], type = "l", main = "Line membership trace", 
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  for (n_states in 2:length(v_names_lines)) {
    lines(0:n_cycles, m_TR[, n_states], col = n_states)
  }
  legend("topright", v_names_lines, col = 1:length(v_names_lines), 
         lty = rep(1, length(v_names_lines)), bty = "n", cex = 0.65)
}


# #### Fit multiple functional forms to survival data ####
# # without generalized gamma density function because it is defined using an incomplete gamma in the denominator and when the values of the parameters 
# # approach the extreme values (e.g. close to 0) then the function fails"
# 
# fit.fun <- function(time, status, data = data , add = FALSE, extrapolate = FALSE, times)  
# {
#   #Extract the right data columns 
#   data$time   <-   data[,   time]  
#   data$status <-   data[, status]  
# 
#     if (extrapolate == TRUE)  {
#     plot.times <- max(times)
#   } else if  (extrapolate == FALSE) {
#     plot.times <- max(data$time)
#   }
#   
#   # Progression free survival  
#   KM.fit     <-     survfit(Surv(time, status) ~ 1, data = data)                         # fit Kaplan-Meier curve 
#   fit.llogis <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "llogis" )       # fit model with loglogistic distribution
#   fit.weib   <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "weibull")       # fit model with Weibull distribution
#   fit.lnorm  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "lnorm"  )       # fit model with lognormal distribution
#   fit.gamma  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gamma"  )       # fit model with gamma distribution 
#   fit.exp    <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "exp"    )       # fit model with exponential distribution
#   fit.gompertz  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gompertz"  ) # fit model with gompertz  
#   fit.rps    <- flexsurvspline(Surv(time, status) ~ 1, data = data, k = 2  ) # fit model with splines  )
#   
#   # extarapolate all models beyond the KM curve
#   if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
#   if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F, mark.time= T)}
#   lines(fit.llogis,   t = times, col = 2, ci = F)
#   lines(fit.weib,     t = times, col = 3, ci = F)
#   lines(fit.lnorm,    t = times, col = 4, ci = F)
#   lines(fit.gamma,    t = times, col = 5, ci = F)
#   lines(fit.exp,      t = times, col = 6, ci = F)
#   lines(fit.gompertz, t = times, col = 7, ci = F)
#   lines(fit.rps,      t = times, col = 8, ci = F)
#   legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gamma", "Exponential", "Gompertz", "rps"), col = 1:8, lty = rep(1, 9), bty="n")
#   
#   # compare AIC values
#   AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
#                Weibull     = AIC(fit.weib), 
#                Lognormal   = AIC(fit.lnorm), 
#                Gamma       = AIC(fit.gamma),
#                Exponential = AIC(fit.exp),
#              Gompertz = AIC(fit.gompertz),
#              Rps = AIC(fit.rps))
#               
#   AIC= round(AIC,3)
#   
#   # compare BIC values
#   BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
#                Weibull     = BIC(fit.weib), 
#                Lognormal   = BIC(fit.lnorm), 
#                Gamma       = BIC(fit.gamma),
#                Exponential = BIC(fit.exp),
#                Gompertz    = BIC(fit.gompertz),
#                Rps = BIC(fit.rps))
#   
#   BIC <- round(BIC,3)
#   
#   res <- list(Loglogistic = fit.llogis,
#               Weibull     = fit.weib,
#               Lognormal   = fit.lnorm, 
#               Gamma       = fit.gamma,
#               Exponential = fit.exp, 
#               Gompertz   = fit.gompertz,
#               Rps        = fit.rps,
#               AIC         = AIC,
#               BIC         = BIC)
#   res
# }
# 

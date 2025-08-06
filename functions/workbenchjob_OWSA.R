# 01 Settings and load required packages, functions and data
  
# 01.0 Clear memory

rm(list = ls())     # clear memory (removes all the variables from the workspace)
options(scipen=999) # disable scientific notation
set.seed(1)         # set the seed for random number generation


## 01.1 Load packages
if (!require('pacman')) install.packages('pacman'); library(pacman) # use this package to conveniently install other packages

# load (install if required) packages from CRAN
p_load("here", "dplyr", "ggplot2", "MASS", "tictoc", "tidyr", "survival", 
       "flexsurv", "survminer", "doParallel", "betareg", "truncnorm", "parallel", "foreach", "Matrix", "splines")

# load (install if required) packages from GitHub
# install_github("DARTH-git/dampack", force = TRUE) Uncomment if there is a newer version
# install_github("DARTH-git/darthtools", force = TRUE) Uncomment if there is a newer version
p_load_gh("DARTH-git/dampack", "DARTH-git/darthtools", "DARTH-git/dectree")


## 01.2 Settings
n_i                <- 50000                        # number of simulated individuals

# Model settings 
min_age            <- 18                           # min age: only include adults
max_age            <- 108                          # max age based on max in Human Mortality Database
cl                 <- 0.5                          # cycle length in years: should be 0.5 because not all 
# input parameters are automatically varied (i.e. TTE Vektis)
n_t                <- (max_age-min_age)/cl         # time horizon in cycles
n_cycles           <- n_t                          # number of cycles 
v_n                <- c("SFAF", "SAF", "D")        # model states names
v_names_states     <- v_n                          # variable needed for plot_trace_microsim
n_states           <- n_s <- length(v_n)           # the number of health states
wtp                <- 20000                        # willingness to pay threshold (based on proportional shortfall AF)
d_rc               <- 0.030                        # discount rate costs
d_re               <- 0.015                        # discount rate effects
d_rc               <- (1+d_rc)^cl-1                # discount rate costs to x-weekly 
d_re               <- (1+d_re)^cl-1                # discount rate effects to x-weekly 
v_dwc              <- 1 / ((1 + d_rc) ^ (0:n_t))   # per period discount weight costs
v_dwe              <- 1 / ((1 + d_re) ^ (0:n_t))   # per period discount weight effects
v_wcc              <- darthtools::gen_wcc(n_cycles = n_t, method = "Simpson1/3")  # within-cycle correction (WCC) 
v_names_lines      <- c("TRT1", "TRT2", "TRT3", "TRT4", "TRT5", "TRT6", "no-treatment", "D") # line names

# Base case or scenario
# The scenario "TTE_2CV" includes the time to event data from the health insurers declaration data (Vektis) of the scenario analysis where AF symptoms are defined as re-ablation, His-ablation, MAZE procedure, switch to amiodarone or at least 2 cardioversions within one year (i.e. single cardioversions without a second cardioversion within a year are not counted as recurrence of AF symptoms). In the base case all cardioversions are counted as recurrence of AF symptoms after which a treatment switch will be incurred in the model.

scenario <- "basecase" # "basecase" or "TTE_2CV"


## 01.3 Load functions
# Load generic functions
source(here::here("functions", "functions.R"))
source(here::here("functions", "model functions.R"))

# 0.8 OWSA
# 0.8.1 Define input parameters

## Define two sequences for which you want to run the OWSA
# CA-AAD-AAD-AAD-CA-CA
TRT1a <- "CA"
TRT2a <- "AAD"
TRT3a <- "AAD"
TRT4a <- "AAD"
TRT5a <- "CA"
TRT6a <- "CA"

# CA-AAD-AAD-AAD-AAD-CA
TRT1b <- "CA"
TRT2b <- "AAD"
TRT3b <- "AAD"
TRT4b <- "AAD"
TRT5b <- "AAD"
TRT6b <- "CA"

l_params_owsa <- l_params_owsa_min <- l_params_owsa_max <- l_params 

# Make a vector with names that you wish to appear in the tornado plot that are easier to understand for readers
#v_names_params <- c() 

# Proportion with which you vary the base case settings
lower <- 0.8
upper <- 1.2

# Patient characteristics: confidence intervals derived from Mol et al. 2022 using normal and beta distributions
l_params_owsa_min$age                      <- age
l_params_owsa_max$age                      <- age
l_params_owsa_min$sd_age                   <- sd_age
l_params_owsa_max$sd_age                   <- sd_age
owsa_min_age                               <- age - qnorm(0.975) * (sd_age/sqrt(n_total))
owsa_max_age                               <- age + qnorm(0.975) * (sd_age/sqrt(n_total))
l_params_owsa_min$v_age                    <- round(rtruncnorm(n_i, min_age, max_age, owsa_min_age, sd_age)) 
l_params_owsa_max$v_age                    <- round(rtruncnorm(n_i, min_age, max_age, owsa_max_age, sd_age)) 

l_params_owsa_min$prop_female              <- prop_female
l_params_owsa_max$prop_female              <- prop_female
owsa_min_prop_female                       <- qbeta(0.025, n_female, n_male)
owsa_max_prop_female                       <- qbeta(0.975, n_female, n_male)
l_params_owsa_min$v_Sex                    <- rbinom(n_i, 1, owsa_min_prop_female)
l_params_owsa_max$v_Sex                    <- rbinom(n_i, 1, owsa_max_prop_female)

l_params_owsa_min$prop_parox               <- qbeta(0.025, n_paroxysmal, n_persistent)
l_params_owsa_max$prop_parox               <- qbeta(0.975, n_paroxysmal, n_persistent)

# Clinical input parameters 
# Confidence intervals as reported by Cochrane Netherlands
l_params_owsa_min$RR_nai                   <- unname(exp(predrema1$ci.lb)) 
l_params_owsa_max$RR_nai                   <- unname(exp(predrema1$ci.ub))
l_params_owsa_min$RR_AAD_exp_parox         <- unname(exp(predrema2$ci.lb)) 
l_params_owsa_max$RR_AAD_exp_parox         <- unname(exp(predrema2$ci.ub))
l_params_owsa_min$RR_AAD_exp_pers          <- unname(exp(predrema3$ci.lb)) 
l_params_owsa_max$RR_AAD_exp_pers          <- unname(exp(predrema3$ci.ub))
l_params_owsa_min$RR_CA_exp                <- 1.61 # based on 1 study
l_params_owsa_max$RR_CA_exp                <- 2.80 # based on 1 study
l_params_owsa_min$prop_AF_CA_nai           <- unname(plogis(predrema1i$ci.lb))
l_params_owsa_max$prop_AF_CA_nai           <- unname(plogis(predrema1i$ci.ub)) 
l_params_owsa_min$prop_AF_CA_AAD_exp_parox <- unname(plogis(predrema2i$ci.lb))
l_params_owsa_max$prop_AF_CA_AAD_exp_parox <- unname(plogis(predrema2i$ci.ub))
l_params_owsa_min$prop_AF_CA_AAD_exp_pers  <- unname(plogis(predrema3i$ci.lb))
l_params_owsa_max$prop_AF_CA_AAD_exp_pers  <- unname(plogis(predrema3i$ci.ub))
l_params_owsa_min$p_disc_AAD_nai           <- 0.067
l_params_owsa_max$p_disc_AAD_nai           <- 0.206
l_params_owsa_min$p_disc_AAD_exp           <- 0.161
l_params_owsa_max$p_disc_AAD_exp           <- 0.305

# Confidence intervals derived from survival curve of Vektis data
l_params_owsa_min$prop_AF_CA_CA_exp_1      <- 0.2439314
l_params_owsa_max$prop_AF_CA_CA_exp_1      <- 0.2717558 
l_params_owsa_min$prop_AF_CA_CA_exp_2      <- 0.2843864
l_params_owsa_max$prop_AF_CA_CA_exp_2      <- 0.3713022  

# Confidence intervals of time to AF symptoms derived from survial curves of Vektis data
l_params_owsa_min$p_SFAF_CA_SAF_L1   <- trans_prob(df_s_ablation1[[1]]$lcl) 
l_params_owsa_max$p_SFAF_CA_SAF_L1   <- trans_prob(df_s_ablation1[[1]]$ucl) 
l_params_owsa_min$p_SFAF_CA_SAF_L2   <- trans_prob(df_s_ablation2[[1]]$lcl) 
l_params_owsa_max$p_SFAF_CA_SAF_L2   <- trans_prob(df_s_ablation2[[1]]$ucl) 
l_params_owsa_min$p_SFAF_CA_SAF_L3   <- trans_prob(df_s_ablation3[[1]]$lcl) 
l_params_owsa_max$p_SFAF_CA_SAF_L3   <- trans_prob(df_s_ablation3[[1]]$ucl) 

# Confidence intervals as reported by from Vinter et al.
l_params_owsa_min$HR_EM <- 1.71
l_params_owsa_max$HR_EM <- 2.36

# Confidence intervals derived from beta distribution of counts reported in ablation register of NHR in 2023
l_params_owsa_min$p_CT_CA      <- qbeta(0.025, n_CT_CA, n_tot_CT_CA-n_CT_CA)                 
l_params_owsa_max$p_CT_CA      <- qbeta(0.975, n_CT_CA, n_tot_CT_CA-n_CT_CA)  
l_params_owsa_min$p_PP_CA      <- qbeta(0.025, n_PP_CA, n_tot_PP_CA-n_PP_CA)                 
l_params_owsa_max$p_PP_CA      <- qbeta(0.975, n_PP_CA, n_tot_PP_CA-n_PP_CA)  
l_params_owsa_min$p_VC_CA      <- qbeta(0.025, n_VC_CA, n_tot_VC_CA-n_VC_CA)                 
l_params_owsa_max$p_VC_CA      <- qbeta(0.975, n_VC_CA, n_tot_VC_CA-n_VC_CA)  

# Vary coefficients productivity 
l_params_owsa_min$beta_reg_attendance_intercept <- beta_reg_attendance_intercept - qnorm(0.975) * 0.04763  
l_params_owsa_max$beta_reg_attendance_intercept <- beta_reg_attendance_intercept + qnorm(0.975) * 0.04763 
l_params_owsa_min$beta_reg_attendance_SAF       <- beta_reg_attendance_SAF - qnorm(0.975) * 0.08072    
l_params_owsa_max$beta_reg_attendance_SAF       <- beta_reg_attendance_SAF + qnorm(0.975) * 0.08072   

l_params_owsa_min$beta_reg_presenteeism_intercept <- beta_reg_presenteeism_intercept - qnorm(0.975) * 0.02500    
l_params_owsa_max$beta_reg_presenteeism_intercept <- beta_reg_presenteeism_intercept + qnorm(0.975) * 0.02500   
l_params_owsa_min$beta_reg_presenteeism_SAF       <- beta_reg_presenteeism_SAF - qnorm(0.975) * 0.04384      
l_params_owsa_max$beta_reg_presenteeism_SAF       <- beta_reg_presenteeism_SAF + qnorm(0.975) * 0.04384   

# Vary informal care coefficients based on SE in de Groot et al.
# Estimate coefficients of regression formula based on sampling from vcov and vary them separately
# Small adjustment to make the vcov "positive definite" to ensure it has certain mathematical properties that are required for multivariate normal sampling. The adjustments needed here were essentially zero (0.00000000000000000001%).
m_vcov_ic_log <- matrix(as.numeric(unlist(vcov_ic_log)), nrow = 5, ncol = 5)
vcov_ic_log_fixed <- as.matrix(nearPD(m_vcov_ic_log)$mat)

coef_samples_ic_log <- as.data.frame(mvrnorm(n = 1000, mu = params_ic_log, Sigma = vcov_ic_log_fixed))
params_ic_log_lower <- apply(coef_samples_ic_log, 2, function(x) {
  quantile(x, probs = 0.025)
})
params_ic_log_upper <- apply(coef_samples_ic_log, 2, function(x) {
  quantile(x, probs = 0.975)
})

coef_samples_ic_hours <- as.data.frame(mvrnorm(n = 1000, mu = params_ic_hours, Sigma = vcov_ic_hours[c(1:4), c(1:4)]))
params_ic_hours_lower <- apply(coef_samples_ic_hours, 2, function(x) {
  quantile(x, probs = 0.025)
})
params_ic_hours_upper <- apply(coef_samples_ic_hours, 2, function(x) {
  quantile(x, probs = 0.975)
})

l_params_owsa_min$params_ic_log_intercept <- params_ic_log_lower[1]
l_params_owsa_min$params_ic_log_female    <- params_ic_log_lower[2]
l_params_owsa_min$params_ic_log_age       <- params_ic_log_lower[3]
l_params_owsa_min$params_ic_log_age2      <- params_ic_log_lower[4]
l_params_owsa_min$params_ic_log_T2D       <- params_ic_log_lower[5]
l_params_owsa_max$params_ic_log_intercept <- params_ic_log_upper[1]
l_params_owsa_max$params_ic_log_female    <- params_ic_log_upper[2]
l_params_owsa_max$params_ic_log_age       <- params_ic_log_upper[3]
l_params_owsa_max$params_ic_log_age2      <- params_ic_log_upper[4]
l_params_owsa_max$params_ic_log_T2D       <- params_ic_log_upper[5]

l_params_owsa_min$params_ic_hours_intercept <- params_ic_hours_lower[1]
l_params_owsa_min$params_ic_hours_female    <- params_ic_hours_lower[2]
l_params_owsa_min$params_ic_hours_age       <- params_ic_hours_lower[3]
l_params_owsa_min$params_ic_hours_T2D       <- params_ic_hours_lower[4]
l_params_owsa_max$params_ic_hours_intercept <- params_ic_hours_upper[1]
l_params_owsa_max$params_ic_hours_female    <- params_ic_hours_upper[2]
l_params_owsa_max$params_ic_hours_age       <- params_ic_hours_upper[3]
l_params_owsa_max$params_ic_hours_T2D       <- params_ic_hours_upper[4]

# Confidence interval derived from standard deviation and sample size in van Dries et al. based on normal distribution
l_params_owsa_min$c_IC_77_AF_2yr            <- c_IC_77_AF_2yr - qnorm(0.975) * (321.76/sqrt(425))
l_params_owsa_max$c_IC_77_AF_2yr            <- c_IC_77_AF_2yr + qnorm(0.975) * (321.76/sqrt(425))

# Confidence interval derived from GLM model for coefficient of disutility of symptoms of AF
du_SAF_CI <- confint(glm_AVATARAF_NL, "sympAF_new", level = 0.95)
l_params_owsa_min$du_SAF <- abs(unname(du_SAF_CI[2]))
l_params_owsa_max$du_SAF <- abs(unname(du_SAF_CI[1]))

l_params_owsa_min$gen_pop_utility <- "lwr"
l_params_owsa_max$gen_pop_utility <- "upr"

# No confidence interval available, so varied +/- 20%
l_params_owsa_min$c_AAD    <- qlnorm(0.025, c_AAD_params$estimate[1], c_AAD_params$estimate[2])
l_params_owsa_max$c_AAD    <- qlnorm(0.975, c_AAD_params$estimate[1], c_AAD_params$estimate[2])
l_params_owsa_min$c_CA_DBC <- qweibull(0.025, shape = c_CA_DBC_params$estimate[1], scale = c_CA_DBC_params$estimate[2])
l_params_owsa_max$c_CA_DBC <- qweibull(0.975, shape = c_CA_DBC_params$estimate[1], scale = c_CA_DBC_params$estimate[2])

l_params_owsa_min$c_d_before_CA <- qnorm(0.025, c_d_before_CA, (c_sd_d_before_CA/sqrt(1495)))
l_params_owsa_max$c_d_before_CA <- qnorm(0.975, c_d_before_CA, (c_sd_d_before_CA/sqrt(1495)))
l_params_owsa_min$c_d_after_CA1 <- qnorm(0.025, c_d_after_CA1, (c_sd_d_after_CA1/sqrt(1495)))
l_params_owsa_max$c_d_after_CA1 <- qnorm(0.975, c_d_after_CA1, (c_sd_d_after_CA1/sqrt(1495)))
l_params_owsa_min$c_d_after_CA2 <- qnorm(0.025, c_d_after_CA2, (c_sd_d_after_CA2/sqrt(1495)))
l_params_owsa_max$c_d_after_CA2 <- qnorm(0.975, c_d_after_CA2, (c_sd_d_after_CA2/sqrt(1495)))
l_params_owsa_min$c_d_SAF       <- qnorm(0.025, c_d_SAF, (c_sd_d_SAF/sqrt(1495)))
l_params_owsa_max$c_d_SAF       <- qnorm(0.975, c_d_SAF, (c_sd_d_SAF/sqrt(1495)))

v_names_params_20perc <- c("prop_work_M", "prop_work_F", "v_hours_per_week_M",
                           "v_hours_per_week_F", "c_hourly_wage_2022", "c_IC_hr",
                           "LE_M_77", "LE_F_77", "du_CT_mo", "du_PP_mo", "du_VC_mo")
l_params_owsa_min[v_names_params_20perc] <- lapply(l_params[v_names_params_20perc], "*", lower) 
l_params_owsa_max[v_names_params_20perc] <- lapply(l_params[v_names_params_20perc], "*", upper) 

l_params_owsa_min$df_FMC[,3:4] <- l_params$df_FMC[,3:4]*lower
l_params_owsa_max$df_FMC[,3:4] <- l_params$df_FMC[,3:4]*upper

# Save base case values
l_params_owsa_bc  <- l_params

# Make a dataframe of the input parameters mean, min and max that are not strings, vectors of length >1 or dataframes or are not varied in OWSA
single_value_names <- v_names_params[!v_names_params %in% c("age", "prop_female", "sd_age", "v_age", "v_Sex", 
                                                            "p_SFAF_CA_SAF_L1", "p_SFAF_CA_SAF_L2", "p_SFAF_CA_SAF_L3", 
                                                             "df_FMC","gen_pop_utility",  "mod_splines_coef")]
mean_values <- min_values <- max_values <- NULL
for (i in 1:length(single_value_names)) {
  param_name <- single_value_names[i]
  mean_values[i] <- l_params[[param_name]]
  min_values[i] <- l_params_owsa_min[[param_name]]
  max_values[i] <- l_params_owsa_max[[param_name]]
}

df_owsa_params <- data.frame(
  parameter = single_value_names,
  mean = mean_values,
  min = min_values,
  max = max_values
)

# Add columns to check the relationships
df_owsa_params$min_less_than_mean <- df_owsa_params$min < df_owsa_params$mean
df_owsa_params$max_greater_than_mean <- df_owsa_params$max > df_owsa_params$mean

#View(df_owsa_params)

# Check values that are in a vector or dataframe
df_p_SFAF_CA_SAF_L1 <- data.frame(mean = l_params$p_SFAF_CA_SAF_L1,
                                  min = l_params_owsa_min$p_SFAF_CA_SAF_L1,
                                  max = l_params_owsa_max$p_SFAF_CA_SAF_L1)
df_p_SFAF_CA_SAF_L1$min_less_than_mean <- df_p_SFAF_CA_SAF_L1$min < df_p_SFAF_CA_SAF_L1$mean
df_p_SFAF_CA_SAF_L1$max_greater_than_mean <- df_p_SFAF_CA_SAF_L1$max > df_p_SFAF_CA_SAF_L1$mean
#View(df_p_SFAF_CA_SAF_L1)

df_p_SFAF_CA_SAF_L2 <- data.frame(mean = l_params$p_SFAF_CA_SAF_L2,
                                  min = l_params_owsa_min$p_SFAF_CA_SAF_L2,
                                  max = l_params_owsa_max$p_SFAF_CA_SAF_L2)
df_p_SFAF_CA_SAF_L2$min_less_than_mean <- df_p_SFAF_CA_SAF_L2$min < df_p_SFAF_CA_SAF_L2$mean
df_p_SFAF_CA_SAF_L2$max_greater_than_mean <- df_p_SFAF_CA_SAF_L2$max > df_p_SFAF_CA_SAF_L2$mean
#View(df_p_SFAF_CA_SAF_L2)

df_p_SFAF_CA_SAF_L3 <- data.frame(mean = l_params$p_SFAF_CA_SAF_L3,
                                  min = l_params_owsa_min$p_SFAF_CA_SAF_L3,
                                  max = l_params_owsa_max$p_SFAF_CA_SAF_L3)
df_p_SFAF_CA_SAF_L3$min_less_than_mean <- df_p_SFAF_CA_SAF_L3$min < df_p_SFAF_CA_SAF_L3$mean
df_p_SFAF_CA_SAF_L3$max_greater_than_mean <- df_p_SFAF_CA_SAF_L3$max > df_p_SFAF_CA_SAF_L3$mean
#View(df_p_SFAF_CA_SAF_L3)

df_FMC_other_y <- data.frame(mean = l_params$df_FMC$other_y,
                            min = l_params_owsa_min$df_FMC$other_y,
                            max = l_params_owsa_max$df_FMC$other_y)
df_FMC_other_y$min_less_than_mean <- df_FMC_other_y$min < df_FMC_other_y$mean
df_FMC_other_y$max_greater_than_mean <- df_FMC_other_y$max > df_FMC_other_y$mean
#View(df_FMC_other_y)

df_FMC_last_y <- data.frame(mean = l_params$df_FMC$last_y,
                            min = l_params_owsa_min$df_FMC$last_y,
                            max = l_params_owsa_max$df_FMC$last_y)
df_FMC_last_y$min_less_than_mean <- df_FMC_last_y$min < df_FMC_last_y$mean
df_FMC_last_y$max_greater_than_mean <- df_FMC_last_y$max > df_FMC_last_y$mean
#View(df_FMC_last_y)

df_c_states_min <- l_params_owsa_max$df_c_states[, 2:4] > l_params$df_c_states[, 2:4]
df_c_states_max <- l_params_owsa_min$df_c_states[, 2:4] < l_params$df_c_states[, 2:4]
# df_c_states_min
# df_c_states_max

## 08.2 Run the OWSA
# Run the base-case analysis
# 
PSA <- T #if FALSE increases number of stored parameters, is slower but gives additional insights in deterministic setting
mainresults <- T # if TRUE only save te_hat and tc_hat
base_case_result_a     <- MicroSim(l_params, n_i, df_X, TRT1 = TRT1a, TRT2 = TRT2a, TRT3 = TRT3a,
                                   TRT4 = TRT4a, TRT5 = TRT5a, TRT6 = TRT6a, seed = 1)
base_case_result_a$NHB <-  base_case_result_a$te_hat - (base_case_result_a$tc_hat/wtp)

base_case_result_b     <- MicroSim(l_params, n_i, df_X, TRT1 = TRT1b, TRT2 = TRT2b, TRT3 = TRT3b,
                                                           TRT4 = TRT4b, TRT5 = TRT5b, TRT6 = TRT6a, seed = 1)
base_case_result_b$NHB <-  base_case_result_b$te_hat - (base_case_result_b$tc_hat/wtp)

d_e <- base_case_result_a$te_hat-base_case_result_b$te_hat
d_c <- base_case_result_a$tc_hat-base_case_result_b$tc_hat
iNHB <- d_e - (d_c/wtp)

# Prepare list object for lapply to use
l_owsa <- list()
l_owsa <- rep(list(l_params),(length(v_names_params)*2))

for (i in 1:length(v_names_params)){
  l_owsa[[i]][[i]]                        <- l_params_owsa_min[[i]]
  l_owsa[[i+length(v_names_params)]][[i]] <- l_params_owsa_max[[i]]  
}

# Setup parallel processing with optimized memory management
num_cores <- min(detectCores() - 1, 4)  # Limit cores to prevent memory overload
cluster <- makeCluster(num_cores, type = "PSOCK")
registerDoParallel(cluster)

# Set memory limits and clean environment
clusterEvalQ(cluster, {
  options(future.globals.maxSize = 4 * 1024^3)  # 4GB per worker
  gc(reset = TRUE)
  library(dplyr) 
  library(darthtools)
  library(splines)
  source(here::here("functions", "functions.R"))
})

clusterExport(cluster, c('MicroSim', 'Create_df_X', 'n_i', 'rtruncnorm', 'v_n', 'n_states',
                         'n_t',  'Costs', 'Effs', 'Probs', 'PrepareProbs', 'PrepareCosts', 
                         'cl', 'v_dwc', 'v_dwe', 'v_wcc', 'min_age', 'max_age', 'PSA', 'mainresults',
                         'TRT1a','TRT2a','TRT3a','TRT4a', 'TRT5a', 'TRT6a', 
                         'TRT1b','TRT2b','TRT3b','TRT4b', 'TRT5b', 'TRT6b', 'df_mort',
                         'beta_reg_attendance', 'beta_reg_presenteeism', 'params_ic_log', 'params_ic_hours','mod_splines'))

# Process parameters in chunks
chunk_size <- ceiling(length(l_owsa) / (3 * num_cores))
l_owsa_chunks <- split(l_owsa, ceiling(seq_along(l_owsa)/chunk_size))

# Initialize results containers
l_results_owsa_a <- vector("list", length(l_owsa))
l_results_owsa_b <- vector("list", length(l_owsa))

tic()
# Process each chunk with clean-up between
for (i in seq_along(l_owsa_chunks)) {
  cat("Processing chunk", i, "of", length(l_owsa_chunks), "\n")
  # Treatment A chunk processing
  chunk_results_a <- parLapply(cluster, l_owsa_chunks[[i]], function(o){
    # Return minimal required data
    result <- MicroSim(o, n_i, df_X, 
                       TRT1 = TRT1a, TRT2 = TRT2a, TRT3 = TRT3a, 
                       TRT4 = TRT4a, TRT5 = TRT5a, TRT6 = TRT6a, 
                       seed = 1)
    list(tc_hat = result$tc_hat, te_hat = result$te_hat)
  })
  # Treatment B chunk processing
  chunk_results_b <- parLapply(cluster, l_owsa_chunks[[i]], function(o){
    result <- MicroSim(o, n_i, df_X, 
                       TRT1 = TRT1b, TRT2 = TRT2b, TRT3 = TRT3b, 
                       TRT4 = TRT4b, TRT5 = TRT5b, TRT6 = TRT6b, 
                       seed = 1)
    list(tc_hat = result$tc_hat, te_hat = result$te_hat)
  })
  # Store chunk results
  chunk_indices <- ((i-1) * chunk_size + 1):min(i * chunk_size, length(l_owsa))
  for (j in seq_along(chunk_indices)) {
    l_results_owsa_a[[chunk_indices[j]]] <- chunk_results_a[[j]]
    l_results_owsa_b[[chunk_indices[j]]] <- chunk_results_b[[j]]
  }
  # Free memory
  rm(chunk_results_a, chunk_results_b)
  gc(full = TRUE)
  # Save progress
  if (i %% 2 == 0 || i == length(l_owsa_chunks)) {
    saveRDS(list(a = l_results_owsa_a, b = l_results_owsa_b), 
            file = here::here("output", paste0("owsa_progress_", i, ".rds")))
  }
}
toc()

# Free cluster resources
stopCluster(cluster)

# Efficient vectorized results processing 
process_results <- function(results_list, param_names, param_count, wtp) {
  # Split results for min and max values
  results_min <- results_list[1:param_count]
  results_max <- results_list[(param_count+1):(2*param_count)]
  # Extract values using vectorized operations
  results_df <- data.frame(
    Parameter = param_names,
    Lower_Bound_c = sapply(results_min, function(x) x$tc_hat),
    Lower_Bound_e = sapply(results_min, function(x) x$te_hat),
    Upper_Bound_c = sapply(results_max, function(x) x$tc_hat),
    Upper_Bound_e = sapply(results_max, function(x) x$te_hat)
  )
  # Vectorized calculations
  results_df$Lower_Bound_NHB <- results_df$Lower_Bound_e - (results_df$Lower_Bound_c/wtp)
  results_df$Upper_Bound_NHB <- results_df$Upper_Bound_e - (results_df$Upper_Bound_c/wtp)
  results_df$UL_Difference_c <- abs(results_df$Upper_Bound_c - results_df$Lower_Bound_c)
  results_df$UL_Difference_e <- abs(results_df$Upper_Bound_e - results_df$Lower_Bound_e)
  results_df$UL_Difference_NHB <- abs(results_df$Upper_Bound_NHB - results_df$Lower_Bound_NHB)
  return(results_df)
}

# Process and save results for both treatments
results_owsa_a <- process_results(l_results_owsa_a, v_names_params, length(v_names_params), wtp)
write.csv(results_owsa_a, here::here("output", "results_owsa_a_rank1vs2.csv"))

results_owsa_b <- process_results(l_results_owsa_b, v_names_params, length(v_names_params), wtp)
write.csv(results_owsa_b, here::here("output", "results_owsa_b_rank1vs2.csv"))

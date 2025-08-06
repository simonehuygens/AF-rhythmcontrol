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

# Model settings 
min_age            <- 18                           # min age: only include adults
max_age            <- 108                          # max age based on max in Human Mortality Database
cl                 <- 0.5                          # cycle length in years: should be 0.5 because not all 
# input parameters are automatically varied (i.e. TTE Vektis)
n_t                <- (max_age-min_age)/cl         # time horizon in cycles
n_cycles           <- n_t                          # number of cycles 
n_i                <- 20000                        # number of simulated individuals
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

# 0.9.7 Test model stability in deterministic setting
#In this section we loop over a change in the sample size and then check the stability of the calculated ICER


# Function to process a single iteration for a given sample size
process_iteration <- function(i, n_i, l_params, df_X) {
  outcomes_trt1 <- MicroSim(l_params, n_i, df_X, 
                            TRT1 = "AAD", TRT2 = "CA", TRT3 = "AAD", 
                            TRT4 = "AAD", TRT5 = "AAD", TRT6 = "AAD", seed = i)
  
  outcomes_trt2 <- MicroSim(l_params, n_i, df_X, 
                            TRT1 = "AAD", TRT2 = "AAD", TRT3 = "CA", 
                            TRT4 = "AAD", TRT5 = "AAD", TRT6 = "AAD", seed = i)
  
  d_e <- outcomes_trt2$te_hat-outcomes_trt1$te_hat
  d_c <- outcomes_trt2$tc_hat-outcomes_trt1$tc_hat
  iNHB <- d_e - (d_c/wtp)
  
  return(iNHB)
}

# Let op: duurt erg lang (>15 uur) door de 120000 en 150000 aan het eind!
# Main analysis with parallel processing
PSA                <- TRUE  # if TRUE stores less variables for 'faster' computing
iterations         <- 20
n_i_options        <- c(1000, 5000, 10000, 15000, 20000, 30000, 50000, 100000, 120000)
df_iNHB            <- matrix(nrow=iterations, ncol=length(n_i_options))
colnames(df_iNHB)  <- n_i_options

# Setup parallel processing
num_cores          <- detectCores() - 1  # Leave one core free for system processes
cluster            <- makeCluster(num_cores)
registerDoParallel(cluster)

# Export necessary objects to the cluster
clusterExport(cluster, c('MicroSim', 'Create_df_X', 'n_i', 'rtruncnorm', 'v_n', 'n_states',
                         'n_t',  'Costs', 'Effs', 'Probs', 'PrepareProbs', 'PrepareCosts', 
                         'cl', 'v_dwc', 'v_dwe', 'v_wcc', 'min_age', 'max_age', 'PSA',
                         'df_mort', 'beta_reg_attendance', 'beta_reg_presenteeism', 'params_ic_log', 'params_ic_hours', 'mod_splines'))

clusterEvalQ(cluster, 
             {library(dplyr) 
               library(darthtools)
               library(splines)
               source(here::here("functions", "functions.R"))
             })

# Loop over sample sizes
for (j in 1:length(n_i_options)) {
  n_i <- n_i_options[j]
  
  set.seed(2)
  v_age                 <- round(rtruncnorm(n_i, min_age, max_age, age, sd_age)) 
  v_Sex                 <- rbinom(n_i, 1, prop_female) 
  l_params$v_age        <- round(rtruncnorm(n_i, min_age, max_age, age, sd_age))
  l_params$v_Sex         <- rbinom(n_i, 1, prop_female)
  
  # Parallel processing of iterations
  results <- foreach(i = 1:iterations, 
                     .combine = 'c',
                     .packages = c("dplyr", "darthtools")) %dopar% {
                       process_iteration(i, n_i, l_params, df_X)
                     }
  
  df_iNHB[,j] <- results
  cat(paste0(round(j / length(n_i_options) * 100), '% completed\n'))
}

# Stop the cluster
stopCluster(cluster)

save(df_iNHB, file = here::here("output", "df_iNHB.RData"))



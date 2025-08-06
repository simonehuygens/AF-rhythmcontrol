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
n_i                <- 50000                        # number of simulated individuals
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

#### Treatment sequences ####
# The model has a maximum of 6 treatment lines
# Not all patients that enter the model will receive all 6 treatment lines, only those with AF recurrence move to the next line
df_trtseq <- expand.grid(line1 = c("AAD", "CA"), 
                         line2 = c("AAD", "CA"),
                         line3 = c("AAD", "CA"),
                         line4 = c("AAD", "CA"),
                         line5 = c("AAD", "CA"),
                         line6 = c("AAD", "CA"), stringsAsFactors = FALSE) 

#write.csv(df_trtseq, here::here("input", "df_trtseq.csv"))

# Create variable with names of the treatment sequences for results output
v_names_str <- c(paste(df_trtseq$line1, df_trtseq$line2, df_trtseq$line3, 
                       df_trtseq$line4, df_trtseq$line5, df_trtseq$line6, sep="-"))
n_str       <- length(v_names_str)      



## 01.3 Load functions
# Load generic functions
source(here::here("functions", "functions.R"))
source(here::here("functions", "model functions.R"))


# Mutiple sequences
# By specifying all the arguments in the MicroSim() function the simulation can be started
PSA <- F
mainresults <- F

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

tic()
results <- foreach(g = 1:nrow(df_trtseq), 
                   .packages = c("darthtools")) %dopar% { # darthtools voor sample v
                     MicroSim(l_params, n_i, df_X, TRT1 = df_trtseq[g,1], TRT2 = df_trtseq[g,2], 
                              TRT3 = df_trtseq[g,3], TRT4 = df_trtseq[g,4], TRT5 = df_trtseq[g,5], 
                              TRT6 = df_trtseq[g,6], seed = 1)
                   }
toc()
stopCluster(cluster)

# Store the mean costs, QALYs, and clinical outcomes of each strategy in new variables
v_C <- v_E <- tLY_hat <- tLY_undisc_hat <- p_L1 <- p_L2 <- p_L3 <- p_L4 <- p_L5 <- p_L6 <- p_noTrt <- 
  t_L1 <- t_L2 <- t_L3 <- t_L4 <- t_L5 <- t_L6 <- t_noTrt  <- NULL

for (i in 1:length(results)){
  v_C[i]               <- results[[i]]$tc_hat
  v_E[i]               <- results[[i]]$te_hat
  tLY_hat[i]           <- results[[i]]$tLY_hat
  tLY_undisc_hat[i]    <- results[[i]]$tLY_undisc_hat
  t_L1[i]              <- results[[i]]$t_L1   
  t_L2[i]              <- results[[i]]$t_L2   
  t_L3[i]              <- results[[i]]$t_L3   
  t_L4[i]              <- results[[i]]$t_L4  
  t_L5[i]              <- results[[i]]$t_L5  
  t_L6[i]              <- results[[i]]$t_L6  
  t_noTrt[i]           <- results[[i]]$t_noTrt
  p_L1[i]              <- results[[i]]$p_L1   
  p_L2[i]              <- results[[i]]$p_L2   
  p_L3[i]              <- results[[i]]$p_L3  
  p_L4[i]              <- results[[i]]$p_L4   
  p_L5[i]              <- results[[i]]$p_L5   
  p_L6[i]              <- results[[i]]$p_L6  
  p_noTrt[i]           <- results[[i]]$p_noTrt
}

# Combine the results in a data frame (without calculate_icers)
results_ce <- data.frame(Treatment_sequence = v_names_str,
                         Cost               = v_C,
                         QALYs              = v_E,
                         LYs                = tLY_hat,
                         LYs_undiscounted   = tLY_undisc_hat,
                         NHB                = v_E-(v_C/wtp),
                         Time_in_line1      = t_L1,
                         Time_in_line2      = t_L2,
                         Time_in_line3      = t_L3,
                         Time_in_line4      = t_L4,
                         Time_in_line5      = t_L5,
                         Time_in_line6      = t_L6,
                         Time_discontinued  = t_noTrt,
                         Prop_in_line1      = p_L1,
                         Prop_in_line2      = p_L2, 
                         Prop_in_line3      = p_L3,
                         Prop_in_line4      = p_L4,
                         Prop_in_line5      = p_L5, 
                         Prop_in_line6      = p_L6, 
                         Prop_discontinued  = p_noTrt
)

results_ce <- cbind(results_ce[, 1], separate(results_ce, col = "Treatment_sequence", c("L1","L2","L3", "L4", "L5", "L6"), "-"))
colnames(results_ce)[colnames(results_ce) == 'results_ce[, 1]'] <- 'Treatment_sequence'

# Save results
write.csv(results_ce, here::here("output", "results_bc_ni50000.csv"))
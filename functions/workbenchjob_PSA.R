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

# PSA settings
n_sim <- 1000  # Number of simulations 
n_i   <- 5000  # Number of patients 

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



# 0.9 PSA
# 0.9.1 Define PSA sequence
# 0.9.2 Define input values


# generate input data
set.seed(071818)    

gen_psa <- function(n_sim = 100, seed = 071818){
  with((l_params), {
    set.seed(seed) # set a seed to be able to reproduce the same results
    
    # Informal care & productivity: estimate coefficients here using the vcov
    # Include the new coefficient in separate vectors per variable below
    params_ic_log_psa            <- mvrnorm(n_sim, params_ic_log, nearPD(as.matrix(vcov_ic_log))$mat) 
    params_ic_hours_psa          <- mvrnorm(n_sim, params_ic_hours, vcov_ic_hours[1:4,1:4]) 
    beta_reg_attendance_psa      <- mvrnorm(n_sim, beta_reg_attendance$coefficients$mean, vcov_beta_reg_attendance[-3,-3])
    beta_reg_presenteeism_psa    <- mvrnorm(n_sim, beta_reg_presenteeism$coefficients$mean, vcov_beta_reg_presenteeism[-3,-3])
    
    df_psa <- data.frame(
      prop_parox                          = rbeta(n_sim, n_paroxysmal, n_persistent), 
      RR_nai                              = rlnorm(n_sim, predrema1$pred, predrema1$se),
      RR_AAD_exp_parox                    = rlnorm(n_sim, predrema2$pred, predrema2$se),
      RR_AAD_exp_pers                     = rlnorm(n_sim, predrema3$pred, predrema3$se),
      RR_CA_exp                           = rlnorm(n_sim, log(RR_CA_exp), (log(2.80)-log(1.61))/(2*1.96)),
      prop_AF_CA_nai                      = plogis(rlogis(n_sim, predrema1i$pred, predrema1i$se)), 
      prop_AF_CA_AAD_exp_parox            = plogis(rlogis(n_sim, predrema2i$pred, predrema2i$se)),
      prop_AF_CA_AAD_exp_pers             = plogis(rlogis(n_sim, predrema3i$pred, predrema3i$se)),
      prop_AF_CA_CA_exp_1                 = rnorm(n_sim, prop_AF_CA_CA_exp_1, (0.2717558-0.2439314) / (2*1.96)), 
      prop_AF_CA_CA_exp_2                 = rnorm(n_sim, prop_AF_CA_CA_exp_2, (0.3713022-0.2843864) / (2*1.96)),  
      HR_EM                               = exp(rnorm(n_sim, log(HR_EM), (log(2.36)-log(1.71))/(2*1.96))),
      p_disc_AAD_nai                      = rbeta(n_sim, 12, 99-12),
      p_disc_AAD_exp                      = plogis(rlogis(n_sim, pred_disc$pred, pred_disc$se)), 
      p_CT_CA                             = rbeta(n_sim, n_CT_CA, n_tot_CT_CA-n_CT_CA),                 
      p_PP_CA                             = rbeta(n_sim, n_PP_CA, n_tot_PP_CA-n_PP_CA),                 
      p_VC_CA                             = rbeta(n_sim, n_VC_CA, n_tot_VC_CA-n_VC_CA), 
      c_AAD                               = rlnorm(n_sim, c_AAD_params$estimate[1], c_AAD_params$estimate[2]), 
      c_CA_DBC                            = rweibull(n_sim, shape = c_CA_DBC_params$estimate[1],               
                                                     scale = c_CA_DBC_params$estimate[2]),  
      c_d_before_CA                       = rnorm(n_sim, c_d_before_CA, (c_sd_d_before_CA/sqrt(1495))), 
      c_d_after_CA1                       = rnorm(n_sim, c_d_after_CA1, (c_sd_d_after_CA1/sqrt(1495))), 
      c_d_after_CA2                       = rnorm(n_sim, c_d_after_CA2, (c_sd_d_after_CA2/sqrt(1495))), 
      c_d_SAF                             = rnorm(n_sim, c_d_SAF, (c_sd_d_SAF/sqrt(1495))),
      prop_work_M                         = prop_work_M        , # not varied in PSA, based on CBS data in whole Dutch population
      prop_work_F                         = prop_work_F        , # not varied in PSA, based on CBS data in whole Dutch population
      v_hours_per_week_M                  = v_hours_per_week_M , # not varied in PSA, based on CBS data in whole Dutch population
      v_hours_per_week_F                  = v_hours_per_week_F , # not varied in PSA, based on CBS data in whole Dutch population
      c_hourly_wage_2022                  = c_hourly_wage_2022 , # not varied in PSA, based on CBS data in whole Dutch population
      params_ic_log_intercept             = params_ic_log_psa[,1],
      params_ic_log_female                = params_ic_log_psa[,2],
      params_ic_log_age                   = params_ic_log_psa[,3],
      params_ic_log_age2                  = params_ic_log_psa[,4],
      params_ic_log_T2D                   = params_ic_log_psa[,5],
      params_ic_hours_intercept           = params_ic_hours_psa[,1],
      params_ic_hours_female              = params_ic_hours_psa[,2], 
      params_ic_hours_age                 = params_ic_hours_psa[,3],
      params_ic_hours_T2D                 = params_ic_hours_psa[,4],
      beta_reg_attendance_intercept       = beta_reg_attendance_psa[,1],
      beta_reg_attendance_SAF             = beta_reg_attendance_psa[,2],
      beta_reg_presenteeism_intercept     = beta_reg_presenteeism_psa[,1],
      beta_reg_presenteeism_SAF           = beta_reg_presenteeism_psa[,2],
      LE_M_77                             = LE_M_77        , # not varied in PSA, based on CBS data in whole Dutch population
      LE_F_77                             = LE_F_77        , # not varied in PSA, based on CBS data in whole Dutch population
      c_IC_hr                             = c_IC_hr        , # not varied in PSA
      c_IC_77_AF_2yr                      = rgamma(n_sim, shape = c_IC_77_AF_2yr_gamma$shape, scale = c_IC_77_AF_2yr_gamma$scale), 
      du_SAF                              = -mvrnorm(n_sim,glm_AVATARAF_NL$coefficients, vcov_glm_AVATARAF_NL)[,2],
      du_CT_mo                            = du_CT_mo,
      du_PP_mo                            = du_PP_mo,
      du_VC_mo                            = du_VC_mo
    )
    
    return(df_psa)
  }
  )
}

# # test run the PSA function
# gen_psa(10) 

# Generate PSA input dataset
df_psa_input <- gen_psa(n_sim = n_sim)

# # First six observations
# head(df_psa_input)

### Create other PSA input not within dataframe ###
set.seed(071818)  

#baseline population
patient_pop <- list()
for(i in 1:n_sim) {
  patient_pop[[i]] <- list(
    set.seed = i,
    v_age = round(rtruncnorm(n_i, min_age, max_age, age, sd_age)),
    v_Sex = rbinom(n_i, 1, prop_female)
  )
}

# Baseline utilities
mod_splines_psa                  <- mvrnorm(n_sim, coef(mod_splines), vcov_splines) 
apply(mod_splines_psa, 1, min)
summary(apply(mod_splines_psa, 1, max))

# Time to event probabilities
times  <- seq(0, 100, 0.5)
df_s_ablation1_PSA <- df_s_ablation2_PSA <- df_s_ablation3_PSA <- data.frame(times = times) # prepare dataframe survival
p_SFAF_CA_SAF_L1_PSA <- p_SFAF_CA_SAF_L2_PSA <- p_SFAF_CA_SAF_L3_PSA <- data.frame(times = times) # prepare dataframe transition prob

# Copy base case survival function, these will be overwritten in de PSA
best_fit_ablation1_PSA <- best_fit_ablation1  
best_fit_ablation2_PSA <- best_fit_ablation2 
best_fit_ablation3_PSA <- best_fit_ablation3 

# Sample parameters of the survival objects for the PSA
best_fit_ablation1_sample <- normboot.flexsurvreg(best_fit_ablation1, B = n_sim) # sample parameters     
best_fit_ablation2_sample <- normboot.flexsurvreg(best_fit_ablation2, B = n_sim) #sample parameters
best_fit_ablation3_sample <- normboot.flexsurvreg(best_fit_ablation3, B = n_sim) #sample parameters

for(i in 1:n_sim){
  best_fit_ablation1_PSA$res.t[,1] <- best_fit_ablation1_sample[i,] 
  best_fit_ablation2_PSA$res.t[,1] <- best_fit_ablation2_sample[i,] 
  best_fit_ablation3_PSA$res.t[,1] <- best_fit_ablation3_sample[i,] 
  
  df_s_ablation1_PSA[,i+1]        <- summary(best_fit_ablation1_PSA, t = times)[[1]][2] # i + 1 because first column is time
  p_SFAF_CA_SAF_L1_PSA[,i+1]      <- c(NA,trans_prob(df_s_ablation1_PSA[,i+1])) # i + 1 because first column is time
  
  df_s_ablation2_PSA[,i+1]        <- summary(best_fit_ablation2_PSA, t = times)[[1]][2] 
  p_SFAF_CA_SAF_L2_PSA[,i+1]      <- c(NA,trans_prob(df_s_ablation2_PSA[,i+1]))
  
  df_s_ablation3_PSA[,i+1]        <- summary(best_fit_ablation3_PSA, t = times)[[1]][2] 
  p_SFAF_CA_SAF_L3_PSA[,i+1]      <- c(NA,trans_prob(df_s_ablation3_PSA[,i+1]))
  
}

p_SFAF_CA_SAF_L1_PSA                    <- p_SFAF_CA_SAF_L1_PSA[,-1]
p_SFAF_CA_SAF_L2_PSA                    <- p_SFAF_CA_SAF_L2_PSA[,-1]
p_SFAF_CA_SAF_L3_PSA                    <- p_SFAF_CA_SAF_L3_PSA[,-1]

p_SFAF_CA_SAF_L1_PSA                    <- p_SFAF_CA_SAF_L1_PSA[-1,]
p_SFAF_CA_SAF_L2_PSA                    <- p_SFAF_CA_SAF_L2_PSA[-1,]
p_SFAF_CA_SAF_L3_PSA                    <- p_SFAF_CA_SAF_L3_PSA[-1,]

# summary(apply(p_SFAF_CA_SAF_L1_PSA, 1, min))
# summary(apply(p_SFAF_CA_SAF_L2_PSA, 1, min))
# summary(apply(p_SFAF_CA_SAF_L3_PSA, 1, min))
# 
# summary(apply(p_SFAF_CA_SAF_L1_PSA, 1, max))
# summary(apply(p_SFAF_CA_SAF_L2_PSA, 1, max))
# summary(apply(p_SFAF_CA_SAF_L3_PSA, 1, max))

# Initialize dataframes with PSA output
# Data frame of costs
df_c             <- as.data.frame(matrix(0, nrow = n_sim,  ncol = n_str))
colnames(df_c)   <- v_names_str

# Data frame of effectiveness
df_e             <- as.data.frame(matrix(0, nrow = n_sim, ncol = n_str))
colnames(df_e)   <- v_names_str

# 0.9.3 Run the PSA
PSA <- T
mainresults <- F


# Setup parallel processing
num_cores          <- detectCores()   # Leave one core free for system processes
cluster            <- makeCluster(num_cores)
registerDoParallel(cluster)

clusterExport(cluster, c('MicroSim', 'Create_df_X', 'n_i', 'rtruncnorm', 'v_n', 'n_states',
                         'n_t',  'Costs', 'Effs', 'Probs', 'PrepareProbs', 'PrepareCosts', 
                         'cl', 'v_dwc', 'v_dwe', 'v_wcc', 'min_age', 'max_age', 'PSA',
                         'gen_pop_utility', 'mod_splines'))

clusterEvalQ(cluster, 
             {library(dplyr) 
               library(darthtools)
               library(splines)
               source(here::here("functions", "functions.R"))
             })

iterations <- nrow(df_trtseq)
pb         <- txtProgressBar(max = iterations, style = 3)
progress   <- function(n) setTxtProgressBar(pb, n)
opts       <- list(progress = progress)
pb = txtProgressBar(min = 0, max = n_sim, initial = 0)

tic()
for(i in 1:n_sim){
  
  l_params                    <- list(
    prop_parox                         = df_psa_input$prop_parox[i]                         ,  
    RR_nai                             = df_psa_input$RR_nai[i]                             ,  
    RR_AAD_exp_parox                   = df_psa_input$RR_AAD_exp_parox[i]                   ,  
    RR_AAD_exp_pers                    = df_psa_input$RR_AAD_exp_pers[i]                    ,  
    RR_CA_exp                          = df_psa_input$RR_CA_exp[i]                          ,  
    prop_AF_CA_nai                     = df_psa_input$prop_AF_CA_nai[i]                     ,  
    prop_AF_CA_AAD_exp_parox           = df_psa_input$prop_AF_CA_AAD_exp_parox[i]           ,  
    prop_AF_CA_AAD_exp_pers            = df_psa_input$prop_AF_CA_AAD_exp_pers[i]            ,  
    prop_AF_CA_CA_exp_1                = df_psa_input$prop_AF_CA_CA_exp_1[i]                ,  
    prop_AF_CA_CA_exp_2                = df_psa_input$prop_AF_CA_CA_exp_2[i]                ,  
    HR_EM                              = df_psa_input$HR_EM[i]                              ,  
    p_disc_AAD_nai                     = df_psa_input$p_disc_AAD_nai[i]                     ,  
    p_disc_AAD_exp                     = df_psa_input$p_disc_AAD_exp[i]                     ,  
    p_CT_CA                            = df_psa_input$p_CT_CA[i]                            ,  
    p_PP_CA                            = df_psa_input$p_PP_CA[i]                            ,  
    p_VC_CA                            = df_psa_input$p_VC_CA[i]                            ,  
    c_AAD                              = df_psa_input$c_AAD[i]                              ,  
    c_CA_DBC                           = df_psa_input$c_CA_DBC[i]                           ,  
    c_d_before_CA                      = df_psa_input$c_d_before_CA[i]                      ,  
    c_d_after_CA1                      = df_psa_input$c_d_after_CA1[i]                      ,  
    c_d_after_CA2                      = df_psa_input$c_d_after_CA2[i]                      ,  
    c_d_SAF                            = df_psa_input$c_d_SAF[i]                            ,  
    prop_work_M                        = df_psa_input$prop_work_M[i]                        ,  
    prop_work_F                        = df_psa_input$prop_work_F[i]                        ,  
    v_hours_per_week_M                 = df_psa_input$v_hours_per_week_M[i]                 ,  
    v_hours_per_week_F                 = df_psa_input$v_hours_per_week_F[i]                 ,  
    c_hourly_wage_2022                 = df_psa_input$c_hourly_wage_2022[i]                 ,  
    params_ic_log_intercept            = df_psa_input$params_ic_log_intercept[i]            ,  
    params_ic_log_female               = df_psa_input$params_ic_log_female[i]               ,  
    params_ic_log_age                  = df_psa_input$params_ic_log_age[i]                  ,  
    params_ic_log_age2                 = df_psa_input$params_ic_log_age2[i]                 ,  
    params_ic_log_T2D                  = df_psa_input$params_ic_log_T2D[i]                  ,  
    params_ic_hours_intercept          = df_psa_input$params_ic_hours_intercept[i]          ,  
    params_ic_hours_female             = df_psa_input$params_ic_hours_female[i]             ,  
    params_ic_hours_age                = df_psa_input$params_ic_hours_age[i]                ,  
    params_ic_hours_T2D                = df_psa_input$params_ic_hours_T2D[i]                ,  
    beta_reg_attendance_intercept      = df_psa_input$beta_reg_attendance_intercept[i]      ,  
    beta_reg_attendance_SAF            = df_psa_input$beta_reg_attendance_SAF[i]            ,  
    beta_reg_presenteeism_intercept    = df_psa_input$beta_reg_presenteeism_intercept[i]    ,  
    beta_reg_presenteeism_SAF          = df_psa_input$beta_reg_presenteeism_SAF[i]          ,  
    LE_M_77                            = df_psa_input$LE_M_77[i]                            ,  
    LE_F_77                            = df_psa_input$LE_F_77[i]                            ,  
    c_IC_hr                            = df_psa_input$c_IC_hr[i]                            ,  
    c_IC_77_AF_2yr                     = df_psa_input$c_IC_77_AF_2yr[i]                     ,  
    du_SAF                             = df_psa_input$du_SAF[i]                             ,  
    du_CT_mo                           = df_psa_input$du_CT_mo[i]                           ,  
    du_PP_mo                           = df_psa_input$du_PP_mo[i]                           ,  
    du_VC_mo                           = df_psa_input$du_VC_mo[i]                           ,  
    mod_splines_coef                   = mod_splines_psa[i,],
    p_SFAF_CA_SAF_L1                   = p_SFAF_CA_SAF_L1_PSA[,i], 
    p_SFAF_CA_SAF_L2                   = p_SFAF_CA_SAF_L2_PSA[,i], 
    p_SFAF_CA_SAF_L3                   = p_SFAF_CA_SAF_L3_PSA[,i],
    v_age                              = patient_pop[[i]]$v_age,
    v_Sex                              = patient_pop[[i]]$v_Sex,
    df_FMC                             = df_FMC
  )       
  
  if(any(is.na(l_params))) {   
    na_params <- names(l_params)[sapply(l_params, function(x) any(is.na(x)))]   
    stop(paste("NA values detected in parameters:", paste(na_params, collapse=", "))) 
  }
  
  tic()
  # Example for your MicroSim function inside the PSA loop
  tryCatch({
    results <- foreach(g = 1:nrow(df_trtseq), 
                       .packages = c("darthtools")) %dopar% {
                         MicroSim(l_params, n_i, df_X, TRT1 = df_trtseq[g,1], TRT2 = df_trtseq[g,2], 
                                  TRT3 = df_trtseq[g,3], TRT4 = df_trtseq[g,4], TRT5 = df_trtseq[g,5], 
                                  TRT6 = df_trtseq[g,6], seed = 1)
                       }
  }, error = function(e) {
    # Log detailed information about the failed iteration
    cat("Error in iteration", i, ":\n")
    cat("Parameter values:\n")
    print(l_params)
    cat("Error message:", e$message, "\n")
    # Return NA or some baseline results instead of stopping
    return(NULL)
  })
  
  # Check if results is NULL before proceeding
  if (!is.null(results)) {
    # Process results normally
  }
  
  
  temp_v_C <- NULL
  for (j in 1:length(results)){
    temp_v_C[j]<- results[[j]]$tc_hat
  }
  
  temp_v_E <- NULL
  for (j in 1:length(results)){
    temp_v_E[j]<- results[[j]]$te_hat
  }
  
  df_c[i, ] <- c(temp_v_C)   # take the cost from the psa run and store in df_c
  df_e[i, ] <- c(temp_v_E)   # take the effect from the psa run in store in a df_e
  
  
  setTxtProgressBar(pb,i) # Display simulation progress                       
}
toc()

# Stop using multiple cores
stopCluster(cluster)

wtp <- 20000

# Save results
psa_seq_results           <- data.frame(cost   = gather(df_c, key = treatment, value = costs),  
                                        effect = gather(df_e, key = treatment, value = effect))
psa_seq_results           <- psa_seq_results[, -3]
colnames(psa_seq_results) <- c("treatment", "costs", "effects")
psa_seq_results$NHB       <- psa_seq_results$effects - (psa_seq_results$costs/wtp) 


# Save results
write.csv(psa_seq_results, here::here("output", "results_psa_ni5000_nsim5000.csv"))

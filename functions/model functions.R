# 02 Define treatment sequences and parameter input values
#### PATIENT CHARACTERISTICS ####
n_survived         <- 30182                                         # de Mol et al. 2022 (NHR 2013-2020) Table 1
n_deceased         <- 15                                            # de Mol et al. 2022 (NHR 2013-2020) Table 1
n_total            <- n_survived+n_deceased                         # de Mol et al. 2022 (NHR 2013-2020) Table 1 
age                <- ((61.4*n_survived)+(65.9*n_deceased))/n_total # de Mol et al. 2022 (NHR 2013-2020) Table 1 
sd_age             <- ((9.8*n_survived)+(7.1*n_deceased))/n_total   # de Mol et al. 2022 (NHR 2013-2020) Table 1 
n_female           <- 9816+7                                        # de Mol et al. 2022 (NHR 2013-2020) Table 1 (survived + deceased)
n_male             <- n_total-n_female                              
prop_female        <- n_female/(n_female+n_male)
n_paroxysmal       <- 19527+11     # proportion paroxysmal AF: de Mol et al. 2022 (NHR 2013-2020) Table 1 (survived + deceased)
n_persistent       <- 7300+3+542+0 # proportion persistent AF: de Mol et al. 2022 (NHR 2013-2020) Table 1 (survived + deceased)
prop_parox         <- n_paroxysmal/(n_paroxysmal+n_persistent)

# Generate patient population (will be varied in PSA)
set.seed(2)
v_age                 <- round(rtruncnorm(n_i, min_age, max_age, age, sd_age)) 
v_Sex                 <- rbinom(n_i, 1, prop_female) 

#### CLINICAL INPUTS ####

#### Probabilities of treatment success ####
load(here::here("input", "meta-analysis-efficacy.RData"))
load(here::here("input", "meta-analyse disc_AAD.RData"))

# Relative risks of AF recurrence AAD vs. CA (Source: meta-analysis Cochrane Nederland) 
RR_nai                 <- exp(predrema1$pred) # 1.76 95% CI 1.34-2.32. Patients naive to AAD and CA, i.e. first line. 
RR_AAD_exp_parox       <- exp(predrema2$pred) # 2.32 95% CI 1.78-3.01. Patients with paroxysmal AF exposed to AAD but naive to CA
RR_AAD_exp_pers        <- exp(predrema3$pred) # 1.62 95% CI 1.4-1.87.  Patients with persistent AF exposed to AAD but naive to CA
RR_CA_exp              <-  2.12 # CI 1.61-2.80, based on 1 study. Patients exposed to CA and naive or exposed to AAD.

# Probability of recurrence of AF after CA after 1 year 
prop_AF_CA_nai           <- plogis(predrema1i$pred) # 0.18 CI 0.08-0.35. After 1st CA in AAD/CA naive patients (Source: meta-analysis)
prop_AF_CA_AAD_exp_parox <- plogis(predrema2i$pred) # 0.26 CI 0.18-0.36. After 1st CA in AAD exposed paroxysmal AF patients (Source: meta-analysis)
prop_AF_CA_AAD_exp_pers  <- plogis(predrema3i$pred) # 0.38 CI 0.33-0.42. After 1st CA in AAD exposed persistent AF patients. (Source: meta-analysis)
prop_AF_CA_CA_exp_1      <- 0.2579740 # After 2nd CA. (Source: Vektis)
prop_AF_CA_CA_exp_2      <- 0.3292506 # After 3+ CA. (Source: Vektis)

#### Probability of recurrence of AF symptoms after CA over time ####
# Load transition probabilities for recurrence of AF symptoms based on data from the health insurers declaration data (Vektis)
if(scenario == "basecase"){
  load(here::here("input", "TTE_AF_Vektis.RData")) # survival model objects (version 17-04-2025)
  
  # Generate transition probabilities from survival model objects (requires flexsurv package 2.3.2)
  df_s_ablation1     <- summary(best_fit_ablation1, t = seq(0, n_t*cl, cl), ci = T, newdata = T)
  p_SFAF_CA_SAF_L1   <- trans_prob(df_s_ablation1[[1]]$est) # probability of recurrence of AF after first CA
  df_s_ablation2     <- summary(best_fit_ablation2, t = seq(0, n_t*cl, cl), , ci = T, newdata = T)
  p_SFAF_CA_SAF_L2   <- trans_prob(df_s_ablation2[[1]]$est) # probability of recurrence of AF after second CA
  df_s_ablation3     <- summary(best_fit_ablation3, t = seq(0, n_t*cl, cl), , ci = T, newdata = T)
  p_SFAF_CA_SAF_L3   <- trans_prob(df_s_ablation3[[1]]$est)
}

if(scenario == "TTE_2CV"){
  load(here::here("input", "TTE_AF_Vektis_scenario.RData")) # survival model objects (version 17-04-2025)
  
  # Generate transition probabilities from survival model objects (requires flexsurv package 2.3.2)
  df_s_ablation1     <- summary(best_fit_ablation1, t = seq(0, n_t*cl, cl), ci = T, newdata = T)
  p_SFAF_CA_SAF_L1   <- trans_prob(df_s_ablation1[[1]]$est) # probability of recurrence of AF after first CA
  df_s_ablation2     <- summary(best_fit_ablation2, t = seq(0, n_t*cl, cl), , ci = T, newdata = T)
  p_SFAF_CA_SAF_L2   <- trans_prob(df_s_ablation2[[1]]$est) # probability of recurrence of AF after second CA
  p_SFAF_CA_SAF_L3   <- p_SFAF_CA_SAF_L2 # in scenario analysis the probability of recurrence of AF after is pooled for second or more CA
}

#### Probability of death ####
# Background mortality in the general population weighted by the sex distribution of the patient population
# We used 2019 to avoid bias from COVID-19
df_mort              <- read.csv(here::here("input", "sterfte_NL_2019.csv"), sep = ";", header = TRUE)
df_mort              <- df_mort[ , -c(1,5)] # remove year and total column
df_mort <- df_mort %>% # convert table from wide to long with a variable for sex
  pivot_longer(
    cols = c(Female, Male),
    names_to = "Sex",
    values_to = "r_mort"
  ) %>%
  mutate(
    Sex = ifelse(Sex == "Female", 1, 0),
    r_mort = as.numeric(r_mort)
  )
df_mort$r_mort_cl  <- df_mort$r_mort*cl # convert annual mortality rate to x-weekly mortality rate
df_mort$p_mort_cl  <- 1-exp(-df_mort$r_mort_cl)     # convert x-weekly mortality rate to probability
df_mort[df_mort$Age == 108, "p_mort_cl"] <- 1 # set mortality rate to 1 at max age to make sure everyone dies in the model

# Excess mortality in patients with AF based on Vinter et al. (2020)
# HR of 'Model adjusted for baseline covariates' of 2001-15
# We did not use 'Model adjusted for time varying covariates' because it was no proportional hazards model
HR_EM <- 2.01 # 95% CI: 1.71 to 2.36 (Vinter et al. 2020)

#### Probability of severe adverse events ####
# Discontinuation of AAD's due to adverse events
p_disc_AAD_nai <- 0.121 # Source: meta-analysis Cochrane Netherlands
p_disc_AAD_exp <- 0.225 # Source: meta-analysis Cochrane Netherlands

# Adverse events CA (Source: event rates from NHR ablation analysis tool: Behandelgroep: Katheterablatie AF, Jaar: 2023) 
n_CT_CA         <- 26   # cardiac tamponade
n_tot_CT_CA     <- 6150
p_CT_CA         <- n_CT_CA/n_tot_CT_CA
n_PP_CA         <- 25   # phrenicus paralysis
n_tot_PP_CA     <- 6211
p_PP_CA         <- n_PP_CA/n_tot_PP_CA 
n_VC_CA         <- 40   # vascular complications
n_tot_VC_CA     <- 6149
p_VC_CA         <- n_VC_CA/n_tot_VC_CA 

#### COST INPUTS ####

#### Health state costs and costs of interventions (Source: Vektis)
df_c_states  <- load(here::here("input", "c_Vektis.RData")) 

#### Future medical costs
# All costs using PAID tool, excluding costs for "Other heart diseases including pulmonary circulation" because we assume these are already included in the model in the cost inputs above. 
df_FMC                <- read.csv(here::here("input", "PAID_AF_FMC_Unrelated_Costs_2025-02-26.csv")) 
colnames(df_FMC)      <- c("Age", "last_y_m", "last_y_f", "other_y_m", "other_y_f")
df_FMC[101:121,]      <- df_FMC[nrow(df_FMC),] # repeat at age 99 for age 100-120
df_FMC$Age[101:121]   <- seq(from = 100, to = 120, by = 1)
df_FMC <- df_FMC %>%
  pivot_longer(
    cols = c(last_y_m, last_y_f, other_y_m, other_y_f),
    names_to = c(".value", "Sex"),
    names_pattern = "(.*)_(.*)$"
  ) %>%
  mutate(
    Sex = ifelse(Sex == "f", 1, 0),
    last_y = as.numeric(last_y),
    other_y = as.numeric(other_y)
  )

#### Productivity costs inputs
# Inputs are used in PrepareCosts to calculate productivity costs
# Load the regression parameters and the variance covariance matrix of relationship between EQ-5D and productivity loss
load(here::here("input", "betareg_productivity_AVATARAF_NL.RData"))

beta_reg_attendance_intercept   <- beta_reg_attendance$coefficients$mean[1]
beta_reg_attendance_SAF         <- beta_reg_attendance$coefficients$mean[2]
beta_reg_presenteeism_intercept <- beta_reg_presenteeism$coefficients$mean[1]
beta_reg_presenteeism_SAF       <- beta_reg_presenteeism$coefficients$mean[2]

# Proportion of general population working, average hours per week and hourly wage
# Calculations performed in PrepareCosts
prop_work_M        <- 0.681 # Source: CBS 2024 Netto arbeidsparticipatie tussen 45-75 jaar.
prop_work_F        <- 0.572 # Source: CBS 2024 Netto arbeidsparticipatie tussen 45-75 jaar.
v_hours_per_week_M <- 35.9  # Source: "Werkzame beroepsbevolking; arbeidsduur" CBS 2024
v_hours_per_week_F <- 27.9  # Source: "Werkzame beroepsbevolking; arbeidsduur" CBS 2024
c_hourly_wage_2022 <- 39.88 # Source: Dutch costing manual Euro 2022. 

#### Informal care costs 
# Informal care costs are based on a regression model of de Groot et al. (2023) that estimates informal care costs according to age and proximity to death
# Regression model for probability of informal care (de Groot et al. 2023)
params_ic_log_intercept   <- -1.451 # SE 0.220  
params_ic_log_female      <-  0.368 # SE 0.120  
params_ic_log_age         <-  0.054 # SE 0.008  
params_ic_log_age2        <-  0.000 # SE 0.000 
params_ic_log_T2D         <- -0.061 # SE 0.015  
params_ic_hours_intercept <-  0.497 # SE 0.198  
params_ic_hours_female    <-  0.112 # SE 0.106  
params_ic_hours_age       <-  0.019 # SE 0.005  
params_ic_hours_T2D       <- -0.034 # SE 0.011 

# Logistic regression for use of informal care
# Coefficients on log-scale: intercept, gender (female = 1), age (age-centered at 70), age2, TTD (in years)
params_ic_log   <- c(params_ic_log_intercept, params_ic_log_female, params_ic_log_age, params_ic_log_age2, params_ic_log_T2D) 
vcov_ic_log     <- read.csv(here::here("input", "vcov_ic_log.csv"), sep = ",", header = FALSE)

# Lineair regression for number of hours of informal care
# Coefficients: intercept, gender (female = 1), age  (age-centered at 70), TTD (in years)
params_ic_hours <- c(params_ic_hours_intercept, params_ic_hours_female, params_ic_hours_age, params_ic_hours_T2D)
vcov_ic_hours   <- read.csv(here::here("input", "vcov_ic_hours.csv"), sep = ",", header = FALSE)

# Inputs to calculate the average time to death and corresponding informal care costs of patients with the age of Van Den Dries et al. (i.e. 77 years)
# Calculations performed in PrepareCosts
LE_M_77              <- 10.3 # Life expectancy males of 77 years, Source: CBS
LE_F_77              <- 11.9 # Life expectancy females of 77 years, Source: CBS
c_IC_hr              <- 18.8 # Hourly costs of informal care (Source: Kostenhandleiding, 2022 Euros)
c_IC_77_AF_2yr       <- 3296.45 #SD = 321.76, 24-month costs of AF patients in Van den Dries et al. (2023) 
c_IC_77_AF_2yr_gamma <- gamma_params(c_IC_77_AF_2yr, (321.76/sqrt(425)))

#### UTILITIES INPUTS ####
#### Disutilities AF and adverse events ####
# Load utilities based on the EQ-5D-5L data in the AVATAR-AF trial using the Dutch tariff
load(here::here("input", "glm_AVATARAF_NL.RData"))

du_SAF       <- abs(unname((glm_AVATARAF_NL$coefficients[2]))) # Re-analysis of AVATAR-AF with EQ-5D-5L NL tariff based on Moss et al. 

# Assumption similar to Akerborg et al. (2012) and Reynolds et al. (2014)
du_CT_mo  <- 0.1 # Assumption disutility of 0.1 for 1 month, corrected to disutility for the cycle length in Effs
du_PP_mo  <- 0.1 # Assumption disutility of 0.1 for 1 month, corrected to disutility for the cycle length in Effs
du_VC_mo  <- 0.1 # Assumption disutility of 0.1 for 1 month, corrected to disutility for the cycle length in Effs

#### General population utilities ####
# Regression model general population by age and sex (based on Versteegh et al. 2016)
load(here::here("input", "splines_HRQoL.RData")) 
gen_pop_utility  <- "fit"
mod_splines_coef <- mod_splines$coefficients

#### STORE INPUT PARAMETERS ####
# Create a vector of variable names
v_names_params <- c("age", "sd_age", "prop_female", "v_age", "v_Sex", "prop_parox",
                    "RR_nai", "RR_AAD_exp_parox", "RR_AAD_exp_pers", "RR_CA_exp", 
                    "prop_AF_CA_nai", "prop_AF_CA_AAD_exp_parox","prop_AF_CA_AAD_exp_pers", 
                    "prop_AF_CA_CA_exp_1", "prop_AF_CA_CA_exp_2",   
                    "p_SFAF_CA_SAF_L1", "p_SFAF_CA_SAF_L2", "p_SFAF_CA_SAF_L3", 
                    "HR_EM", "p_disc_AAD_nai", "p_disc_AAD_exp",
                    "p_CT_CA", "p_PP_CA", "p_VC_CA",
                    "c_AAD", "c_CA_DBC", "c_d_before_CA", "c_d_after_CA1", "c_d_after_CA2", "c_d_SAF", "df_FMC", 
                    "prop_work_M", "prop_work_F", "v_hours_per_week_M", "v_hours_per_week_F",
                    "c_hourly_wage_2022", "params_ic_log_intercept", "params_ic_log_female", 
                    "params_ic_log_age", "params_ic_log_age2", 'params_ic_log_T2D',  "params_ic_hours_intercept",
                    "params_ic_hours_female", "params_ic_hours_age", "params_ic_hours_T2D",   
                    "beta_reg_attendance_intercept", "beta_reg_attendance_SAF", 
                    "beta_reg_presenteeism_intercept", "beta_reg_presenteeism_SAF",      
                    "LE_M_77", "LE_F_77", "c_IC_hr", "c_IC_77_AF_2yr",
                    "du_SAF", "du_CT_mo", "du_PP_mo", "du_VC_mo", 
                    "gen_pop_utility", "mod_splines_coef"
)

# Store the parameters into a list
l_params <- mget(v_names_params)
#View(l_params)

# 06 Microsimulation model

## 06.1 df_X: dataframe with individual level characteristics
Create_df_X <- function(l_params){
  with((l_params), {
    
    set.seed(1)
    
    # v_age and v_Sex are created in the previous section to enable variation in PSA
    v_Trt      <- rep("AAD", n_i) # Placeholder that will be overwritten in first section of MicroSim
    v_Line     <- rep(1, n_i)
    v_disc_AAD <- rep(0, n_i)
    v_CT_CA    <- rep(0, n_i)  # event counter for cardiac tamponade
    v_PP_CA    <- rep(0, n_i)  # event counter for phrenicus paralysis
    v_VC_CA    <- rep(0, n_i)  # event counter for vascular complications
    v_TH       <- rep(0, n_i)  # Treatment history: 1 if current treatment is first month of CA (for side effects)
    v_AAD_exp  <- rep(0, n_i)  # Exposed to AAD
    v_CA_exp_1 <- rep(0, n_i)  # Exposed to CA
    v_CA_exp_2 <- rep(0, n_i)  # Exposed to at least 2 CAs
    v_death    <- rep(0, n_i)  # Death
    v_time     <- rep(1, n_i)  # Time (for extracting the first transition probability from the K-M curves)
    
    df_X       <- data.frame(ID = 1:n_i, Age = v_age, Age_cl = v_age, Age_start = v_age, Sex = v_Sex, time = v_time,
                             curTrt = v_Trt, Line = v_Line, 
                             disc_AAD = v_disc_AAD, 
                             CT_CA = v_CT_CA, PP_CA = v_PP_CA, VC_CA = v_VC_CA, 
                             TH = v_TH, AAD_exp = v_AAD_exp, CA_exp_1 = v_CA_exp_1, CA_exp_2 = v_CA_exp_2, 
                             death = v_death)
    return(df_X)
  })
}

## 06.2 Probabilities
### 06.2.1 Prepare probabilities 

PrepareProbs <- function(l_params){
  
  with((l_params),{
    
    # Calculate weighted average for prop_AF_CA_AAD_exp and RR_AAD_exp
    prop_AF_CA_AAD_exp     <- expit((logit(prop_AF_CA_AAD_exp_parox)*prop_parox)+(logit(prop_AF_CA_AAD_exp_pers)*(1-prop_parox)))  
    
    # Probability of recurrence of AF after AAD after 1 year 
    RR_AAD_exp             <- exp((log(RR_AAD_exp_parox)*prop_parox)+(log(RR_AAD_exp_pers)*(1-prop_parox))) # weighted RR 
    
    # Based on proportion of recurrence of AF after CA multiplied with relative risk of AAD vs. CA
    prop_AF_AAD_nai        <- prop_AF_CA_nai*RR_nai    
    prop_AF_AAD_AAD_exp    <- prop_AF_CA_AAD_exp*RR_AAD_exp
    prop_AF_AAD_CA_exp_1   <- prop_AF_CA_CA_exp_1*RR_CA_exp
    prop_AF_AAD_CA_exp_2   <- prop_AF_CA_CA_exp_2*RR_CA_exp
    
    # To prevent probabilities > 1 (not needed for deterministic analyses, only occurs in some iterations of the PSA)
    prop_AF_CA_nai        <- ifelse(prop_AF_CA_nai       > 1, 1, prop_AF_CA_nai)       
    prop_AF_CA_AAD_exp    <- ifelse(prop_AF_CA_AAD_exp   > 1, 1, prop_AF_CA_AAD_exp)   
    prop_AF_CA_CA_exp_1   <- ifelse(prop_AF_CA_CA_exp_1  > 1, 1, prop_AF_CA_CA_exp_1)  
    prop_AF_CA_CA_exp_2   <- ifelse(prop_AF_CA_CA_exp_2  > 1, 1, prop_AF_CA_CA_exp_2)  
    prop_AF_AAD_nai       <- ifelse(prop_AF_AAD_nai      > 1, 1, prop_AF_AAD_nai     ) 
    prop_AF_AAD_AAD_exp   <- ifelse(prop_AF_AAD_AAD_exp  > 1, 1, prop_AF_AAD_AAD_exp ) 
    prop_AF_AAD_CA_exp_1  <- ifelse(prop_AF_AAD_CA_exp_1 > 1, 1, prop_AF_AAD_CA_exp_1) 
    prop_AF_AAD_CA_exp_2  <- ifelse(prop_AF_AAD_CA_exp_2 > 1, 1, prop_AF_AAD_CA_exp_2) 
    
    # Convert annual probabilities to cycle length
    p_SFAF_CA_nai          <- 1-convert_probability_to_cl(prop_AF_CA_nai,       1, cl) 
    p_SFAF_CA_AAD_exp      <- 1-convert_probability_to_cl(prop_AF_CA_AAD_exp,   1, cl) 
    p_SFAF_CA_CA_exp_1     <- 1-convert_probability_to_cl(prop_AF_CA_CA_exp_1,  1, cl) 
    p_SFAF_CA_CA_exp_2     <- 1-convert_probability_to_cl(prop_AF_CA_CA_exp_2,  1, cl) 
    p_SFAF_AAD_nai         <- 1-convert_probability_to_cl(prop_AF_AAD_nai,      1, cl) 
    p_SFAF_AAD_AAD_exp     <- 1-convert_probability_to_cl(prop_AF_AAD_AAD_exp,  1, cl) 
    p_SFAF_AAD_CA_exp_1    <- 1-convert_probability_to_cl(prop_AF_AAD_CA_exp_1, 1, cl) 
    p_SFAF_AAD_CA_exp_2    <- 1-convert_probability_to_cl(prop_AF_AAD_CA_exp_2, 1, cl) 
    
    l_inputs_probs_names <- c("RR_AAD_exp", "p_SFAF_CA_nai", "p_SFAF_CA_AAD_exp", "p_SFAF_CA_CA_exp_1", "p_SFAF_CA_CA_exp_2",
                              "p_SFAF_AAD_nai", "p_SFAF_AAD_AAD_exp", "p_SFAF_AAD_CA_exp_1", "p_SFAF_AAD_CA_exp_2")
    l_inputs_probs <- mget(l_inputs_probs_names)
    
    return(l_inputs_probs)
  })
}

### 06.2.2 Probability function
Probs <- function(l_params_all, M_t, df_X, t) { 
  # Arguments:
  # M_t:   health state occupied by individual i at cycle t (character variable)
  # df_X:  data frame with individual characteristics data 
  # t:     current cycle 
  # Returns: 
  # transition probabilities for that cycle
  
  with((l_params_all),{
    
    # Create matrix of state transition probabilities  
    m_p_t           <- matrix(0, nrow = n_states, ncol = n_i) 
    rownames(m_p_t) <-  v_n  # give the state names to the rows
    
    # Look up probability of dying based on current age
    df_p_D          <- inner_join(df_X, df_mort, by = c("Age", "Sex")) 
    
    # Multiply with HR for excess mortality and cycle length and convert back to probability
    p_D             <- 1 - (1-df_p_D$p_mort_cl)^HR_EM 
    
    # Success rates
    p_SAF_SFAF                                                                <- NULL
    p_SAF_SFAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 0] <- p_SFAF_AAD_nai
    p_SAF_SFAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 0] <- p_SFAF_AAD_AAD_exp
    p_SAF_SFAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 1] <- p_SFAF_AAD_CA_exp_1
    p_SAF_SFAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 1] <- p_SFAF_AAD_CA_exp_1
    p_SAF_SFAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_2 == 1] <- p_SFAF_AAD_CA_exp_2
    p_SAF_SFAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_2 == 1] <- p_SFAF_AAD_CA_exp_2
    
    p_SAF_SFAF[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 0]  <- p_SFAF_CA_nai
    p_SAF_SFAF[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 0]  <- p_SFAF_CA_AAD_exp
    p_SAF_SFAF[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 1]  <- p_SFAF_CA_CA_exp_1
    p_SAF_SFAF[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 1]  <- p_SFAF_CA_CA_exp_1
    p_SAF_SFAF[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_2 == 1]  <- p_SFAF_CA_CA_exp_2
    p_SAF_SFAF[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_2 == 1]  <- p_SFAF_CA_CA_exp_2
    p_SAF_SFAF[df_X$curTrt == "no-treatment"]                                 <- 0 # keep symptoms when no rhythm control treatment
    
    # Recurrence rates at time t
    p_SFAF_SAF                                                                 <- NULL
    
    p_SFAF_SAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 0]  <- ifelse(p_SFAF_CA_SAF_L1[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 0]]*RR_nai     >1, 1, 
                                                                                         p_SFAF_CA_SAF_L1[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 0]]*RR_nai    )
    p_SFAF_SAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 0]  <- ifelse(p_SFAF_CA_SAF_L1[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 0]]*RR_AAD_exp >1, 1, 
                                                                                         p_SFAF_CA_SAF_L1[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 0]]*RR_AAD_exp)
    p_SFAF_SAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 1]  <- ifelse(p_SFAF_CA_SAF_L2[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 1]]*RR_CA_exp  >1, 1, 
                                                                                         p_SFAF_CA_SAF_L2[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 1]]*RR_CA_exp )
    p_SFAF_SAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 1]  <- ifelse(p_SFAF_CA_SAF_L2[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 1]]*RR_CA_exp  >1, 1, 
                                                                                         p_SFAF_CA_SAF_L2[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 1]]*RR_CA_exp )
    p_SFAF_SAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_2 == 1]  <- ifelse(p_SFAF_CA_SAF_L3[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_2 == 1]]*RR_CA_exp  >1, 1, 
                                                                                         p_SFAF_CA_SAF_L3[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 0 & df_X$CA_exp_2 == 1]]*RR_CA_exp )
    p_SFAF_SAF[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_2 == 1]  <- ifelse(p_SFAF_CA_SAF_L3[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_2 == 1]]*RR_CA_exp  >1, 1, 
                                                                                         p_SFAF_CA_SAF_L3[df_X$time[df_X$curTrt == "AAD" & df_X$AAD_exp == 1 & df_X$CA_exp_2 == 1]]*RR_CA_exp )
    
    p_SFAF_SAF[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 0]  <- p_SFAF_CA_SAF_L1[df_X$time[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 0]]
    p_SFAF_SAF[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 0]  <- p_SFAF_CA_SAF_L1[df_X$time[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 0]]
    p_SFAF_SAF[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 1]  <- p_SFAF_CA_SAF_L2[df_X$time[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_1 == 1]]
    p_SFAF_SAF[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 1]  <- p_SFAF_CA_SAF_L2[df_X$time[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_1 == 1]]
    p_SFAF_SAF[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_2 == 1]  <- p_SFAF_CA_SAF_L3[df_X$time[df_X$curTrt == "CA" & df_X$AAD_exp == 0 & df_X$CA_exp_2 == 1]]
    p_SFAF_SAF[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_2 == 1]  <- p_SFAF_CA_SAF_L3[df_X$time[df_X$curTrt == "CA" & df_X$AAD_exp == 1 & df_X$CA_exp_2 == 1]]
    
    p_SFAF_SAF[df_X$curTrt == "no-treatment"]                             <- 1-p_D[df_X$curTrt == "no-treatment"] # get symptoms when no rhythm control treatment
    
    # Fill the transition probability matrix with the appropriate probabilities
    m_p_t[, M_t == "SFAF"]     <- rbind((1-p_SFAF_SAF[M_t == "SFAF"]-p_D[M_t == "SFAF"]),
                                        p_SFAF_SAF[M_t == "SFAF"], 
                                        p_D[M_t == "SFAF"])
    m_p_t[, M_t == "SAF"]      <- rbind(p_SAF_SFAF[M_t == "SAF"], 
                                        1-p_SAF_SFAF[M_t == "SAF"]-p_D[M_t == "SAF"], 
                                        p_D[M_t == "SAF"])
    m_p_t[, M_t == "D"]        <- rbind(0, 0, 1)
    
    if(any(is.na(m_p_t))) {
      warning("NA values found in transition probability matrix")
      print(which(is.na(m_p_t), arr.ind = TRUE))
    }
    
    return(t(m_p_t))
  }) # End of with l_params
}  


## 06.3 Costs
### 06.3.1 Prepare costs
PrepareCosts <- function(l_params){
  
  with((l_params),{
    
    #### Health state costs #### (Source: Vektis)
    c_before_CA   <- adjust_inflation(c_d_before_CA, "2021")
    c_after_CA    <- adjust_inflation(c_d_after_CA1, "2021") + adjust_inflation(c_d_after_CA2, "2021")
    c_d_SAF       <- adjust_inflation(c_d_SAF, "2021")
    
    #### Intervention costs ####
    # Based on average AAD costs in patients with symptoms of AF (12-6 months before CA)
    c_AAD   <- adjust_inflation(c_AAD, "2021")
    
    # Intervention costs + additional  health care consumption 6 months before and 12 months after CA
    c_CA    <- adjust_inflation(c_CA_DBC, "2021")+c_before_CA+c_after_CA 
    c_CA    <- unname(c_CA)
    
    #### Future medical costs ####
    df_FMC$other_y <- df_FMC$other_y*cl # Adapt to the cycle length
    df_FMC$last_y  <- adjust_inflation(df_FMC$last_y, "2017")
    df_FMC$other_y <- adjust_inflation(df_FMC$other_y, "2017")
    df_FMC         <- as.data.frame(df_FMC)
    
    #### Productivity costs ####
    # Productivity costs consist of three components in this model, recovery time after CA (1 week), work attendance and presenteeism
    
    # Work attendance and presenteeism are estimated with regression models from the EQ-5D-5L to 3L converted health states, age (work attendance) and sex (presenteeism) (Source: Krol et al. 2014). We used the AVATAR-AF data to estimate work attendance and presenteeism for each AF patient in the trial. Then we used a beta-regression to estimate how symptom status influenced work attendance and presenteeism.
    
    # Probability of work attendance based on symptoms. Source: AVATAR-AF and Krol et al. (2014)
    p_SFAF_workattendance <- exp(beta_reg_attendance_intercept)/(1+exp(beta_reg_attendance_intercept))
    p_SAF_workattendance  <- exp((beta_reg_attendance_intercept+beta_reg_attendance_SAF))/
      (1+exp(beta_reg_attendance_intercept+beta_reg_attendance_SAF))
    
    # Proportion at work based on general population corrected for age and sex distribution in AF population
    prop_work             <- (prop_work_M*(1-prop_female)) + prop_work_F*prop_female
    
    # Correct proportion at work for work attendance of AF patients with/without symptoms
    prop_work_SAF         <- prop_work*p_SAF_workattendance # Proportion of patients at work with AF symptoms
    prop_work_SFAF        <- prop_work*p_SFAF_workattendance # Proportion of patients at work without AF symptoms
    
    # Probability of presenteeism (i.e. productivity at work)
    p_SFAF_presenteeism   <- exp(beta_reg_presenteeism_intercept)/(1+exp(beta_reg_presenteeism_intercept))
    p_SAF_presenteeism    <- exp((beta_reg_presenteeism_intercept+beta_reg_presenteeism_SAF))/
      (1+exp(beta_reg_presenteeism_intercept+beta_reg_presenteeism_SAF))
    
    # Productivity costs based on average work duration and wage
    v_hours_per_week      <- v_hours_per_week_F*prop_female+(v_hours_per_week_M*(1-prop_female)) 
    v_hours_per_cycle     <- cl*(365.25/7)*v_hours_per_week # Work hours in a cycle
    c_hourly_wage         <- adjust_inflation(c_hourly_wage_2022, "2022") 
    c_prod_cycle          <- v_hours_per_cycle*c_hourly_wage
    
    # Total productivity: work attendance corrected for presenteeism multiplied with costs
    c_prod_SAF_tot        <- prop_work_SAF*p_SAF_presenteeism*c_prod_cycle # Total productivity costs patients with AF symptoms
    c_prod_SFAF_tot       <- prop_work_SFAF*c_prod_cycle*p_SFAF_presenteeism # Total productivity costs patients without AF symptoms
    c_prod_SAF            <- unname(c_prod_SFAF_tot-c_prod_SAF_tot) # Total productivity costs attributable to having AF symptoms
    
    # Productivity costs of recovery after CA: 1 week absent from work
    c_prod_CA_recovery    <- v_hours_per_week*c_hourly_wage*prop_work_SAF
    
    #### Informal care costs ####
    # Informal care costs are based on a regression model of de Groot et al. (2023) that estimates informal care costs according to age and proximity to death
    # Weighted average life expectancy AF gender distribution, life expectancy is equal to time to death
    T2D              <- (LE_M_77*(1-prop_female))+(LE_F_77*prop_female) 
    
    # Calculate informal care costs for a 77-year old from the general population 
    # Age is centered at the mean in the analysis and 70 is the mean age in de Groot et al. (2023) 
    p_IC_log         <- params_ic_log_intercept + params_ic_log_female*prop_female + params_ic_log_age*(77-70) + 
      params_ic_log_age2*((77-70)*(77-70)) + params_ic_log_T2D*T2D # Propotion use of informal care
    p_IC             <- exp(p_IC_log)/(1+exp(p_IC_log)) # Convert to probability
    v_hr_IC_log      <- params_ic_hours_intercept + params_ic_hours_female*prop_female + params_ic_hours_age*(77-70) + params_ic_hours_T2D*T2D 
    v_hr_day_IC      <- exp(v_hr_IC_log) # Convert from log scale to hours per day
    c_IC_77_cl       <- p_IC*(v_hr_day_IC*365.25*cl)*c_IC_hr # Total costs of informal care in 77-year old in a cycle
    c_IC_77_AF_cl    <- (c_IC_77_AF_2yr/2)*cl # correct for cycle length
    c_IC_SAF         <- c_IC_77_AF_cl-c_IC_77_cl # Informal care costs attributable to symptomatic AF
    
    l_inputs_costs_names <- c("c_d_SAF", "c_AAD", "c_CA", "df_FMC", "c_prod_SAF", "c_prod_CA_recovery", "c_IC_SAF")
    
    l_inputs_costs <- mget(l_inputs_costs_names)
    
    return(l_inputs_costs)
  })
}



### 06.3.2 Cost function
Costs <- function (l_params_all, M_t, df_X) {
  # M_t: health state occupied by individual i at cycle t (character variable)
  
  with((l_params_all),{
    
    # Objects for future medical costs
    c_FMC   <- inner_join(df_X, df_FMC, by = c("Age", "Sex")) # Look up future medical costs based on current age and sex
    
    # Assign costs based on health states (informal care costs will be added at the end of the MicroSim)
    c_t                                                 <- NULL 
    
    c_t[M_t == "SFAF" & df_X$curTrt == "AAD"]           <- c_AAD + c_FMC$other_y[M_t == "SFAF" & df_X$curTrt == "AAD"] 
    
    c_t[M_t == "SFAF" & df_X$curTrt == "CA"]            <- c_CA * df_X$TH[M_t == "SFAF" & df_X$curTrt == "CA"] +  
      c_FMC$other_y[M_t == "SFAF" & df_X$curTrt == "CA"] +
      (c_prod_CA_recovery * df_X$TH[M_t == "SFAF" & df_X$curTrt == "CA"]) * 
      (df_X$Age[M_t == "SFAF"  & df_X$curTrt == "CA"]<68)
    
    c_t[M_t == "SFAF" & df_X$curTrt == "no-treatment"]  <- c_FMC$other_y[M_t == "SFAF" & df_X$curTrt == "no-treatment"] 
    
    c_t[M_t == "SAF"  & df_X$curTrt == "AAD"]           <- c_AAD + c_d_SAF + c_FMC$other_y[M_t == "SAF" & df_X$curTrt == "AAD"] +
      c_prod_SAF*(df_X$Age[M_t == "SAF"  & df_X$curTrt == "AAD"]<68) 
    
    c_t[M_t == "SAF"  & df_X$curTrt == "CA"]            <- c_CA + c_d_SAF + c_FMC$other_y[M_t == "SAF" & df_X$curTrt == "CA"] +
      c_prod_SAF*(df_X$Age[M_t == "SAF"  & df_X$curTrt == "CA"]<68) + 
      (c_prod_CA_recovery * df_X$TH[M_t == "SAF" & df_X$curTrt == "CA"]) * 
      (df_X$Age[M_t == "SAF"  & df_X$curTrt == "CA"]<68) 
    
    c_t[M_t == "SAF"  & df_X$curTrt == "no-treatment"]  <- c_d_SAF + c_FMC$other_y[M_t == "SAF" & df_X$curTrt == "no-treatment"] +
      c_prod_SAF*(df_X$Age[M_t == "SAF"  & df_X$curTrt == "no-treatment"]<68) 
    
    c_t[M_t == "D"]                                     <- c_FMC$last_y[M_t == "D"]*df_X$death[M_t == "D"]
    
    return(c_t)        		                      # return the costs
  }) # end of with l_params_all
}

## 06.4 Health outcome function

Effs <- function (l_params_all, M_t, df_X, t, cl = cl) {
  # M_t: health state occupied by individual i at cycle t (character variable)
  # df_X: data frame with individual characteristics data 
  
  with((l_params_all),{
    
    # Determine distutilites of adverse events after CA
    du_CT    <- du_CT_mo*cl # Assumption disutility of 0.1 for 1 month, corrected to disutility for the cycle length
    du_PP    <- du_PP_mo*cl # Assumption disutility of 0.1 for 1 month, corrected to disutility for the cycle length
    du_VC    <- du_VC_mo*cl # Assumption disutility of 0.1 for 1 month, corrected to disutility for the cycle length
    
    # Create vector with total disutility per patient based on df_X
    du_AE_CA <- NULL
    du_AE_CA <- (df_X$CT_CA*du_CT)+(df_X$PP_CA*du_PP)+(df_X$VC_CA*du_VC)   
    
    # Determine baseline utility based on general population estimates
    mod_splines$coefficients <- mod_splines_coef # included for use in PSA
    
    if(gen_pop_utility == "fit"){ # For use in base case
      u_t <- unname(predict(mod_splines, newdata = df_X)) 
    }
    
    if(gen_pop_utility == "lwr"){ # For use in OWSA
      u_t <- unname(predict(mod_splines, newdata = df_X, interval = "confidence")[, "lwr"]) 
    }
    
    if(gen_pop_utility == "upr"){ # For use in OWSA
      u_t <- unname(predict(mod_splines, newdata = df_X, interval = "confidence")[, "upr"]) 
    }
    
    # Assign utilities to health states
    u_t[M_t == "SFAF"] <- u_t[M_t == "SFAF"] - du_AE_CA[M_t == "SFAF"] 
    u_t[M_t == "SAF"]  <- u_t[M_t == "SAF"] - du_AE_CA[M_t == "SAF"] - du_SAF 
    u_t[M_t == "D"]    <- 0
    
    QALYs <-  u_t * cl  # calculate the QALYs during cycle t
    
    return(QALYs)       # return the QALYs
    
  }) # end of with l_params_all
}

## 06.5 Microsimulation

MicroSim <- function(l_params, n_i, df_X, TRT1 = TRT1, TRT2 = TRT2, TRT3 = TRT3, 
                     TRT4 = TRT4, TRT5 = TRT5, TRT6 = TRT6, seed = 1) {
  
  # Arguments:  
  # n_i:     number of individuals
  # df_X     data frame with individual characteristics data 
  # TRT1-3:  The treatments in the sequence
  # seed:    default is 1
  
  #### Set up starting values ####
  df_X        <- Create_df_X(l_params)             # Create dataframe with patient characteristics
  df_X$curTrt <- TRT1                              # Assign the first treatment to all patients
  v_M_init    <- rep("SAF", n_i)                   # All patients start with AF symptoms 
  df_X$TH     <- ifelse(df_X$curTrt == "CA", 1, 0) # If a patient starts on CA, treatment history (TH) is set to 1 
  
  # Calculate input parameters and combine them in l_params_all
  l_inputs_Probs <- PrepareProbs(l_params)
  l_inputs_Costs <- PrepareCosts(l_params)
  l_params_all <- c(l_params, l_inputs_Probs, l_inputs_Costs)
  l_params_all <- l_params_all[!duplicated(names(l_params_all), fromLast = T)]
  
  with((l_params_all), {
    set.seed(seed) # set the seed
    n_states <- length(v_n) # the number of health states
    
    #### Create matrices ####
    # create matrices with number of rows equal to the n_i, the number of columns equal to n_t  
    # (the initial state and all the n_t cycles)
    # m_M is used to store the health state information over time for every individual
    # m_C is used to store the costs information over time for every individual
    # m_E is used to store the effects information over time for every individual
    # m_L is used to store the line information over time for every individual
    # m_TH is used to store the treatment history of CA over time for every individual
    # m_AE is used to store the adverse events after CA over time for every individual
    # m_D is used to store the death status over time for every individual; used to calculate future medical costs
    # m_C_IC_cl is used to store the informal care costs over time for every individual
    # m_curTrt is used to store the current treatment over time for every individual
    # m_Age is used to store the current age over time for every individual
    
    m_M <- m_C <- m_E <- m_L <- m_TH <- m_AE <- m_D <- m_C_IC_cl <- m_curTrt <- m_Age <- matrix(nrow = n_i, ncol = n_t + 1, 
                                                                                                dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                                                                                paste("cycle", 0:n_t, sep = " ")))  
    
    m_M [, 1]      <- v_M_init    # initial health state at cycle 0 for individual i
    m_TH[, 1]      <- df_X$TH     # initial treatment history at cycle 0 for individual i
    m_L [, 1]      <- rep(1, n_i) # initial line at cycle 0 for individual i
    m_D [, 1]      <- df_X$death  # initial death status at cycle 0 for individual i
    m_C_IC_cl[, 1] <- rep(0, n_i) # initial costs of informal care
    m_curTrt[, 1]  <- df_X$curTrt # initial treatment
    m_Age[, 1]     <- df_X$Age    # initial age
    
    #### Cycle 0 ####
    # Costs and QALYs in cycle 1
    m_C[, 1]  <- Costs(l_params_all, m_M[, 1], df_X)     
    m_E[, 1]  <- Effs (l_params_all, m_M[, 1], df_X, t = 1, cl = cl) 
    
    #### Start loop cycle 1 to n_t ####
    # Open a loop for time running cycles 1 to n_t 
    for (t in 1:n_t) {
      
      # To remove variability due to random draw procedure (seed) but keep variation between cycles (+ t)
      set.seed(seed + t)
      
      #### Switch health states ####
      # Calculate the transition probabilities for the cycle based on  health state t
      m_P <- Probs(l_params_all, m_M[, t], df_X, t)             
      
      # Sample the current health state based on the transition probabilities and store that state in matrix m_M 
      m_M[, t + 1]  <- samplev(m_P, 1)   
      
      #### Adverse events ####
      # Discontinue AAD due to adverse events
      df_X$disc_AAD[df_X$curTrt == "AAD" & df_X$AAD_exp == 0] <- rbinom(sum(df_X$curTrt == "AAD" & df_X$AAD_exp == 0), 1, p_disc_AAD_nai)
      df_X$disc_AAD[df_X$curTrt == "AAD" & df_X$AAD_exp == 1] <- rbinom(sum(df_X$curTrt == "AAD" & df_X$AAD_exp == 1), 1, p_disc_AAD_exp)
      
      # Adverse events during the first cycle after CA 
      df_X$PP_CA[df_X$curTrt == "CA" & df_X$TH == 1] <- rbinom(sum(df_X$curTrt == "CA" & df_X$TH == 1), 1, p_PP_CA)
      df_X$CT_CA[df_X$curTrt == "CA" & df_X$TH == 1] <- rbinom(sum(df_X$curTrt == "CA" & df_X$TH == 1), 1, p_CT_CA)
      df_X$VC_CA[df_X$curTrt == "CA" & df_X$TH == 1] <- rbinom(sum(df_X$curTrt == "CA" & df_X$TH == 1), 1, p_VC_CA)
      
      m_AE[, t+ 1] <- df_X$PP_CA + df_X$CT_CA +df_X$VC_CA 
      
      #### Switching lines ####
      # Switch lines when symptoms of AF or discontinuation of AAD
      # If maximum number of lines has been reached, switch to no rhythm control treatment
      df_X$Line[m_M[ , t + 1] == "SFAF" & df_X$disc_AAD == 0]                   <- df_X$Line[m_M[ , t + 1] == "SFAF" & df_X$disc_AAD == 0] 
      df_X$Line[m_M[ , t + 1] == "SFAF" & df_X$disc_AAD == 1 & df_X$Line >= 6]  <- 555 # no rhythm control treatment 
      df_X$Line[m_M[ , t + 1] == "SFAF" & df_X$disc_AAD == 1 & df_X$Line < 6]   <- df_X$Line[m_M[ , t + 1] == "SFAF"& df_X$disc_AAD == 1 & df_X$Line < 6] + 1  
      df_X$Line[m_M[ , t + 1] == "SAF" & df_X$Line >= 6]                        <- 555 # no rhythm control treatment
      df_X$Line[m_M[ , t + 1] == "SAF" & df_X$Line < 6]                         <- df_X$Line[m_M[ , t + 1] == "SAF" & df_X$Line < 6] + 1
      df_X$Line[m_M[ , t + 1] == "D"]                                           <- 999 # dead
      
      # Update matrix with current line number
      m_L[ , t + 1] <- df_X$Line
      
      # Update current treatment for the next cycle
      df_X$curTrt[df_X$Line == 1]   <- TRT1
      df_X$curTrt[df_X$Line == 2]   <- TRT2
      df_X$curTrt[df_X$Line == 3]   <- TRT3
      df_X$curTrt[df_X$Line == 4]   <- TRT4
      df_X$curTrt[df_X$Line == 5]   <- TRT5
      df_X$curTrt[df_X$Line == 6]   <- TRT6
      df_X$curTrt[df_X$Line == 555] <- "no-treatment"
      df_X$curTrt[df_X$Line == 999] <- "D"
      
      m_curTrt[, t + 1]  <- df_X$curTrt
      
      # Update the evaluation of first month of ablation
      df_X$TH       <- ifelse(df_X$curTrt == "CA" & m_L[,t] != m_L[,t+1], 1, 0)
      m_TH[, t + 1] <- df_X$TH
      
      # Update variables for AAD or CA exposure
      # Exposed to AAD
      df_X$AAD_exp[TRT1 == "AAD" & df_X$Line > 1] <- 1
      df_X$AAD_exp[TRT2 == "AAD" & df_X$Line > 2] <- 1
      df_X$AAD_exp[TRT3 == "AAD" & df_X$Line > 3] <- 1
      df_X$AAD_exp[TRT4 == "AAD" & df_X$Line > 4] <- 1
      df_X$AAD_exp[TRT5 == "AAD" & df_X$Line > 5] <- 1
      
      # Exposed to CA at least once
      df_X$CA_exp_1[TRT1 == "CA" & df_X$Line > 1]   <- 1
      df_X$CA_exp_1[TRT2 == "CA" & df_X$Line > 2]   <- 1
      df_X$CA_exp_1[TRT3 == "CA" & df_X$Line > 3]   <- 1
      df_X$CA_exp_1[TRT4 == "CA" & df_X$Line > 4]   <- 1
      df_X$CA_exp_1[TRT5 == "CA" & df_X$Line > 5]   <- 1
      
      # Exposed to CA at least twice
      df_X$CA_exp_2[TRT1 == "CA" & TRT2 == "CA" & df_X$Line > 2] <- 1
      df_X$CA_exp_2[TRT1 == "CA" & TRT3 == "CA" & df_X$Line > 3] <- 1
      df_X$CA_exp_2[TRT1 == "CA" & TRT4 == "CA" & df_X$Line > 4] <- 1
      df_X$CA_exp_2[TRT1 == "CA" & TRT5 == "CA" & df_X$Line > 5] <- 1
      df_X$CA_exp_2[TRT2 == "CA" & TRT3 == "CA" & df_X$Line > 3] <- 1
      df_X$CA_exp_2[TRT2 == "CA" & TRT4 == "CA" & df_X$Line > 4] <- 1
      df_X$CA_exp_2[TRT2 == "CA" & TRT5 == "CA" & df_X$Line > 5] <- 1
      df_X$CA_exp_2[TRT3 == "CA" & TRT4 == "CA" & df_X$Line > 4] <- 1
      df_X$CA_exp_2[TRT3 == "CA" & TRT5 == "CA" & df_X$Line > 5] <- 1
      df_X$CA_exp_2[TRT4 == "CA" & TRT5 == "CA" & df_X$Line > 5] <- 1
      
      # Set adverse events status back to zero for everyone 
      df_X$disc_AAD <- df_X$CT_CA <- df_X$PP_CA <- df_X$VC_CA <- 0 
      
      # Update the age of individuals that are alive
      df_X$Age_cl[m_M[, t + 1] != "D"]  <- df_X$Age_cl[m_M[, t + 1] != "D"] + cl
      df_X$Age[m_M[, t + 1] != "D"]     <- round_age(df_X$Age_cl[m_M[, t + 1] != "D"]) # rounded for use with background mortality and FMC
      
      m_Age[, t + 1] <- df_X$Age
      
      # Update death status in m_D and in df_X to capture the costs of last year of life in the Costs function
      m_D[, t + 1] <- ifelse(m_M[, t + 1] == "D" , 1, 0)
      df_X$death <- ifelse(m_D[, t] == 0 & m_D[, t + 1] == 1, 1, 0) # i.e. only 1 if not dead in previous cycle
      
      #### Calculate costs and QALYs ####
      # Calculate costs per individual during cycle t + 1
      m_C[, t + 1]  <- Costs(l_params_all, m_M[, t + 1], df_X)         
      
      # Calculate QALYs per individual during cycle t + 1
      m_E[, t + 1]  <- Effs(l_params_all, m_M[, t + 1], df_X, t, cl = cl) 
      
      # Update time on treatment (used in Probs)    
      df_X$time <- ifelse( m_L[ , t] !=  m_L[ , t + 1], 1, df_X$time + 1) # If previous line was a different treatment, restart at 1   
      
      # Display simulation progress
      if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
        cat('\r', paste(t/n_t * 100, "% done", sep = " "))
      }
      
    } # close the loop for the time points 
    
    #### Calculate informal care costs based on T2D ####
    # Calculate the time to death per individual
    T2D_data_temp    <- ifelse(m_M == "D", 0, 1) # Replace all 'A' for 1 and all 'D' for 0
    T2D_data         <- unname(rowSums(T2D_data_temp)*cl) # Calculate the sums of every row i.e. the time to death per individual, expressed in years
    m_C_IC_a         <- matrix(nrow = n_i, ncol = (n_t + 1)*cl)
    m_C_IC_a[, 1]    <- rep(0, n_i) # Matrix to capture informal care costs with time to death
    
    # Start loop over patients
    for(i in 1:n_i){ 
      T2D <- T2D_data[i] # determine time to death of individual i
      # Start loop over time to death for individual i
      for(t in 1:T2D){ # loop over start simulation until death of individual i
        
        # Calculate informal care costs for individual i at time t
        p_care_use    <- l_params$params_ic_log_intercept + l_params$params_ic_log_female*df_X$Sex[i] +
          l_params$params_ic_log_age*((df_X$Age_start[i]+t)-70) + 
          l_params$params_ic_log_age2*(((df_X$Age_start[i]+t)-70)*((df_X$Age_start[i]+t)-70)) + l_params$params_ic_log_T2D*(T2D-t)
        v_hour_care   <- l_params$params_ic_hours_intercept + l_params$params_ic_hours_female*df_X$Sex[i] +
          l_params$params_ic_hours_age*((df_X$Age_start[i]+t)-70) + l_params$params_ic_hours_T2D*(T2D-t) #estimates hours per day
        m_C_IC_a[i,t] <- (exp(p_care_use)/(1+exp(p_care_use)))*exp(v_hour_care)*c_IC_hr*365.25  # matrix with annual costs
      }
    }
    
    # Adjust the caregiver costs annual matrix to  cycle length
    for (h in 1:((n_t + 1)*cl)) {
      # Each original column value is split evenly between two new columns
      m_C_IC_cl[, 2*h-1] <- m_C_IC_a[, h]*cl
      m_C_IC_cl[, 2*h]   <- m_C_IC_a[, h]*cl
    }
    
    # Adjust the caregiver matrix for those with and without symptoms
    m_C_IC_cl[m_M == "SAF"] <- m_C_IC_cl[m_M == "SAF"] + c_IC_SAF # Add the AF specific caregiver burden costs
    
    #### Calculate and discount lifetime costs and effects ####
    m_C_IC_cl             <- ifelse(is.na(m_C_IC_cl), 0, m_C_IC_cl) # replace NAs with zero in the cycles where patients are dead
    m_C                   <- m_C + m_C_IC_cl                        # combine regular costs with informal care costs
    m_LY                  <- ifelse(m_M=="D", 0, cl)                # undiscounted life years (i.e. not corrected for quality of life)
    
    tc                    <- m_C %*% (v_dwc * v_wcc)  # total discounted cost per individual
    te                    <- m_E %*% (v_dwe * v_wcc)  # total discounted QALYs per individual 
    tLY                   <- m_LY %*% (v_dwe * v_wcc) # total discounted LYs per individual 
    tLY_undisc            <- m_LY %*% (v_wcc)         # total undiscounted LYs per individual 
    
    tc_hat                <- mean(tc)                 # average discounted cost 
    te_hat                <- mean(te)                 # average discounted QALYs
    tLY_hat               <- mean(tLY)                # average disounted LYs
    tLY_undisc_hat        <- mean(tLY_undisc)         # average undisounted LYs
    
    
    #### Store the results from the simulation in a list ####
    if(PSA == F){
      
      # Determine time on treatment in each line
      m_ToT             <- m_L 
      m_ToT[m_M == "D"] <- NA # if you are dead, remove treatment line
      m_ToT             <- m_ToT[,-1] # remove cycle 0
      
      # Determine proportion on treatment line
      df_PoT <- data.frame(p_L1 = rep(0, n_i),
                           p_L2 = rep(0, n_i),
                           p_L3 = rep(0, n_i),
                           p_L4 = rep(0, n_i),
                           p_L5 = rep(0, n_i),
                           p_L6 = rep(0, n_i), 
                           p_noTrt = rep(0, n_i))
      
      df_PoT[ , 1] <- as.numeric(ifelse(rowSums(m_ToT=="1", na.rm = T)==0, NA, rowSums(m_ToT=="1", na.rm = T)))  
      df_PoT[ , 2] <- as.numeric(ifelse(rowSums(m_ToT=="2", na.rm = T)==0, NA, rowSums(m_ToT=="2", na.rm = T)))
      df_PoT[ , 3] <- as.numeric(ifelse(rowSums(m_ToT=="3", na.rm = T)==0, NA, rowSums(m_ToT=="3", na.rm = T)))
      df_PoT[ , 4] <- as.numeric(ifelse(rowSums(m_ToT=="4", na.rm = T)==0, NA, rowSums(m_ToT=="4", na.rm = T)))
      df_PoT[ , 5] <- as.numeric(ifelse(rowSums(m_ToT=="5", na.rm = T)==0, NA, rowSums(m_ToT=="5", na.rm = T)))
      df_PoT[ , 6] <- as.numeric(ifelse(rowSums(m_ToT=="6", na.rm = T)==0, NA, rowSums(m_ToT=="6", na.rm = T)))
      df_PoT[ , 7] <- as.numeric(ifelse(rowSums(m_ToT=="555", na.rm = T)==0, NA, rowSums(m_ToT=="555", na.rm = T)))
      
      t_L1    <- mean(df_PoT[, 1], na.rm = T)*cl 
      t_L2    <- mean(df_PoT[, 2], na.rm = T)*cl 
      t_L3    <- mean(df_PoT[, 3], na.rm = T)*cl 
      t_L4    <- mean(df_PoT[, 4], na.rm = T)*cl 
      t_L5    <- mean(df_PoT[, 5], na.rm = T)*cl 
      t_L6    <- mean(df_PoT[, 6], na.rm = T)*cl 
      t_noTrt <- mean(df_PoT[, 7], na.rm = T)*cl 
      
      # Proportion receiving treatment
      p_L1    <- 100 #everybody starts on line 1
      p_L2    <- (sum(!is.na(df_PoT[, 2]))/n_i)*100
      p_L3    <- (sum(!is.na(df_PoT[, 3]))/n_i)*100
      p_L4    <- (sum(!is.na(df_PoT[, 4]))/n_i)*100
      p_L5    <- (sum(!is.na(df_PoT[, 5]))/n_i)*100
      p_L6    <- (sum(!is.na(df_PoT[, 6]))/n_i)*100
      p_noTrt <- (sum(!is.na(df_PoT[, 7]))/n_i)*100
      
      results <- list(m_M = m_M, m_C = m_C, m_E = m_E, m_L = m_L, 
                      #m_ToT = m_ToT, m_TH = m_TH, m_AE = m_AE, m_C_IC_cl = m_C_IC_cl,
                      #m_C_IC_a = m_C_IC_a,T2D_data = T2D_data,T2D_data_temp = T2D_data_temp, m_curTrt = m_curTrt, m_D = m_D, m_Age = m_Age,
                      tc = tc, te = te, tLY = tLY, tLY_undisc = tLY_undisc, 
                      tc_hat = tc_hat, te_hat = te_hat, tLY_hat = tLY_hat, tLY_undisc_hat = tLY_undisc_hat, 
                      p_L1 = p_L1, p_L2 = p_L2, p_L3 = p_L3, p_L4 = p_L4, p_L5 = p_L5, p_L6 = p_L6, p_noTrt = p_noTrt, 
                      t_L1= t_L1, t_L2 = t_L2, t_L3 = t_L3, t_L4 = t_L4,  t_L5 = t_L5, t_L6 = t_L6, t_noTrt = t_noTrt)
    }
    
    if(PSA == T){
      results <- list(tc_hat = tc_hat, te_hat = te_hat)
      }
    
    if(mainresults == T){
      results <- list(tc_hat = tc_hat, te_hat = te_hat)
    }
    
    return(results)  # return the results
  }) # end of with(l_params)
  
} # end of the MicroSim function  
library(readr)
library(survival)
library(ggplot2)
library(dplyr)
library(survminer)
library(flexsurv)
library(here)
library(pacman)
library(tidyr)
p_load_gh("DARTH-git/darthtools")
source(here::here("functions", "functions.R"))

rm(list = ls())
options(scipen = 999)
df <- read.csv(here::here("input", "KA_TTE_20250417.csv"), sep = ",", header = TRUE)
colnames(df) <- c("verz_nr_","KA_n", "dag_KA_n", "jr_KA_n", "type_z", "jr_z", "event", "dgn_na_KA_1", "dgn_na_KA_n", "jr_na_KA_n", "sel_90_d")

# Recode KA_n variable to numeric variable
df$KA_n <- as.numeric(ifelse(df$KA_n == "KA_1", 1,
                             ifelse(df$KA_n == "KA_2", 2,
                                    ifelse(df$KA_n == "KA_3", 3,
                                           ifelse(df$KA_n == "KA_4", 4,
                                                  ifelse(df$KA_n == "KA_5", 5, df$KA_n))))))


km_data_combined <- df %>%
  # Remove events that are not end of follow-up within the blanking period
  filter(!(dgn_na_KA_n <= 90 & type_z != "00_eind_fu")) %>%
  # Keep one event per patient ID and KA_n
  group_by(verz_nr_, KA_n) %>%
  arrange(desc(event), dgn_na_KA_n > 90, dgn_na_KA_n, desc(type_z == "1_kath_abl")) %>%
  slice(1) %>%
  ungroup() %>%
  # Group by patient and arrange
  group_by(verz_nr_) %>%
  arrange(verz_nr_, KA_n) %>%
  # Update dag_KA_n based on previous catheter ablation
  mutate(
    dag_KA_n = case_when(
      lag(type_z) == "1_kath_abl" ~ lag(dgn_na_KA_1),
      TRUE ~ dag_KA_n
    ),
    # Assign sequential numbers from 1 to n within each group
    KA_n = row_number(),
    # Ensure KA_n = 1 always has dag_KA_n = 0
    dag_KA_n = case_when(
      KA_n == 1 ~ 0,
      TRUE ~ dag_KA_n
    )
  ) %>%
  ungroup() %>%
  # Calculate time differences
  mutate(
    dgn_na_KA_n = dgn_na_KA_1 - dag_KA_n,
    jr_na_KA_n = dgn_na_KA_n/365.25
  )

# Create survival object with intervention as group
surv_obj_combined <- survfit(Surv(jr_na_KA_n, event) ~ KA_n, data = km_data_combined)

# Create combined Kaplan-Meier plot
ggsurvplot(surv_obj_combined,
           data = km_data_combined,
           title = "Time to AF recurrence",
           xlab = "Years",  
           ylab = "Proportion without AF recurrence",
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.height = 0.30,
           censor = FALSE,
           surv.scale = "percent",
           legend.title = "Catheter ablation",
           legend.labs = c("First", "Second", "Third", "Fourth", "Fifth"),
           pval = FALSE,
           ggtheme = theme_bw())


# Follow-up 
fu_df <- km_data_combined %>%
  group_by(verz_nr_) %>%
  slice_max(order_by = KA_n, n = 1, with_ties = FALSE) %>%
  ungroup()

summary(km_data_combined$dgn_na_KA_1/365.25)
sd(km_data_combined$dgn_na_KA_1/365.25)

# Pool ablation 3-5 due to overlapping K-M curves
km_data_combined_pool3to5 <- km_data_combined
km_data_combined_pool3to5$KA_n <- ifelse(km_data_combined$KA_n>3, 3, km_data_combined$KA_n)

# Create survival object with intervention as group
surv_obj_combined_pool3to5 <- survfit(Surv(jr_na_KA_n, event) ~ KA_n, data = km_data_combined_pool3to5)

# Extract survival probabilities
summary(surv_obj_combined_pool3to5, times = 1)$surv
summary(surv_obj_combined_pool3to5, times = 2)$surv
summary(surv_obj_combined_pool3to5, times = 3)$surv
summary(surv_obj_combined_pool3to5, times = 4)$surv
summary(surv_obj_combined_pool3to5, times = 5)$surv

# 1-year probabilities of recurrence of AF after second and third CA
1-summary(surv_obj_combined_pool3to5, times = 1)$surv

# Confidence interval (let op: lower is upper en vice versa vanwege 1-surv)
1-summary(surv_obj_combined_pool3to5, times = 1)$lower
1-summary(surv_obj_combined_pool3to5, times = 1)$upper

# Create combined Kaplan-Meier plot
ggsurvplot(surv_obj_combined_pool3to5,
           data = km_data_combined_pool3to5,
           title = "Time to AF recurrence",
           xlab = "Years",  
           ylab = "Proportion without AF recurrence",
           conf.int = FALSE,
           risk.table = TRUE,
           tables.height = 0.30,
           censor = FALSE,
           surv.scale = "percent",
           legend.title = "Catheter ablation",
           legend.labs = c("First", "Second", "Third or higher"),
           surv.table.title = "Time to AF recurrence",
           pval = FALSE,
           ggtheme = theme_bw())

# Distribution of events
count_table <- table(km_data_combined_pool3to5$KA_n, km_data_combined_pool3to5$type_z)
prop_table <- prop.table(count_table[,-1], margin = 1)
df_prop <- as.data.frame.matrix(prop_table)
df_prop*100

# Fit parametric functions
times <- seq(0, 90, 1)

km_data_combined <- as.data.frame(km_data_combined)
km_data_combined1 <- km_data_combined[km_data_combined$KA_n == 1,]
km_data_combined2 <- km_data_combined[km_data_combined$KA_n == 2,]
km_data_combined3_4_5 <- km_data_combined[km_data_combined$KA_n == 3|km_data_combined$KA_n == 4|km_data_combined$KA_n == 5,]

#### DATA ANALYSIS WITHOUT FIRST HALF YEAR: used in model ####
times <- seq(0, 100, 0.5) #n_t = 180 half years so 90 might be enough but just to be sure 100 years

km_data_combined <- as.data.frame(km_data_combined)
km_data_excl_1st_6m <- km_data_combined[km_data_combined$dgn_na_KA_n>182.625, ]
km_data_excl_1st_6m$dgn_na_KA_n <- km_data_excl_1st_6m$dgn_na_KA_n-182.625
km_data_excl_1st_6m$jr_na_KA_n <- km_data_excl_1st_6m$dgn_na_KA_n/365.25
km_data_excl_1st_6m$jr_na_KA_n <- as.numeric(km_data_excl_1st_6m$jr_na_KA_n)

km_data_excl_1st_6m_KA1 <- km_data_excl_1st_6m[km_data_excl_1st_6m$KA_n == 1,]
km_data_excl_1st_6m_KA2 <- km_data_excl_1st_6m[km_data_excl_1st_6m$KA_n == 2,]
km_data_excl_1st_6m_KA3_4_5 <- km_data_excl_1st_6m[km_data_excl_1st_6m$KA_n == 3|km_data_excl_1st_6m$KA_n == 4|km_data_excl_1st_6m$KA_n == 5,]

fit_ablation1  <- fit.fun(time = 'jr_na_KA_n', status  = 'event', data = km_data_excl_1st_6m_KA1, extrapolate = TRUE, times = times) 
# Check AIC/BIC of each model to assess goodness-of-fit
GoF_ablation1      <- data.frame(AIC = fit_ablation1$AIC,  BIC = fit_ablation1$BIC) 
# Select best-fitting models
choose_ablation1   <- rownames(GoF_ablation1)[which.min(GoF_ablation1$AIC)]
# Save best-fitting model
best_fit_ablation1 <- fit_ablation1[[choose_ablation1]]

fit_ablation2  <- fit.fun(time = "jr_na_KA_n", status  = "event" , data = km_data_excl_1st_6m_KA2, extrapolate = TRUE, times = times) 
# Check AIC/BIC of each model to assess goodness-of-fit
GoF_ablation2      <- data.frame(AIC = fit_ablation2$AIC,  BIC = fit_ablation2$BIC) 
# Select best-fitting models
choose_ablation2   <- rownames(GoF_ablation2)[which.min(GoF_ablation2$AIC)]
# Save best-fitting model
best_fit_ablation2 <- fit_ablation2[[choose_ablation2]]

fit_ablation3  <- fit.fun(time = "jr_na_KA_n", status  = "event" , data = km_data_excl_1st_6m_KA3_4_5, extrapolate = TRUE, times = times) 
# Check AIC/BIC of each model to assess goodness-of-fit
GoF_ablation3      <- data.frame(AIC = fit_ablation3$AIC,  BIC = fit_ablation3$BIC) 
# Select best-fitting models
choose_ablation3   <- rownames(GoF_ablation3)[which.min(GoF_ablation3$AIC)]
# Save best-fitting model
best_fit_ablation3 <- fit_ablation3[[choose_ablation3]]

save(best_fit_ablation1, best_fit_ablation2, best_fit_ablation3, file = here::here("input", "TTE_AF_Vektis_update.RData"))

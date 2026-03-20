library(dplyr)
library(tidyr)
library(car)
library(nlme)

cd4_data <- read.table("CD4_data_Toread.txt", header = T, sep = "")
colnames(cd4_data) <- c("SubjectID", "Treatment", "Age", "Gender", "Week", "LogCD4")
cd4_data$Treatment <- factor(cd4_data$Treatment)
cd4_data$SubjectID <- factor(cd4_data$SubjectID)

# The raw 'Week' is continuous (e.g., 7.5714). We need to align these to the 
# scheduled visits: 0, 8, 16, 24, 32, 40.
# We round the weeks to the nearest multiple of 8.
cd4_data$ScheduledWeek <- round(cd4_data$Week / 8) * 8

# Filter to keep only the weeks of interest (0 to 40)
cd4_data <- cd4_data %>% filter(ScheduledWeek <= 40)

# Handle duplicates: If a patient has multiple measures in the same bin, 
# typically we take the first one or the one closest to the scheduled date.
# Here we take the first occurrence for simplicity.
cd4_data_clean <- cd4_data %>%
  group_by(SubjectID, ScheduledWeek) %>%
  slice(1) %>%
  ungroup()

# 1. Prepare data (ensure it is in long format)
cd4_long <- cd4_data_clean %>%
  mutate(Treatment = factor(Treatment),
         Gender = factor(Gender),
         Time = ScheduledWeek)


# 2. Fit the General Linear Model using REML
library(nlme)
library(dplyr)
gls_model <- gls(LogCD4 ~ Treatment + Treatment:Time + Age + Gender,
                 data = cd4_long,
                 correlation = corSymm(form = ~ 1 | SubjectID), 
                 weights = varIdent(form = ~ 1 | Time),     
                 method = "REML",
                 na.action = na.omit)

summary(gls_model)$tTable
summary(gls_model)$modelStruct
corMatrix(gls_model$modelStruct$corStruct)[[1]]
getVarCov(gls_model)

library(nlme)
library(dplyr)

cd4_long <- cd4_data_clean %>%
  mutate(Treatment = factor(Treatment),
         Gender = factor(Gender),
         Time_cat = factor(ScheduledWeek),
         Time_cont = ScheduledWeek)       

# Model A: Saturated Model (Profile Analysis)
mod_saturated <- gls(LogCD4 ~ Treatment * Time_cat + Age + Gender,
                 data = cd4_long,
                 correlation = corSymm(form = ~ 1 | SubjectID), 
                 weights = varIdent(form = ~ 1 | Time_cat),     
                 method = "ML", # MUST BE ML FOR MEAN SELECTION
                 na.action = na.omit)

# Model B: Linear Trend Model
mod_linear <- gls(LogCD4 ~ Treatment * Time_cont + Age + Gender,
                 data = cd4_long,
                 correlation = corSymm(form = ~ 1 | SubjectID), 
                 weights = varIdent(form = ~ 1 | Time_cat),     
                 method = "ML", 
                 na.action = na.omit)

# Model C: Quadratic Trend Model
mod_quadratic <- gls(LogCD4 ~ Treatment * (Time_cont + I(Time_cont^2)) + Age + Gender,
                 data = cd4_long,
                 correlation = corSymm(form = ~ 1 | SubjectID), 
                 weights = varIdent(form = ~ 1 | Time_cat),     
                 method = "ML", 
                 na.action = na.omit)

anova(mod_linear, mod_quadratic, mod_saturated)

# Model 1: Unstructured Covariance (UN)
mod_cov_un <- gls(LogCD4 ~ Treatment * Time_cat + Age + Gender,
                 data = cd4_long,
                 correlation = corSymm(form = ~ 1 | SubjectID), 
                 weights = varIdent(form = ~ 1 | Time_cat),     
                 method = "REML", # MUST BE REML FOR COVARIANCE SELECTION
                 na.action = na.omit)

# Model 2: Compound Symmetry (CS)
mod_cov_cs <- gls(LogCD4 ~ Treatment * Time_cat + Age + Gender,
                 data = cd4_long,
                 correlation = corCompSymm(form = ~ 1 | SubjectID), 
                 method = "REML", 
                 na.action = na.omit)

# Model 3: Autoregressive Order 1 (AR1)
mod_cov_ar1 <- gls(LogCD4 ~ Treatment * Time_cat + Age + Gender,
                 data = cd4_long,
                 correlation = corAR1(form = ~ 1 | SubjectID), 
                 method = "REML", 
                 na.action = na.omit)

anova(mod_cov_un, mod_cov_cs, mod_cov_ar1)

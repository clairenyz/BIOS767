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

library(nlme)
library(dplyr)

library(nlme)
library(dplyr)

# Part 3: Final Model
mod_final <- gls(LogCD4 ~ Treatment * Time_cat + Age + Gender,
                 data = cd4_long,
                 correlation = corSymm(form = ~ 1 | SubjectID),
                 weights = varIdent(form = ~ 1 | Time_cat),
                 method = "REML",
                 na.action = na.omit)

# 1. Fixed Effects Table
fixef_table <- as.data.frame(summary(mod_final)$tTable)
fixef_table$Parameter <- rownames(fixef_table)
rownames(fixef_table) <- NULL

fixef_table <- fixef_table %>%
  select(Parameter, Value, Std.Error, `t-value`, `p-value`)

# 95% CI (protected against Hessian failure)
fixef_ci <- tryCatch({
  ci <- as.data.frame(intervals(mod_final, which = "coef")$coef)
  ci$Parameter <- rownames(ci)
  rownames(ci) <- NULL
  ci
}, error = function(e) {
  message("intervals() failed: ", e$message)
  data.frame(Parameter = fixef_table$Parameter,
             lower = NA, est. = NA, upper = NA)
})

fixed_results <- fixef_table %>%
  left_join(fixef_ci, by = "Parameter") %>%
  rename(
    Estimate = Value,
    SE = Std.Error,
    CI_lower = lower,
    CI_upper = upper
  ) %>%
  select(Parameter, Estimate, SE, `t-value`, `p-value`, CI_lower, CI_upper)

cat("\n================ Fixed Effects Results ================\n")
print(fixed_results %>% mutate(across(where(is.numeric), ~ round(., 4))))

# 2. Variance Structure
sigma_hat <- mod_final$sigma

var_weights_coef <- coef(mod_final$modelStruct$varStruct, unconstrained = FALSE)

time_levels <- levels(cd4_long$Time_cat)
ref_time    <- time_levels[1]

rel_sd <- rep(1, length(time_levels))
names(rel_sd) <- time_levels
if (length(var_weights_coef) > 0) {
  rel_sd[names(var_weights_coef)] <- var_weights_coef
}

var_by_time <- data.frame(
  Time        = time_levels,
  Relative_SD = rel_sd,
  SD          = sigma_hat * rel_sd,
  Variance    = (sigma_hat * rel_sd)^2
)

cat("\n================ Estimated Variance by Time Point ================\n")
print(var_by_time %>% mutate(across(where(is.numeric), ~ round(., 4))))

# ============================================================
# 3. Correlation & Covariance Matrix (robust construction)
cor_mat <- corMatrix(mod_final$modelStruct$corStruct)[[1]]

sd_vec  <- sigma_hat * rel_sd                      # SD at each time point
cov_mat <- cor_mat * outer(sd_vec, sd_vec)         # Cov = Cor * SD_i * SD_j
rownames(cov_mat) <- time_levels
colnames(cov_mat) <- time_levels

cat("\n================ Within-Subject Correlation Matrix ================\n")
print(round(cor_mat, 4))

cat("\n================ Estimated Within-Subject Covariance Matrix ================\n")
print(round(cov_mat, 4))

# 4. Upper Triangle (variances on diagonal, covariances above)
cov_df <- as.data.frame(as.table(cov_mat))
colnames(cov_df) <- c("Time1", "Time2", "Value")

cov_df <- cov_df %>%
  mutate(
    i = match(Time1, time_levels),
    j = match(Time2, time_levels)
  ) %>%
  filter(i <= j) %>%                               # diagonal = variances, above = covariances
  mutate(
    Type  = ifelse(i == j, "Variance", "Covariance"),
    Value = round(Value, 4)
  ) %>%
  select(Time1, Time2, Type, Value)

cat("\n================ Variances and Covariances (upper triangle) ================\n")
print(cov_df)

# 5. Final Report Table
fixed_report <- fixed_results %>%
  mutate(
    Estimate  = round(Estimate, 3),
    SE        = round(SE, 3),
    `t-value` = round(`t-value`, 3),
    `p-value` = round(`p-value`, 4),
    CI        = paste0("(", round(CI_lower, 3), ", ", round(CI_upper, 3), ")")
  ) %>%
  select(Parameter, Estimate, SE, CI, `t-value`, `p-value`)

cat("\n================ Fixed Effects Report Table ================\n")
print(fixed_report)

write.csv(fixed_report, "fixed_effect_results.csv",   row.names = FALSE)
write.csv(var_by_time,  "variance_by_time.csv",        row.names = FALSE)
write.csv(cov_df,       "covariance_parameters.csv",   row.names = FALSE)

library(dplyr)
library(ggplot2)
library(nlme)

# Part 4: Visualization
cd4_long <- cd4_long %>%
  mutate(
    Treatment = factor(Treatment),
    Time_cat = factor(Time_cat, levels = c(0, 8, 16, 24, 32, 40)),
    Time_num = as.numeric(as.character(Time_cat))
  )


obs_summary <- cd4_long %>%
  group_by(Treatment, Time_cat, Time_num) %>%
  summarise(
    n = sum(!is.na(LogCD4)),
    mean_logCD4 = mean(LogCD4, na.rm = TRUE),
    sd_logCD4 = sd(LogCD4, na.rm = TRUE),
    se_logCD4 = sd_logCD4 / sqrt(n),
    .groups = "drop"
  )


mean_age <- mean(cd4_long$Age, na.rm = TRUE)
ref_gender <- levels(cd4_long$Gender)[1]

pred_data <- expand.grid(
  Treatment = levels(cd4_long$Treatment),
  Time_cat = levels(cd4_long$Time_cat)
)

pred_data <- pred_data %>%
  mutate(
    Age = mean_age,
    Gender = ref_gender,
    Time_num = as.numeric(as.character(Time_cat))
  )

pred_data$fitted_logCD4 <- predict(mod_final, newdata = pred_data)

p <- ggplot() +
  geom_point(
    data = obs_summary,
    aes(x = Time_num, y = mean_logCD4, shape = Treatment),
    size = 2.8
  ) +
  geom_errorbar(
    data = obs_summary,
    aes(
      x = Time_num,
      ymin = mean_logCD4 - se_logCD4,
      ymax = mean_logCD4 + se_logCD4,
      group = Treatment
    ),
    width = 1.5,
    alpha = 0.6
  ) +
  geom_line(
    data = pred_data,
    aes(x = Time_num, y = fitted_logCD4, linetype = Treatment, group = Treatment),
    linewidth = 1
  ) +
  scale_x_continuous(breaks = c(0, 8, 16, 24, 32, 40)) +
  labs(
    title = "Observed Group Means and Fitted Trajectories of Log CD4 Over Time",
    x = "Week",
    y = "Log CD4",
    shape = "Treatment",
    linetype = "Treatment"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p)

ggsave("part4_cd4_trajectory_plot.png", plot = p, width = 8, height = 5, dpi = 300)

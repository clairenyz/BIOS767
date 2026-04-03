library(nlme)
library(dplyr)

# 1. Load and Prepare Data (Using your standard processing)
cd4_data <- read.table("CD4_data_Toread.txt", header = TRUE, sep = "")
colnames(cd4_data) <- c("SubjectID", "Treatment", "Age", "Gender", "Week", "LogCD4")

cd4_long <- cd4_data %>%
  mutate(Treatment = factor(Treatment),
         Gender = factor(Gender),
         Time = round(Week / 8) * 8) %>%
  filter(Time <= 40) %>%
  group_by(SubjectID, Time) %>%
  dplyr::slice(1) %>%      
  ungroup() %>%
  na.omit()

# ---------------------------------------------------------
# PART 1: FIT CANDIDATE LINEAR MIXED MODELS
# ---------------------------------------------------------

# Model 1: Random Intercept Only
# Assumes subjects start at different baseline CD4 levels, but change at the same rate.
m_intercept <- lme(LogCD4 ~ Treatment * Time + Age + Gender,
                   random = ~ 1 | SubjectID,
                   data = cd4_long,
                   method = "REML")

# Model 2: Random Intercept and Slope
# Assumes subjects start at different baseline levels AND change at different rates over time.
m_slope <- lme(LogCD4 ~ Treatment * Time + Age + Gender,
               random = ~ 1 + Time | SubjectID,
               data = cd4_long,
               method = "REML")

# ---------------------------------------------------------
# PART 2: EVALUATING RANDOM EFFECTS
# ---------------------------------------------------------

# 1. Testing the Random Effects (Likelihood Ratio Test)
# Compare the random intercept model to the random slope model
cat("\n--- Likelihood Ratio Test for Random Effects ---\n")
anova_results <- anova(m_intercept, m_slope)
print(anova_results)

# 2. Extract Fixed Effects and Variance Components for preferred model
cat("\n--- Fixed Effects Estimates (m_slope) ---\n")
print(summary(m_slope)$tTable)

cat("\n--- Estimated Variance Components (D Matrix) ---\n")
print(getVarCov(m_slope))

cat("\n--- Estimated Residual Variance (Within-Subject) ---\n")
print(m_slope$sigma^2)


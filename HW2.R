library(dplyr)
library(tidyr)
library(car)

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

# Reshape data from Long to Wide format
# We need one column per time point (Week0, Week8, etc.)
cd4_wide <- cd4_data_clean %>%
  select(SubjectID, Treatment, ScheduledWeek, LogCD4) %>%
  pivot_wider(names_from = ScheduledWeek, 
              values_from = LogCD4, 
              names_prefix = "Week")

cd4_wide_complete <- na.omit(cd4_wide)
print(paste("Number of subjects with complete data:", nrow(cd4_wide_complete)))
response_matrix <- as.matrix(cd4_wide_complete[, grep("Week", names(cd4_wide_complete))])

### Part 3: MANOVA
lm_model <- lm(response_matrix ~ Treatment, data = cd4_wide_complete)
time_points <- c("0", "8", "16", "24", "32", "40")
idata <- data.frame(Time = factor(time_points, levels = time_points))
manova_results <- Anova(lm_model, idata = idata, idesign = ~Time, type = "III")
print(manova_results)

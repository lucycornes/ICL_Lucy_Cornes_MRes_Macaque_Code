# Load required libraries
library(tidyverse)
library(lme4)        # for mixed effects models
library(lmerTest) 
library(dplyr)
install.packages("tidyverse")
# for p-values in mixed models

# =============================================================================
# DATA PREPARATION AND MERGING (FOLLOWING WORKING APPROACH)
# =============================================================================

# Read the data
maternal_pca <- read.csv("new_maternal_data_with_PCA_and_DavidScore.csv")
male_behavior <- read.csv("male_behaviour_master_sheet2025.csv")
mother_data <- read.csv("mother_master_sheet2025.csv")

# Check data structures
cat("Maternal PCA data:", nrow(maternal_pca), "rows\n")
cat("Male behavior data:", nrow(male_behavior), "rows\n") 
cat("Mother data:", nrow(mother_data), "rows\n")

# Get mothers with maternal care data
mothers_with_care <- unique(maternal_pca$Actor)
cat("Mothers with maternal care data:", length(mothers_with_care), "\n")

# Create the analysis dataset - selecting only needed columns to avoid conflicts
# Step 1: Start with mother-son links, select only needed columns
analysis_data <- mother_data %>%
  dplyr::select(MotherID, MotherID_unique, SonID, BirthSeason, TotalSons) %>%
  filter(!is.na(SonID)) %>%
  filter(MotherID %in% mothers_with_care) %>%
  # Step 2: Add maternal care data
  left_join(maternal_pca, by = c("MotherID" = "Actor")) %>%
  # Step 3: Add male behavior data - select only needed columns to avoid duplicates
  left_join(male_behavior %>% dplyr::select(ID, adjusted_SSB, Scan_Count), 
            by = c("SonID" = "ID")) %>%
  # Step 4: Remove sons without SSB data
  filter(!is.na(adjusted_SSB))

cat("Dataset before scan filtering:", nrow(analysis_data), "mother-son pairs\n")

# Apply scan threshold of 5
threshold <- 5
analysis_data <- analysis_data %>% 
  filter(Scan_Count >= threshold) %>%
  # Calculate proportion of sons included per mother
  group_by(MotherID) %>%
  mutate(
    sons_in_study = n(),
    proportion_sons = sons_in_study / TotalSons
  ) %>%
  ungroup()

cat("Final sample after scan threshold >=", threshold, ":\n")
cat("Sample size:", nrow(analysis_data), "mother-son pairs\n")
cat("Mothers:", length(unique(analysis_data$MotherID)), "\n")
cat("Sons:", length(unique(analysis_data$SonID)), "\n")

# Check if we have enough data for analysis
if (length(unique(analysis_data$MotherID)) < 5) {
  stop("Too few mothers for mixed effects model. Consider lowering threshold.")
}

# Set weights using proportion of sons
analysis_data$weights <- analysis_data$proportion_sons

# Check for missing values in key variables
key_vars <- c("Maternality_Score", "adjusted_SSB", "MotherID", "proportion_sons")
missing_summary <- sapply(analysis_data[key_vars], function(x) sum(is.na(x)))
cat("\nMissing values in key variables:\n")
print(missing_summary)

# Remove any remaining missing values
analysis_data_clean <- analysis_data[complete.cases(analysis_data[key_vars]), ]
cat("Clean dataset rows:", nrow(analysis_data_clean), "\n")

# =============================================================================
# EXPLORATORY DATA ANALYSIS  
# =============================================================================

cat("\n=== EXPLORATORY DATA ANALYSIS ===\n")

# Examine the distribution of adjusted SSB
hist(analysis_data_clean$adjusted_SSB, 
     main = "Distribution of Sons' Adjusted SSB", 
     xlab = "Adjusted SSB", 
     col = "lightblue")

# Summary statistics
cat("Summary of Adjusted SSB:\n")
print(summary(analysis_data_clean$adjusted_SSB))

# Check for normality
if(nrow(analysis_data_clean) >= 3 && nrow(analysis_data_clean) <= 5000) {
  shapiro_test <- shapiro.test(analysis_data_clean$adjusted_SSB)
  cat("Shapiro-Wilk test for SSB normality:\n")
  cat("W =", round(shapiro_test$statistic, 4), ", p =", round(shapiro_test$p.value, 4), "\n")
}

# Plot relationship
plot(analysis_data_clean$Maternality_Score, analysis_data_clean$adjusted_SSB,
     xlab = "Maternality Score", 
     ylab = "Sons' Adjusted SSB",
     main = "Maternality Score vs Sons' Adjusted SSB")

# Check data structure for random effects
mother_counts <- table(analysis_data_clean$MotherID)
cat("Sons per mother distribution:\n")
print(table(mother_counts))

# =============================================================================
# MODEL FITTING - START WITH LINEAR MIXED EFFECTS
# =============================================================================

cat("\n=== FITTING LINEAR MIXED EFFECTS MODEL ===\n")

# Linear mixed effects model with random intercept for mother
model_lmer <- lmer(adjusted_SSB ~ Maternality_Score + (1|MotherID), 
                   data = analysis_data_clean,
                   weights = weights)

summary(model_lmer)

# =============================================================================
# MODEL DIAGNOSTICS FOR FAMILY SELECTION
# =============================================================================

cat("\n=== MODEL DIAGNOSTICS ===\n")

# Check residuals
residuals_model <- residuals(model_lmer)
fitted_values <- fitted(model_lmer)

# Residual plots
par(mfrow = c(2, 2))

# 1. Residuals vs fitted
plot(fitted_values, residuals_model,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# 2. Q-Q plot
qqnorm(residuals_model, main = "Normal Q-Q Plot")
qqline(residuals_model, col = "red")

# 3. Scale-location plot
plot(fitted_values, sqrt(abs(residuals_model)),
     xlab = "Fitted Values", ylab = "√|Residuals|",
     main = "Scale-Location Plot")

# 4. Histogram of residuals
hist(residuals_model, main = "Histogram of Residuals", 
     xlab = "Residuals", col = "lightgray")

par(mfrow = c(1, 1))

# Formal normality test on residuals
shapiro_residuals <- shapiro.test(residuals_model)
cat("Shapiro-Wilk test for residual normality:\n")
cat("W =", round(shapiro_residuals$statistic, 4), ", p =", round(shapiro_residuals$p.value, 4), "\n")

if(shapiro_residuals$p.value < 0.05) {
  cat("Residuals are NOT normally distributed (p < 0.05)\n")
  cat("Consider GLMM with appropriate family\n")
} else {
  cat("Residuals are approximately normal (p > 0.05)\n")
  cat("Linear mixed model is appropriate\n")
}

# =============================================================================
# GLMM ALTERNATIVES (IF RESIDUALS ARE NON-NORMAL)
# =============================================================================

# You'll need the statmod package for Tweedie
library(statmod)
library(glmmTMB)

model_tweedie <- glmmTMB(
  adjusted_SSB ~ Maternality_Score + DavidScore + (1|MotherID),
  weights = weights,
  family = tweedie(),
  data = analysis_data
)

summary(model_tweedie)


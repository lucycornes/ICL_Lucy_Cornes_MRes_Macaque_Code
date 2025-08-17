# Load required libraries
library(tidyverse)
library(lme4)        # for mixed effects models
library(lmerTest)    # for p-values in mixed models
library(glmmTMB)     # for Tweedie models

# =============================================================================
# DATA PREPARATION AND MERGING (FOLLOWING WORKING APPROACH)
# =============================================================================

# Read the data
social_pca <- read.csv("social_data_with_PCA_and_DavidScore.csv")
male_behavior <- read.csv("male_behaviour_master_sheet2025.csv")
mother_data <- read.csv("mother_master_sheet2025.csv")

# Check data structures
cat("Social PCA data:", nrow(social_pca), "rows\n")
cat("Male behavior data:", nrow(male_behavior), "rows\n") 
cat("Mother data:", nrow(mother_data), "rows\n")

# Get mothers with social connectedness data
mothers_with_social <- unique(social_pca$Actor)
cat("Mothers with social connectedness data:", length(mothers_with_social), "\n")

# Create the analysis dataset - selecting only needed columns to avoid conflicts
# Step 1: Start with mother-son links, select only needed columns
analysis_data <- mother_data %>%
  select(MotherID, MotherID_unique, SonID, BirthSeason, TotalSons) %>%
  filter(!is.na(SonID)) %>%
  filter(MotherID %in% mothers_with_social) %>%
  # Step 2: Add social connectedness data
  left_join(social_pca, by = c("MotherID" = "Actor")) %>%
  # Step 3: Add male behavior data - select only needed columns to avoid duplicates
  left_join(male_behavior %>% select(ID, adjusted_SSB, Scan_Count), 
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
key_vars <- c("Social_Connectedness_Score", "adjusted_SSB", "MotherID", "proportion_sons")
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
plot(analysis_data_clean$Social_Connectedness_Score, analysis_data_clean$adjusted_SSB,
     xlab = "Social Connectedness Score", 
     ylab = "Sons' Adjusted SSB",
     main = "Social Connectedness Score vs Sons' Adjusted SSB")

# Check data structure for random effects
mother_counts <- table(analysis_data_clean$MotherID)
cat("Sons per mother distribution:\n")
print(table(mother_counts))

# =============================================================================
# MODEL FITTING - START WITH LINEAR MIXED EFFECTS
# =============================================================================

cat("\n=== FITTING LINEAR MIXED EFFECTS MODEL ===\n")

# Linear mixed effects model with random intercept for mother
model_lmer <- glmmTMB(adjusted_SSB ~ Social_Connectedness_Score + (1|MotherID), 
                   data = analysis_data_clean,
                   weights = weights)

summary(model_lmer)

library(glmmTMB)

# Tweedie model with random intercept for mother
model_tweedie <- glmmTMB(adjusted_SSB ~ Social_Connectedness_Score + (1|MotherID), 
                         data = analysis_data_clean,
                         family = tweedie(link = "log"),
                         weights = weights)

summary(model_tweedie)

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

cat("\n=== FITTING TWEEDIE GLMM ===\n")

# Fit Tweedie GLMM (handles zeros and continuous positive values)
model_tweedie <- glmmTMB(
  adjusted_SSB ~ Social_Connectedness_Score + DavidScore + (1|MotherID),
  weights = weights,
  family = tweedie(),
  data = analysis_data_clean
)

summary(model_tweedie)

# =============================================================================
# MODEL COMPARISON AND SELECTION
# =============================================================================

cat("\n=== MODEL COMPARISON ===\n")

# Compare models using AIC
models_list <- list("Linear" = model_lmer, "Tweedie" = model_tweedie)

# Extract AIC values
aic_values <- sapply(models_list, AIC)
cat("AIC values:\n")
print(sort(aic_values))

# Select best model
best_model_name <- names(which.min(aic_values))
best_model <- models_list[[best_model_name]]

cat("\nBest model based on AIC:", best_model_name, "\n")

# =============================================================================
# FINAL MODEL RESULTS
# =============================================================================

cat("\n=== FINAL MODEL SUMMARY ===\n")
summary(best_model)

# Get results based on model type
if(best_model_name == "Tweedie") {
  # For Tweedie GLMM
  coef_social <- fixef(best_model)$cond["Social_Connectedness_Score"]
  se_social <- sqrt(diag(vcov(best_model)$cond))["Social_Connectedness_Score"]
  z_value <- coef_social / se_social
  p_value <- 2 * (1 - pnorm(abs(z_value)))
  
  cat("\n=== TWEEDIE GLMM RESULTS ===\n")
  cat("Coefficient (log scale):", round(coef_social, 4), "\n")
  cat("Standard Error:", round(se_social, 4), "\n")
  cat("Z-value:", round(z_value, 3), "\n")
  cat("P-value:", round(p_value, 4), "\n")
  
  # Effect size on original scale
  multiplicative_effect <- exp(coef_social)
  cat("Multiplicative effect:", round(multiplicative_effect, 3), "\n")
  cat("For each 1-unit increase in Social Connectedness Score,")
  cat("sons' SSB multiplied by", round(multiplicative_effect, 3), "\n")
  cat("This represents a", round((multiplicative_effect - 1) * 100, 1), "% change\n")
  
} else {
  # For linear models
  coef_social <- fixef(best_model)["Social_Connectedness_Score"]
  se_social <- sqrt(diag(vcov(best_model)))["Social_Connectedness_Score"] 
  t_value <- coef_social / se_social
  df_val <- df.residual(best_model)
  p_value <- 2 * (1 - pt(abs(t_value), df_val))
  
  cat("\n=== LINEAR MODEL RESULTS ===\n")
  cat("Coefficient:", round(coef_social, 4), "\n")
  cat("Standard Error:", round(se_social, 4), "\n")
  cat("t-value:", round(t_value, 3), "\n")
  cat("P-value:", round(p_value, 4), "\n")
  
  cat("For each 1-unit increase in Social Connectedness Score,")
  cat("sons' SSB changes by", round(coef_social, 4), "units\n")
}

# Significance
if(p_value < 0.05) {
  cat("*** SIGNIFICANT EFFECT: Social Connectedness Score affects sons' SSB ***\n")
} else {
  cat("No significant effect of Social Connectedness Score on sons' SSB\n")
}

# Random effects
random_effects <- VarCorr(best_model)
cat("\nRandom Effects (Mother variance):\n")
print(random_effects)

# =============================================================================
# EFFECT SIZE INTERPRETATION
# =============================================================================

cat("\n=== EFFECT SIZE INTERPRETATION ===\n")
cat("Social Connectedness Score coefficient:", round(coef_social, 4), "\n")

# Show range of Social Connectedness Scores for context
cat("\nSocial Connectedness Score range:\n")
cat("Min:", round(min(analysis_data_clean$Social_Connectedness_Score), 2), "\n")
cat("Max:", round(max(analysis_data_clean$Social_Connectedness_Score), 2), "\n")
cat("Range:", round(max(analysis_data_clean$Social_Connectedness_Score) - min(analysis_data_clean$Social_Connectedness_Score), 2), "\n")

# Calculate effect across the observed range
if(best_model_name == "Tweedie") {
  range_effect <- exp(coef_social * (max(analysis_data_clean$Social_Connectedness_Score) - min(analysis_data_clean$Social_Connectedness_Score)))
  cat("Moving from lowest to highest social connectedness would multiply sons' SSB by:", round(range_effect, 3), "\n")
} else {
  range_effect <- coef_social * (max(analysis_data_clean$Social_Connectedness_Score) - min(analysis_data_clean$Social_Connectedness_Score))
  cat("Moving from lowest to highest social connectedness would change sons' SSB by:", round(range_effect, 3), "units\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Final model:", best_model_name, "\n")
cat("Sample size:", nrow(analysis_data_clean), "mother-son pairs\n")
cat("Number of mothers:", length(unique(analysis_data_clean$MotherID)), "\n")
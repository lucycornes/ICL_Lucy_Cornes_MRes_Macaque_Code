# SOCIAL CONNECTEDNESS ~ FEMALE FECUNDITY HYPOTHESIS

# Load required packages
library(readr)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(stringr)

cat("=== LOADING EXISTING SOCIAL CONNECTEDNESS DATA ===\n")

analysis_data <- read_csv("social_connectedness.csv")

#1. DATA PREPARATION
# ================
# Create log_age variable if it doesn't exist
if (!"log_age" %in% names(analysis_data)) {
  if ("age" %in% names(analysis_data)) {
    analysis_data$log_age <- log(analysis_data$age)
    cat("✓ Created log_age variable from age\n")
  } else {
    cat("ERROR: Neither 'log_age' nor 'age' column found in data\n")
    cat("Available columns:", paste(names(analysis_data), collapse = ", "), "\n")
    stop("Required age variable missing")
  }
}

# Check for required columns
required_social_vars <- c("rate_comfort_contact", "rate_proximity", "rate_grooming", 
                          "rate_food_sharing", "rate_muzzle_contact", "rate_initiating_contact", 
                          "rate_passive_body_contact", "rate_approaching", "rate_stationary_proximity", 
                          "rate_teeth_chatter", "rate_reconciliation", "rate_lip_smacking", 
                          "rate_food_sharing_give", "rate_food_sharing_receive")

missing_vars <- setdiff(required_social_vars, names(analysis_data))
if (length(missing_vars) > 0) {
  cat("WARNING: Missing social behavior variables:", paste(missing_vars, collapse = ", "), "\n")
  # Remove missing variables from analysis
  required_social_vars <- intersect(required_social_vars, names(analysis_data))
}

# Remove rows with missing values in key variables
analysis_data <- analysis_data %>%
  filter(!is.na(offspring_count), !is.na(log_age)) %>%
  filter_at(vars(all_of(required_social_vars)), all_vars(!is.na(.)))

cat("Analysis dataset contains", nrow(analysis_data), "females after removing missing values\n")

# 2. EXPLORATORY DATA ANALYSIS
# ============================
cat("Analysis dataset contains", nrow(analysis_data), "females\n")

# Check distribution of offspring count
hist(analysis_data$offspring_count, 
     main = "Distribution of Offspring Count (All Females)", 
     xlab = "Offspring Count")

# Density plot
p1 <- ggplot(analysis_data, aes(x = offspring_count)) +
  geom_density(fill = "lightgreen", alpha = 0.7) +
  labs(title = "Distribution of Offspring Count", 
       x = "Offspring Count", y = "Density") +
  theme_minimal()
print(p1)

# Test for normality
if (nrow(analysis_data) <= 5000) {  # Shapiro-Wilk has sample size limits
  shapiro_test <- shapiro.test(analysis_data$offspring_count)
  cat("Shapiro-Wilk p-value:", round(shapiro_test$p.value, 4), "\n")
} else {
  cat("Sample too large for Shapiro-Wilk test\n")
}

# Check for overdispersion with Poisson
poisson_test <- glmmTMB(offspring_count ~ rate_comfort_contact + offset(log_age), 
                        family = poisson(), 
                        data = analysis_data)

# Function to check overdispersion
check_overdispersion <- function(model) {
  rp <- residuals(model, type = "pearson")
  rdf <- df.residual(model)
  dispersion <- sum(rp^2) / rdf
  
  cat("Dispersion ratio:", round(dispersion, 2), "\n")
  
  if (dispersion < 1.5) {
    cat("No strong evidence of overdispersion - Poisson appropriate\n")
  } else if (dispersion < 2) {
    cat("Some overdispersion — consider Negative Binomial\n")
  } else {
    cat("Strong overdispersion — use Negative Binomial or Gaussian\n")
  }
  
  return(dispersion)
}

cat("Testing dispersion with Poisson model:\n")
dispersion_result <- check_overdispersion(poisson_test)

# 3. CHOOSE APPROPRIATE MODEL FAMILY
# ==================================
# Based on dispersion test results, choose appropriate family
if (dispersion_result < 1.5) {
  model_family <- poisson()
  cat("Using Poisson family\n")
} else if (dispersion_result < 2) {
  model_family <- nbinom2()  # Negative binomial
  cat("Using Negative Binomial family\n")
} else {
  model_family <- gaussian()
  cat("Using Gaussian family due to overdispersion\n")
}

# 4. SOCIAL CONNECTEDNESS MODEL ANALYSIS
# ======================================

# Full model with available social connectedness categories
formula_full <- as.formula(paste("offspring_count ~", 
                                 paste(required_social_vars, collapse = " + "), 
                                 "+ offset(log_age)"))

full_social_model <- glm(formula_full, family = model_family, data = analysis_data)

summary(full_social_model)

# 5. STEPWISE MODEL SELECTION
# ============================

# Function to perform stepwise elimination for glmmTMB using AIC
stepwise_elimination <- function(data, variables, family_to_use) {
  current_vars <- variables
  
  while(length(current_vars) > 1) {
    cat("\n=== Testing removal of variables ===\n")
    
    aic_values <- numeric()
    
    # Fit current full model
    formula_str_full <- paste("offspring_count ~", paste(current_vars, collapse = " + "), "+ offset(log_age)")
    full_model <- glmmTMB(as.formula(formula_str_full), family = family_to_use, data = data)
    full_aic <- AIC(full_model)
    
    for(var in current_vars) {
      # Create reduced model without this variable
      reduced_vars <- setdiff(current_vars, var)
      formula_str <- paste("offspring_count ~", paste(reduced_vars, collapse = " + "), "+ offset(log_age)")
      
      reduced_model <- glmmTMB(as.formula(formula_str), family = family_to_use, data = data)
      reduced_aic <- AIC(reduced_model)
      
      # Store AIC difference (positive = worse model when removing variable)
      aic_values[var] <- reduced_aic - full_aic
      
      cat("Remove", var, "- AIC change:", round(aic_values[var], 2), "\n")
    }
    
    # Find variable with smallest AIC increase (least important)
    min_aic_var <- names(which.min(aic_values))
    min_aic_change <- min(aic_values)
    
    # If AIC increases by less than 2 when removing variable, remove it
    # (AIC difference < 2 indicates no meaningful difference)
    if(min_aic_change < 2) {
      cat("Removing", min_aic_var, "(AIC change =", round(min_aic_change, 2), ")\n")
      current_vars <- setdiff(current_vars, min_aic_var)
    } else {
      cat("All remaining variables improve model substantially - stopping\n")
      break
    }
  }
  
  # Fit final model
  formula_str_final <- paste("offspring_count ~", paste(current_vars, collapse = " + "), "+ offset(log_age)")
  final_model <- glmmTMB(as.formula(formula_str_final), family = family_to_use, data = data)
  
  return(list(final_model = final_model, final_vars = current_vars))
}

# Perform stepwise elimination with available social connectedness categories
stepwise_result <- stepwise_elimination(analysis_data, required_social_vars, model_family)
final_social_model <- stepwise_result$final_model

cat("\n=== FINAL SOCIAL CONNECTEDNESS MODEL ===\n")
summary(final_social_model)

# 6. MODEL VALIDATION
# ====================

# Check model with DHARMa
if(require(DHARMa, quietly = TRUE)) {
  residuals_dharma <- simulateResiduals(final_social_model)
  plot(residuals_dharma, main = "DHARMa Residual Diagnostics")
} else {
  cat("Install DHARMa package for residual diagnostics: install.packages('DHARMa')\n")
}

# 7. MODEL VISUALIZATION
# ======================

# Get significant predictors
significant_vars <- stepwise_result$final_vars

# Create plots for significant predictors
for(var in significant_vars) {
  # Create prediction data
  pred_data <- data.frame(
    x = seq(min(analysis_data[[var]], na.rm = TRUE), 
            max(analysis_data[[var]], na.rm = TRUE), length.out = 100)
  )
  names(pred_data)[1] <- var
  
  # Add other variables at their means
  for(other_var in setdiff(significant_vars, var)) {
    pred_data[[other_var]] <- mean(analysis_data[[other_var]], na.rm = TRUE)
  }
  pred_data$log_age <- mean(analysis_data$log_age, na.rm = TRUE)
  
  # Get predictions
  pred_data$predicted <- predict(final_social_model, newdata = pred_data, type = "response")
  
  # Create plot using modern ggplot2 syntax
  behavior_name <- gsub("rate_", "", var)
  behavior_name <- gsub("_", " ", behavior_name)
  behavior_name <- str_to_title(behavior_name)
  
  p <- ggplot() +
    geom_point(data = analysis_data, 
               aes(x = .data[[var]], y = offspring_count), 
               alpha = 0.6, color = "gray50") +
    geom_line(data = pred_data, 
              aes(x = .data[[var]], y = predicted), 
              color = "darkgreen", linewidth = 1) +
    labs(x = paste("Rate of", behavior_name, "(events/hour)"),
         y = "Offspring Count",
         title = paste("Effect of", behavior_name, "on Offspring Count")) +
    theme_minimal()
  
  print(p)
}

# 8. SAVE RESULTS
# ================

# Create summary table of results
if(length(significant_vars) > 0) {
  model_summary <- summary(final_social_model)
  coef_table <- model_summary$coefficients$cond
  results_table <- data.frame(
    Variable = rownames(coef_table),
    Estimate = round(coef_table[,1], 4),
    SE = round(coef_table[,2], 4),
    Z_value = round(coef_table[,3], 4),
    P_value = round(coef_table[,4], 4),
    Significant = ifelse(coef_table[,4] < 0.05, "Yes", "No")
  )
  write.csv(results_table, "social_connectedness_model_results.csv", row.names = FALSE)
  
  # Print results summary
  cat("\n=== MODEL RESULTS SUMMARY ===\n")
  print(results_table)
  
  # Print model diagnostics
  cat("\n=== MODEL DIAGNOSTICS ===\n")
  cat("AIC:", AIC(final_social_model), "\n")
  cat("Number of observations:", nrow(analysis_data), "\n")
  cat("Significant predictors:", paste(significant_vars, collapse = ", "), "\n")
}

# Save final dataset with predictions
analysis_data$predicted_offspring <- predict(final_social_model, type = "response")
write.csv(analysis_data, "social_connectedness_analysis_complete.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to: social_connectedness_model_results.csv\n")
cat("Complete dataset saved to: social_connectedness_analysis_complete.csv\n")
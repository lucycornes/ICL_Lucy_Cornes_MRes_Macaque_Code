# Load required libraries
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)  # for VIF testing
library(corrplot)

# ==============================================================================
# 1. LOAD DATA AND PREPARE (using your existing code structure)
# ==============================================================================


# Load the existing maternal dataset
maternal_data <- read_csv("maternal_behaviour_rate.csv")
male_behavior <- read_csv("male_behaviour_master_sheet2025.csv") 
mother_data <- read_csv("mother_master_sheet2025.csv")

cat("Maternal behavior data:", nrow(maternal_data), "mothers\n")
cat("Male behavior data:", nrow(male_behavior), "rows\n")
cat("Mother data:", nrow(mother_data), "rows\n")

# ==============================================================================
# 2. IDENTIFY MATERNAL BEHAVIORS WITH SUFFICIENT VARIATION
# ==============================================================================

cat("\n=== IDENTIFYING MATERNAL BEHAVIORS WITH SUFFICIENT VARIATION ===\n")

# Get all rate columns representing maternal behaviors
maternal_columns <- names(maternal_data)[grepl("^rate_", names(maternal_data))]
cat("Found", length(maternal_columns), "maternal behavior rate columns:\n")
cat(paste(maternal_columns, collapse = ", "), "\n")

# Check which behaviors have variation and sufficient data
behavior_summary <- maternal_data %>%
  select(all_of(maternal_columns)) %>%
  summarise(across(everything(), list(
    mean = ~mean(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE),
    zeros = ~sum(. == 0, na.rm = TRUE),
    nonzeros = ~sum(. > 0, na.rm = TRUE),
    max = ~max(., na.rm = TRUE),
    min = ~min(., na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("behavior", "statistic"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(has_variation = sd > 0 & nonzeros >= 3)  # Need at least 3 non-zero observations

cat("\nMaternal behavioral summary:\n")
print(behavior_summary)

# Keep only behaviors with variation
good_behaviors <- behavior_summary$behavior[behavior_summary$has_variation]
cat("\nBehaviors with sufficient variation (", length(good_behaviors), "):\n")
cat(paste(good_behaviors, collapse = ", "), "\n")

if (length(good_behaviors) == 0) {
  stop("No maternal behaviors have sufficient variation for analysis")
}

# ==============================================================================
# 3. LINK MATERNAL DATA TO SONS' SSB DATA
# ==============================================================================

cat("\n=== LINKING MATERNAL DATA TO SONS' SSB DATA ===\n")

# Get mothers with maternal behavior data
mothers_with_data <- unique(maternal_data$Actor)
cat("Mothers with maternal behavior data:", length(mothers_with_data), "\n")

# Create mother-son links and add behavioral data
analysis_data <- mother_data %>%
  select(MotherID, MotherID_unique, SonID, BirthSeason, TotalSons) %>%
  filter(!is.na(SonID)) %>%
  filter(MotherID %in% mothers_with_data) %>%
  left_join(maternal_data, by = c("MotherID" = "Actor")) %>%
  left_join(male_behavior %>% select(ID, adjusted_SSB, Scan_Count), by = c("SonID" = "ID")) %>%
  filter(!is.na(adjusted_SSB)) 

cat("Dataset before scan filtering:", nrow(analysis_data), "mother-son pairs\n")

# Apply scan threshold
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

# Check if we have enough variation for analysis
if (length(unique(analysis_data$MotherID)) < 5) {
  stop("Too few mothers for mixed effects model. Consider lowering threshold.")
}

# ==============================================================================
# 4. STANDARDIZE MATERNAL BEHAVIOR VARIABLES
# ==============================================================================

cat("\n=== STANDARDIZING MATERNAL BEHAVIOR VARIABLES ===\n")

# Create standardized versions of all good behavior variables
good_behaviors_in_data <- intersect(good_behaviors, names(analysis_data))

cat("Standardizing", length(good_behaviors_in_data), "maternal behavior variables:\n")
cat(paste(good_behaviors_in_data, collapse = ", "), "\n")

# Before standardization - show original scales
cat("\nOriginal scales (first 3 behaviors):\n")
for(i in 1:min(3, length(good_behaviors_in_data))) {
  col <- good_behaviors_in_data[i]
  if(col %in% names(analysis_data)) {
    values <- analysis_data[[col]]
    cat(col, ": Range =", round(min(values, na.rm = TRUE), 4), "to", 
        round(max(values, na.rm = TRUE), 4), 
        "| Mean =", round(mean(values, na.rm = TRUE), 4),
        "| SD =", round(sd(values, na.rm = TRUE), 4), "\n")
  }
}

# Standardize all good behavior columns
analysis_data <- analysis_data %>%
  mutate(across(all_of(good_behaviors_in_data), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_std"))

# Get list of standardized column names
standardized_behaviors <- paste0(good_behaviors_in_data, "_std")

cat("\nAfter standardization (first 3 behaviors):\n")
for(i in 1:min(3, length(standardized_behaviors))) {
  col <- standardized_behaviors[i]
  if(col %in% names(analysis_data)) {
    values <- analysis_data[[col]]
    cat(col, ": Range =", round(min(values, na.rm = TRUE), 2), "to", 
        round(max(values, na.rm = TRUE), 2), 
        "| Mean =", round(mean(values, na.rm = TRUE), 2),
        "| SD =", round(sd(values, na.rm = TRUE), 2), "\n")
  }
}

cat("\n*** All maternal behavior variables now standardized (mean=0, SD=1) ***\n")

# ==============================================================================
# 5. CHECK FOR MULTICOLLINEARITY
# ==============================================================================

cat("\n=== CHECKING FOR MULTICOLLINEARITY ===\n")

# Create correlation matrix of standardized behaviors
correlation_matrix <- cor(analysis_data[standardized_behaviors], use = "complete.obs")

cat("Correlation matrix of standardized maternal behaviors:\n")
print(round(correlation_matrix, 3))

# Visualize correlation matrix
corrplot(correlation_matrix, method = "color", type = "upper", 
         order = "hclust", tl.cex = 0.7, tl.col = "black",
         title = "Correlation Matrix: Standardized Maternal Behaviors",
         mar = c(0,0,1,0))

# Check for high correlations (> 0.7 or < -0.7)
high_corr_pairs <- which(abs(correlation_matrix) > 0.7 & correlation_matrix != 1, arr.ind = TRUE)
if(nrow(high_corr_pairs) > 0) {
  cat("\nHigh correlations detected (|r| > 0.7):\n")
  for(i in 1:nrow(high_corr_pairs)) {
    row_idx <- high_corr_pairs[i, 1]
    col_idx <- high_corr_pairs[i, 2]
    var1 <- rownames(correlation_matrix)[row_idx]
    var2 <- colnames(correlation_matrix)[col_idx]
    corr_val <- correlation_matrix[row_idx, col_idx]
    cat(paste(var1, "vs", var2, ":", round(corr_val, 3), "\n"))
  }
  cat("Consider removing highly correlated variables to avoid multicollinearity\n")
} else {
  cat("\nNo problematic correlations detected (all |r| <= 0.7)\n")
}

# ==============================================================================
# 6. FIT MULTIVARIATE MODEL WITH ALL BEHAVIORS
# ==============================================================================

cat("\n=== FITTING MULTIVARIATE MODEL WITH ALL MATERNAL BEHAVIORS ===\n")

# Create formula for multivariate model
behavior_formula <- paste(standardized_behaviors, collapse = " + ")
full_formula <- as.formula(paste("adjusted_SSB ~", behavior_formula, "+ (1|MotherID)"))

cat("Model formula:\n")
cat(deparse(full_formula), "\n")

cat("\nNumber of predictors:", length(standardized_behaviors), "\n")
cat("Sample size:", nrow(analysis_data), "\n")
cat("Observations per predictor:", round(nrow(analysis_data) / length(standardized_behaviors), 1), "\n")

# Check if we have enough observations (rule of thumb: 10-15 observations per predictor)
if(nrow(analysis_data) / length(standardized_behaviors) < 10) {
  cat("WARNING: Low observations-to-predictor ratio. Consider model selection.\n")
}

# Fit the multivariate Tweedie model
cat("\nFitting multivariate Tweedie model...\n")

multivariate_model <- NULL
tryCatch({
  multivariate_model <- glmmTMB(
    formula = full_formula,
    weights = proportion_sons,
    family = tweedie(),
    data = analysis_data
  )
  cat("✓ Multivariate Tweedie model fitted successfully\n")
}, error = function(e) {
  cat("✗ Multivariate model failed:", e$message, "\n")
})

# ==============================================================================
# 7. MODEL DIAGNOSTICS AND RESULTS
# ==============================================================================

if (!is.null(multivariate_model)) {
  cat("\n=== MULTIVARIATE MODEL RESULTS ===\n")
  
  # Model summary
  model_summary <- summary(multivariate_model)
  print(model_summary)
  
  # Extract coefficients for maternal behaviors
  coef_table <- model_summary$coefficients$cond
  behavior_coefs <- coef_table[standardized_behaviors, , drop = FALSE]
  
  cat("\n=== MATERNAL BEHAVIOR EFFECTS (CONTROLLING FOR ALL OTHER BEHAVIORS) ===\n")
  cat("*** Effects shown as log-scale coefficients for standardized predictors ***\n")
  cat("*** Each effect controls for all other maternal behaviors ***\n")
  
  # Create results table
  results_df <- data.frame(
    Behavior = gsub("_std$", "", rownames(behavior_coefs)),
    Estimate = behavior_coefs[, "Estimate"],
    SE = behavior_coefs[, "Std. Error"],
    Z_value = behavior_coefs[, "z value"],
    P_value = behavior_coefs[, "Pr(>|z|)"],
    Significant = behavior_coefs[, "Pr(>|z|)"] < 0.05
  )
  
  # Calculate percentage changes for interpretation
  results_df$Percent_Change <- (exp(results_df$Estimate) - 1) * 100
  
  # Clean behavior names
  results_df$Behavior_Clean <- gsub("rate_", "", results_df$Behavior)
  results_df$Behavior_Clean <- gsub("_", " ", results_df$Behavior_Clean)
  results_df$Behavior_Clean <- tools::toTitleCase(results_df$Behavior_Clean)
  
  # Sort by p-value
  results_df <- results_df[order(results_df$P_value), ]
  
  print(results_df[, c("Behavior_Clean", "Estimate", "SE", "Percent_Change", "P_value", "Significant")])
  
  # Highlight significant effects
  cat("\n=== SIGNIFICANT MATERNAL BEHAVIOR EFFECTS (p < 0.05) ===\n")
  significant_effects <- results_df[results_df$Significant, ]
  
  if(nrow(significant_effects) > 0) {
    for(i in 1:nrow(significant_effects)) {
      row <- significant_effects[i, ]
      direction <- ifelse(row$Percent_Change > 0, "INCREASES", "DECREASES")
      cat(paste0(i, ". ", row$Behavior_Clean, ": ", 
                 round(row$Percent_Change, 1), "% change per 1 SD increase (p = ", 
                 round(row$P_value, 4), ") - ", direction, " sons' SSB\n"))
    }
  } else {
    cat("No significant maternal behavior effects found in multivariate model\n")
  }
  
  # ==============================================================================
  # 8. MODEL DIAGNOSTICS
  # ==============================================================================
  
  cat("\n=== MODEL DIAGNOSTICS ===\n")
  
  # Basic model information
  cat("AIC:", round(AIC(multivariate_model), 2), "\n")
  cat("Number of observations:", nobs(multivariate_model), "\n")
  cat("Number of groups (mothers):", length(unique(analysis_data$MotherID)), "\n")
  
  # Residual diagnostics
  tryCatch({
    cat("\nRunning residual diagnostics...\n")
    sim_residuals <- simulateResiduals(multivariate_model)
    
    # Basic diagnostic tests
    disp_test <- testDispersion(sim_residuals)
    uniform_test <- testUniformity(sim_residuals)
    zero_test <- testZeroInflation(sim_residuals)
    
    cat("Dispersion test p-value:", round(disp_test$p.value, 4), "\n")
    cat("Uniformity test p-value:", round(uniform_test$p.value, 4), "\n")
    cat("Zero-inflation test p-value:", round(zero_test$p.value, 4), "\n")
    
    if(disp_test$p.value > 0.05 & uniform_test$p.value > 0.05 & zero_test$p.value > 0.05) {
      cat("✓ Model diagnostics look good!\n")
    } else {
      cat("⚠ Some diagnostic issues detected\n")
    }
    
    # Plot diagnostics
    plot(sim_residuals, main = "Multivariate Maternal Behavior Model Diagnostic Plot")
    
  }, error = function(e) {
    cat("Diagnostic testing failed:", e$message, "\n")
  })
  
  # ==============================================================================
  # 9. VARIANCE INFLATION FACTORS (VIF) TEST
  # ==============================================================================
  
  cat("\n=== CHECKING FOR MULTICOLLINEARITY (VIF) ===\n")
  
  tryCatch({
    # Extract fixed effects for VIF calculation
    # Note: VIF calculation may not work directly with glmmTMB, so we'll use a workaround
    
    # Fit equivalent model with lm for VIF testing
    vif_model <- lm(adjusted_SSB ~ . - MotherID - proportion_sons, 
                    data = analysis_data[, c("adjusted_SSB", standardized_behaviors, "MotherID", "proportion_sons")])
    
    vif_values <- vif(vif_model)
    cat("Variance Inflation Factors (VIF):\n")
    print(round(vif_values, 2))
    
    high_vif <- vif_values[vif_values > 5]
    if(length(high_vif) > 0) {
      cat("\nVariables with high VIF (> 5):\n")
      print(round(high_vif, 2))
      cat("Consider removing these variables due to multicollinearity\n")
    } else {
      cat("\nNo problematic multicollinearity detected (all VIF <= 5)\n")
    }
    
  }, error = function(e) {
    cat("VIF calculation failed:", e$message, "\n")
  })
  
  # ==============================================================================
  # 10. BACKWARDS STEPWISE ELIMINATION FOR MODEL SELECTION
  # ==============================================================================
  
  cat("\n=== BACKWARDS STEPWISE ELIMINATION ===\n")
  
  # Function to perform backwards stepwise elimination for glmmTMB
  backwards_stepwise_maternal <- function(full_model, data, alpha = 0.05) {
    
    cat("Starting backwards stepwise elimination with alpha =", alpha, "\n")
    cat("Full model AIC:", round(AIC(full_model), 2), "\n\n")
    
    # Get initial predictors (exclude intercept and random effects)
    current_predictors <- standardized_behaviors
    current_model <- full_model
    current_aic <- AIC(current_model)
    
    step_counter <- 0
    elimination_log <- data.frame(
      Step = integer(),
      Removed_Variable = character(),
      P_value = numeric(),
      AIC_before = numeric(),
      AIC_after = numeric(),
      Delta_AIC = numeric(),
      Action = character(),
      stringsAsFactors = FALSE
    )
    
    while(length(current_predictors) > 1) {
      step_counter <- step_counter + 1
      cat("=== STEP", step_counter, "===\n")
      cat("Current predictors:", length(current_predictors), "\n")
      
      # Get p-values for all current predictors
      model_summary <- summary(current_model)
      coef_table <- model_summary$coefficients$cond
      
      # Extract p-values for current predictors
      predictor_pvals <- coef_table[current_predictors, "Pr(>|z|)", drop = FALSE]
      
      # Find the predictor with highest p-value
      max_pval_idx <- which.max(predictor_pvals)
      max_pval_var <- rownames(predictor_pvals)[max_pval_idx]
      max_pval <- predictor_pvals[max_pval_idx, 1]
      
      cat("Highest p-value:", round(max_pval, 4), "for", gsub("_std", "", max_pval_var), "\n")
      
      # Check if we should remove this variable
      if(max_pval > alpha) {
        cat("P-value >", alpha, "- attempting removal...\n")
        
        # Create new predictor list without this variable
        new_predictors <- current_predictors[current_predictors != max_pval_var]
        
        # Refit model without this predictor
        if(length(new_predictors) > 0) {
          new_formula <- as.formula(paste("adjusted_SSB ~", 
                                          paste(new_predictors, collapse = " + "), 
                                          "+ (1|MotherID)"))
          
          tryCatch({
            new_model <- glmmTMB(
              formula = new_formula,
              weights = proportion_sons,
              family = tweedie(),
              data = analysis_data  # Use analysis_data directly
            )
            
            new_aic <- AIC(new_model)
            delta_aic <- new_aic - current_aic
            
            cat("AIC change:", round(delta_aic, 2), "\n")
            
            # Log this step
            elimination_log <- rbind(elimination_log, data.frame(
              Step = step_counter,
              Removed_Variable = gsub("_std", "", max_pval_var),
              P_value = max_pval,
              AIC_before = current_aic,
              AIC_after = new_aic,
              Delta_AIC = delta_aic,
              Action = "REMOVED"
            ))
            
            # Update current model
            current_model <- new_model
            current_predictors <- new_predictors
            current_aic <- new_aic
            
            cat("✓ Variable removed. New AIC:", round(current_aic, 2), "\n")
            
          }, error = function(e) {
            cat("✗ Failed to fit model without", max_pval_var, ":", e$message, "\n")
            
            # Log failed attempt
            elimination_log <<- rbind(elimination_log, data.frame(
              Step = step_counter,
              Removed_Variable = gsub("_std", "", max_pval_var),
              P_value = max_pval,
              AIC_before = current_aic,
              AIC_after = NA,
              Delta_AIC = NA,
              Action = "FAILED_REMOVAL"
            ))
            
            break  # Exit if model fitting fails
          })
          
        } else {
          cat("Cannot remove - would result in no predictors\n")
          break
        }
        
      } else {
        cat("All remaining variables significant (p <=", alpha, ") - stopping\n")
        
        # Log final state
        elimination_log <- rbind(elimination_log, data.frame(
          Step = step_counter,
          Removed_Variable = "NONE",
          P_value = max_pval,
          AIC_before = current_aic,
          AIC_after = current_aic,
          Delta_AIC = 0,
          Action = "STOPPED"
        ))
        
        break
      }
      
      cat("\n")
    }
    
    return(list(
      final_model = current_model,
      final_predictors = current_predictors,
      elimination_log = elimination_log,
      final_aic = current_aic
    ))
  }
  
  # Perform backwards stepwise elimination
  stepwise_results_maternal <- backwards_stepwise_maternal(multivariate_model, analysis_data, alpha = 0.05)
  
  # ==============================================================================
  # STEPWISE ELIMINATION RESULTS
  # ==============================================================================
  
  cat("\n=== STEPWISE ELIMINATION SUMMARY ===\n")
  
  # Show elimination log
  cat("Elimination steps:\n")
  print(stepwise_results_maternal$elimination_log)
  
  # Compare models
  cat("\n=== MODEL COMPARISON ===\n")
  cat("Full model AIC:", round(AIC(multivariate_model), 2), "\n")
  cat("Final model AIC:", round(stepwise_results_maternal$final_aic, 2), "\n")
  cat("AIC improvement:", round(AIC(multivariate_model) - stepwise_results_maternal$final_aic, 2), "\n")
  
  # Show final model
  final_maternal_model <- stepwise_results_maternal$final_model
  final_maternal_predictors <- stepwise_results_maternal$final_predictors
  
  cat("\nFinal model predictors (", length(final_maternal_predictors), "):\n")
  cat(paste(gsub("_std", "", final_maternal_predictors), collapse = ", "), "\n")
  
  cat("\n=== FINAL STEPWISE MATERNAL MODEL RESULTS ===\n")
  final_maternal_summary <- summary(final_maternal_model)
  print(final_maternal_summary)
  
  # Extract final model coefficients
  final_maternal_coef_table <- final_maternal_summary$coefficients$cond
  final_maternal_behavior_coefs <- final_maternal_coef_table[final_maternal_predictors, , drop = FALSE]
  
  cat("\n=== FINAL MODEL: MATERNAL BEHAVIOR EFFECTS ===\n")
  
  # Create final results table
  final_maternal_results_df <- data.frame(
    Behavior = gsub("_std$", "", rownames(final_maternal_behavior_coefs)),
    Estimate = final_maternal_behavior_coefs[, "Estimate"],
    SE = final_maternal_behavior_coefs[, "Std. Error"],
    Z_value = final_maternal_behavior_coefs[, "z value"],
    P_value = final_maternal_behavior_coefs[, "Pr(>|z|)"],
    Significant = final_maternal_behavior_coefs[, "Pr(>|z|)"] < 0.05
  )
  
  # Calculate percentage changes
  final_maternal_results_df$Percent_Change <- (exp(final_maternal_results_df$Estimate) - 1) * 100
  
  # Clean behavior names
  final_maternal_results_df$Behavior_Clean <- gsub("rate_", "", final_maternal_results_df$Behavior)
  final_maternal_results_df$Behavior_Clean <- gsub("_", " ", final_maternal_results_df$Behavior_Clean)
  final_maternal_results_df$Behavior_Clean <- tools::toTitleCase(final_maternal_results_df$Behavior_Clean)
  
  # Sort by p-value
  final_maternal_results_df <- final_maternal_results_df[order(final_maternal_results_df$P_value), ]
  
  print(final_maternal_results_df[, c("Behavior_Clean", "Estimate", "SE", "Percent_Change", "P_value", "Significant")])
  
  
  # ==============================================================================
  # FINAL DIAGNOSTICS FOR STEPWISE MODEL
  # ==============================================================================
  
  cat("\n=== FINAL STEPWISE MATERNAL MODEL DIAGNOSTICS ===\n")
  
  tryCatch({
    final_maternal_sim_residuals <- simulateResiduals(final_maternal_model)
    
    final_maternal_disp_test <- testDispersion(final_maternal_sim_residuals)
    final_maternal_uniform_test <- testUniformity(final_maternal_sim_residuals)
    final_maternal_zero_test <- testZeroInflation(final_maternal_sim_residuals)
    
    cat("Final maternal model diagnostics:\n")
    cat("Dispersion test p-value:", round(final_maternal_disp_test$p.value, 4), "\n")
    cat("Uniformity test p-value:", round(final_maternal_uniform_test$p.value, 4), "\n")
    cat("Zero-inflation test p-value:", round(final_maternal_zero_test$p.value, 4), "\n")
    
    if(final_maternal_disp_test$p.value > 0.05 & final_maternal_uniform_test$p.value > 0.05 & final_maternal_zero_test$p.value > 0.05) {
      cat("✓ Final maternal model diagnostics look good!\n")
    } else {
      cat("⚠ Some diagnostic issues detected in final maternal model\n")
    }
    
    # Plot final model diagnostics
    plot(final_maternal_sim_residuals, main = "Final Stepwise Maternal Model Diagnostic Plot")
    
  }, error = function(e) {
    cat("Final maternal model diagnostic testing failed:", e$message, "\n")
  })
  
  # ==============================================================================
  # 11. SUMMARY AND INTERPRETATION
  # ==============================================================================
  
  cat("\n=== MULTIVARIATE MATERNAL BEHAVIOR MODEL INTERPRETATION ===\n")
  cat("This model tests the independent effects of each maternal behavior\n")
  cat("while controlling for all other maternal behaviors simultaneously.\n")
  cat("This approach accounts for the intercorrelated nature of maternal behaviors\n")
  cat("and identifies which behaviors have unique associations with sons' SSB.\n\n")
  
  cat("Key advantages of this multivariate approach:\n")
  cat("1. Controls for correlations between different maternal behaviors\n")
  cat("2. Identifies independent effects of each behavior\n")
  cat("3. Reduces risk of false positives from multiple testing\n")
  cat("4. Provides a more comprehensive view of maternal effects on sons' SSB\n\n")
  
  cat("Results should be interpreted as:\n")
  cat("'The effect of a 1 standard deviation increase in [maternal behavior],\n")
  cat("holding all other maternal behaviors constant'\n")
  
  cat("\n=== STEPWISE ELIMINATION COMPLETE ===\n")
  cat("Final maternal model saved as 'final_maternal_model'\n")
  cat("Use summary(final_maternal_model) to see detailed results\n")
  
} else {
  cat("Multivariate maternal model failed to converge. Consider:\n")
  cat("1. Reducing the number of predictors\n")
  cat("2. Checking for perfect multicollinearity\n")
  cat("3. Using a different modeling approach\n")
}
  
  
summary(final_maternal_model)
  
  # ==============================================================================
  # 10. SUMMARY AND INTERPRETATION
  # ==============================================================================
  # Create Nature-style table for maternal care multivariate model
  library(gt)
  library(dplyr)
  
  # Create the data frame for maternal care behaviors
  maternal_care_results <- data.frame(
    Predictor = c("Intercept", 
                  "Nursing behaviors",
                  "Carrying behaviors", 
                  "Proximity behaviors",
                  "Active grooming",
                  "Protective interventions",
                  "Retrieving offspring",
                  "Comfort contact",
                  "Food sharing",
                  "Stationary touch",
                  "Muzzle contact"),
    Estimate = c(0.727, -0.001, -0.077, 0.384, 0.264, 0.066, -0.368, -0.419, 0.249, -0.581, 0.934),
    SE = c(0.280, 0.529, 0.534, 0.366, 0.749, 0.548, 0.553, 0.547, 0.407, 0.399, 0.370),
    z_value = c(2.592, -0.002, -0.144, 1.049, 0.352, 0.121, -0.665, -0.767, 0.611, -1.456, 2.525),
    P_value = c("<0.001", "0.999", "0.885", "0.294", "0.724", "0.904", "0.506", "0.443", "0.541", "0.145", "0.012")
  )
  
  # Create Nature-style gt table
  maternal_care_gt_table <- maternal_care_results %>%
    gt() %>%
    
    # Table title and subtitle
    tab_header(
      title = "Table 2. Multivariate model of maternal care behaviors predicting same-sex sexual behaviour",
      subtitle = "Effects of standardized maternal care behaviors on male offspring SSB"
    ) %>%
    
    # Column labels
    cols_label(
      Predictor = "Fixed Effects",
      Estimate = "Estimate", 
      SE = "SE",
      z_value = "z-value",
      P_value = "P-value"
    ) %>%
    
    # Format numbers
    fmt_number(
      columns = c(Estimate, SE, z_value),
      decimals = 3
    ) %>%
    
    # Bold significant rows (Intercept and Muzzle contact)
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(rows = c(1, 11))  # Intercept and Muzzle contact
    ) %>%
    
    # Highlight significant P-values
    tab_style(
      style = cell_text(weight = "bold", color = "black"),
      locations = cells_body(columns = P_value, rows = c(1, 11))
    ) %>%
    
    # Add footnote
    tab_footnote(
      footnote = "Generalized linear mixed model with Tweedie distribution and log link function. All behavioral predictors standardized to z-scores (mean = 0, SD = 1). Models weighted by proportion of sons sampled per mother. Random intercepts for MotherID included (variance = 8.25e-11). Dispersion parameter = 0.737. AIC = 54.0. n = 21 sons from 16 mothers."
    ) %>%
    
    # Style the table
    tab_options(
      table.font.size = 11,
      heading.title.font.size = 12,
      heading.subtitle.font.size = 10,
      column_labels.font.weight = "bold",
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      column_labels.border.bottom.style = "solid"
    ) %>%
    
    # Add source note
    tab_source_note("Significance: P < 0.001 (***), P < 0.01 (**), P < 0.05 (*)")
  
  # Display the table
  maternal_care_gt_table
  
  # Export using same method as social table
  library(gt)
  library(webshot2)
  gtsave(maternal_care_gt_table, "table2_maternal_care.png", expand = 10)
  # Export as HTML (this always works)
  gtsave(maternal_care_gt_table, "table2_maternal_care.html")


cat("\n=== ANALYSIS COMPLETE ===\n")






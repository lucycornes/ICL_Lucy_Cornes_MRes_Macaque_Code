# H2C: Maternal Social Connectedness and Sons' SSB - Tweedie Analysis
# Hypothesis: Males whose mothers were more socially connected—i.e., higher affiliative interaction rates—exhibit higher frequencies of SSB
# Using existing social_connectedness.csv dataset

# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(glmmTMB)  # For Tweedie mixed models
library(DHARMa)   # For residual diagnostics
library(ggeffects) # For prediction plots
library(performance) # For model comparison
library(tidyr)

# ==============================================================================
# 1. LOAD EXISTING SOCIAL CONNECTEDNESS DATASET AND MALE BEHAVIOR DATA
# ==============================================================================

cat("=== LOADING DATASETS ===\n")

# Load the existing social connectedness dataset
maternal_data <- read_csv("social_connectedness.csv")
male_behavior <- read_csv("male_behaviour_master_sheet2025.csv") 
mother_data <- read_csv("mother_master_sheet2025.csv")

cat("Maternal social connectedness data:", nrow(maternal_data), "mothers\n")
cat("Male behavior data:", nrow(male_behavior), "rows\n")
cat("Mother data:", nrow(mother_data), "rows\n")

# Check the structure of maternal data
cat("\nMaternal social connectedness dataset structure:\n")
print(names(maternal_data))
cat("\nSample of maternal data:\n")
print(head(maternal_data, 3))

# ==============================================================================
# 2. IDENTIFY AFFILIATIVE BEHAVIOR VARIABLES
# ==============================================================================

cat("\n=== IDENTIFYING AFFILIATIVE BEHAVIORS ===\n")

# Get all rate columns representing affiliative/social behaviors
affiliative_columns <- names(maternal_data)[grepl("^rate_", names(maternal_data))]
cat("Found", length(affiliative_columns), "affiliative behavior rate columns:\n")
cat(paste(affiliative_columns, collapse = ", "), "\n")

# Check which behaviors have variation and sufficient data
behavior_summary <- maternal_data %>%
  select(all_of(affiliative_columns)) %>%
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

cat("\nAffiliate behavioral summary:\n")
print(behavior_summary)

# Keep only behaviors with variation
good_behaviors <- behavior_summary$behavior[behavior_summary$has_variation]
cat("\nBehaviors with sufficient variation (", length(good_behaviors), "):\n")
cat(paste(good_behaviors, collapse = ", "), "\n")

if (length(good_behaviors) == 0) {
  stop("No affiliative behaviors have sufficient variation for analysis")
}

# ==============================================================================
# 3. LINK MATERNAL SOCIAL DATA TO SONS' SSB DATA
# ==============================================================================

cat("\n=== LINKING MATERNAL SOCIAL CONNECTEDNESS TO SONS' SSB DATA ===\n")

# Get mothers with social connectedness data
mothers_with_social_data <- unique(maternal_data$Actor)
cat("Mothers with social connectedness data:", length(mothers_with_social_data), "\n")

# Create mother-son links and add behavioral data
analysis_data <- mother_data %>%
  select(MotherID, MotherID_unique, SonID, BirthSeason, TotalSons) %>%
  filter(!is.na(SonID)) %>%
  filter(MotherID %in% mothers_with_social_data) %>%
  left_join(maternal_data, by = c("MotherID" = "Actor")) %>%
  left_join(male_behavior, by = c("SonID" = "ID")) %>%
  filter(!is.na(adjusted_SSB)) 


cat("Dataset before scan filtering:", nrow(analysis_data), "mother-son pairs\n")

# AFTER applying the threshold, recalculate proportion_sons:
threshold <- 5
analysis_data <- analysis_data %>% 
  filter(Scan_Count >= threshold) %>%
  group_by(MotherID) %>%
  mutate(
    sons_in_study = n(),                         # Recalculate on filtered data
    proportion_sons = sons_in_study / TotalSons  # Correct weights for models
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
# 4. STANDARDIZE AFFILIATIVE BEHAVIOR VARIABLES
# ==============================================================================

cat("\n=== STANDARDIZING AFFILIATIVE BEHAVIOR VARIABLES ===\n")

# Create standardized versions of all rate_ variables
rate_columns_in_data <- names(analysis_data)[grepl("^rate_", names(analysis_data))]

cat("Standardizing", length(rate_columns_in_data), "rate variables:\n")
cat(paste(rate_columns_in_data, collapse = ", "), "\n")

# Before standardization - show original scales
cat("\nOriginal scales (first 3 behaviors):\n")
for(i in 1:min(3, length(rate_columns_in_data))) {
  col <- rate_columns_in_data[i]
  if(col %in% names(analysis_data)) {
    values <- analysis_data[[col]]
    cat(col, ": Range =", round(min(values, na.rm = TRUE), 4), "to", 
        round(max(values, na.rm = TRUE), 4), 
        "| Mean =", round(mean(values, na.rm = TRUE), 4),
        "| SD =", round(sd(values, na.rm = TRUE), 4), "\n")
  }
}

# Standardize all rate_ columns
analysis_data <- analysis_data %>%
  mutate(across(all_of(rate_columns_in_data), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_std"))

# Get list of standardized column names
standardized_columns <- paste0(rate_columns_in_data, "_std")

cat("\nAfter standardization (first 3 behaviors):\n")
for(i in 1:min(3, length(standardized_columns))) {
  col <- standardized_columns[i]
  if(col %in% names(analysis_data)) {
    values <- analysis_data[[col]]
    original_col <- rate_columns_in_data[i]
    cat(col, ": Range =", round(min(values, na.rm = TRUE), 2), "to", 
        round(max(values, na.rm = TRUE), 2), 
        "| Mean =", round(mean(values, na.rm = TRUE), 2),
        "| SD =", round(sd(values, na.rm = TRUE), 2), "\n")
  }
}

cat("\n*** All rate variables now standardized (mean=0, SD=1) ***\n")
cat("*** Results will be interpreted as 'per 1 standard deviation increase' ***\n")

# ==============================================================================
# 5. PREPARE INDIVIDUAL AFFILIATIVE BEHAVIORS FOR ANALYSIS (STANDARDIZED)
# ==============================================================================

cat("\n=== PREPARING STANDARDIZED AFFILIATIVE BEHAVIORS ===\n")

# Update good_behaviors to use standardized versions
good_behaviors_std <- paste0(good_behaviors, "_std")

# Check which standardized behaviors still have variation in the linked dataset
linked_behavior_summary <- analysis_data %>%
  select(all_of(good_behaviors_std)) %>%
  summarise(across(everything(), list(
    mean = ~mean(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE),
    min = ~min(., na.rm = TRUE),
    max = ~max(., na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("behavior", "statistic"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(has_variation = sd > 0.01)  # Very small threshold since standardized

print(linked_behavior_summary)

# Final list of behaviors to analyze (standardized versions)
final_good_behaviors <- linked_behavior_summary$behavior[linked_behavior_summary$has_variation]
cat("\nStandardized behaviors with variation (", length(final_good_behaviors), "):\n")
for(i in 1:length(final_good_behaviors)) {
  behavior_name <- gsub("_std$", "", final_good_behaviors[i])
  behavior_name <- gsub("rate_", "", behavior_name)
  behavior_name <- gsub("_", " ", behavior_name)
  behavior_name <- tools::toTitleCase(behavior_name)
  cat(paste0(i, ". ", behavior_name, " (standardized)\n"))
}

# ==============================================================================
# 6. EXPLORATORY DATA ANALYSIS
# ==============================================================================

cat("\n=== SSB DISTRIBUTION ANALYSIS ===\n")
cat("SSB Score distribution:\n")
print(summary(analysis_data$adjusted_SSB))

# Check for zeros
zero_prop <- mean(analysis_data$adjusted_SSB == 0)
cat("Number of zeros:", sum(analysis_data$adjusted_SSB == 0), "out of", 
    nrow(analysis_data), "\n")
cat("Proportion of zeros:", round(zero_prop, 3), "\n")

# Check variance-to-mean ratio for overdispersion
var_mean_ratio <- var(analysis_data$adjusted_SSB) / mean(analysis_data$adjusted_SSB)
cat("Variance-to-mean ratio:", round(var_mean_ratio, 2), "\n")

# Visualize SSB distribution
hist_plot <- ggplot(analysis_data, aes(x = adjusted_SSB)) +
  geom_histogram(bins = 20, alpha = 0.7, fill = "lightblue", color = "black") +
  labs(title = "Distribution of Adjusted SSB Scores", 
       subtitle = paste0("Proportion of zeros: ", round(zero_prop, 3), 
                         " | Variance/Mean ratio: ", round(var_mean_ratio, 2)),
       x = "Adjusted SSB Score", y = "Count") +
  theme_minimal()
print(hist_plot)

cat("Tweedie family selected - handles zeros, continuous values, and overdispersion naturally\n")

# ==============================================================================
# 7. FUNCTION TO FIT TWEEDIE MODELS FOR EACH STANDARDIZED BEHAVIOR
# ==============================================================================

fit_affiliative_behavior_tweedie <- function(behavior_name, data) {
  cat("\n=== ANALYZING", toupper(gsub("_", " ", gsub("_std$", "", behavior_name))), "(STANDARDIZED) ===\n")
  
  # Check if behavior has sufficient variation
  behavior_values <- data[[behavior_name]]
  if (sd(behavior_values, na.rm = TRUE) < 0.01 || sum(!is.na(behavior_values)) < 10) {
    cat("Insufficient variation in", behavior_name, "- skipping\n")
    return(NULL)
  }
  
  cat("Range:", round(min(behavior_values, na.rm = TRUE), 2), "to", 
      round(max(behavior_values, na.rm = TRUE), 2), "SD units\n")
  cat("Non-missing observations:", sum(!is.na(behavior_values)), "\n")
  cat("Mean:", round(mean(behavior_values, na.rm = TRUE), 2), "SD units\n")
  cat("Standard deviation:", round(sd(behavior_values, na.rm = TRUE), 2), "SD units\n")
  
  # Fit Tweedie model
  tweedie_model <- NULL
  
  tryCatch({
    tweedie_model <- glmmTMB(
      formula = as.formula(paste("adjusted_SSB ~", behavior_name, "+ (1|MotherID)")),
      weights = proportion_sons,
      family = tweedie(),
      data = data
    )
    cat("✓ Tweedie model fitted successfully\n")
  }, error = function(e) {
    cat("✗ Tweedie model failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(tweedie_model)) {
    return(NULL)
  }
  
  # Model diagnostics
  cat("\n--- MODEL DIAGNOSTICS ---\n")
  tryCatch({
    suppressWarnings({
      sim_residuals <- simulateResiduals(tweedie_model)
      
      # Tests
      disp_test <- testDispersion(sim_residuals)
      uniform_test <- testUniformity(sim_residuals)
      zero_test <- testZeroInflation(sim_residuals)
    })
    
    cat("Dispersion test p-value:", round(disp_test$p.value, 4), "\n")
    cat("Uniformity test p-value:", round(uniform_test$p.value, 4), "\n")
    cat("Zero-inflation test p-value:", round(zero_test$p.value, 4), "\n")
    cat("Note: Diagnostics performed without weights (DHARMa limitation)\n")
    
    if(disp_test$p.value > 0.05 & uniform_test$p.value > 0.05 & zero_test$p.value > 0.05) {
      cat("✓ Model diagnostics look good!\n")
    } else {
      cat("⚠ Some diagnostic issues detected\n")
    }
    
  }, error = function(e) cat("Diagnostics failed:", e$message, "\n"))
  
  # Extract results
  coef_table <- summary(tweedie_model)$coefficients$cond
  
  # Get affiliative behavior effect
  if (behavior_name %in% rownames(coef_table)) {
    effect_estimate <- coef_table[behavior_name, "Estimate"]
    effect_se <- coef_table[behavior_name, "Std. Error"]
    effect_p <- coef_table[behavior_name, "Pr(>|z|)"]
    
    # For Tweedie with log-link, calculate percentage change
    percent_change <- (exp(effect_estimate) - 1) * 100
    
    cat("\nRESULTS:\n")
    cat("Effect estimate:", round(effect_estimate, 4), "\n")
    cat("Percent change in SSB per 1 SD increase:", round(percent_change, 1), "%\n")
    cat("P-value:", round(effect_p, 4), "\n")
    
    # 95% CI
    ci_lower <- effect_estimate - 1.96 * effect_se
    ci_upper <- effect_estimate + 1.96 * effect_se
    ci_lower_percent <- (exp(ci_lower) - 1) * 100
    ci_upper_percent <- (exp(ci_upper) - 1) * 100
    cat("95% CI for percent change: [", round(ci_lower_percent, 1), "%, ", 
        round(ci_upper_percent, 1), "%]\n")
    
    # Significance
    if (effect_p < 0.05) {
      if (effect_estimate > 0) {
        cat("*** SIGNIFICANT POSITIVE EFFECT ***\n")
      } else {
        cat("*** SIGNIFICANT NEGATIVE EFFECT ***\n")
      }
    } else {
      cat("Non-significant effect\n")
    }
    
    return(list(
      behavior = behavior_name,
      model = tweedie_model,
      effect_estimate = effect_estimate,
      effect_se = effect_se,
      effect_p = effect_p,
      percent_change = percent_change
    ))
  }
  
  return(NULL)
}

# ==============================================================================
# 8. ANALYZE ALL STANDARDIZED AFFILIATIVE BEHAVIORS INDIVIDUALLY
# ==============================================================================

cat("\n=== ANALYZING ALL STANDARDIZED AFFILIATIVE BEHAVIORS INDIVIDUALLY ===\n")

# Run analysis for each standardized behavior with sufficient variation
results_list <- list()
significant_behaviors <- c()

for (behavior in final_good_behaviors) {
  result <- fit_affiliative_behavior_tweedie(behavior, analysis_data)
  if (!is.null(result)) {
    results_list[[behavior]] <- result
    
    # Check if significant
    if (result$effect_p < 0.05) {
      significant_behaviors <- c(significant_behaviors, behavior)
    }
  }
}

# ==============================================================================
# 9. SUMMARY OF ALL RESULTS
# ==============================================================================

cat("\n=== SUMMARY OF ALL STANDARDIZED AFFILIATIVE BEHAVIOR EFFECTS ===\n")

if (length(results_list) > 0) {
  # Create summary table
  summary_results <- data.frame()
  
  for (behavior in names(results_list)) {
    result <- results_list[[behavior]]
    
    # Clean behavior name for display (remove _std suffix)
    clean_behavior_name <- gsub("_std$", "", behavior)
    clean_behavior_name <- gsub("rate_", "", clean_behavior_name)
    clean_behavior_name <- gsub("_", " ", clean_behavior_name)
    clean_behavior_name <- tools::toTitleCase(clean_behavior_name)
    
    summary_results <- rbind(summary_results, data.frame(
      Behavior = clean_behavior_name,
      Model = "Tweedie",
      Estimate = result$effect_estimate,
      SE = result$effect_se,
      Percent_Change = result$percent_change,
      P_value = result$effect_p,
      Significant = result$effect_p < 0.05,
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by significance and effect size
  summary_results <- summary_results[order(summary_results$P_value), ]
  
  cat("Results for all standardized affiliative behaviors:\n")
  cat("*** Effects shown as percent change per 1 standard deviation increase ***\n")
  cat("*** All models use Tweedie family (handles zeros and continuous values naturally) ***\n")
  print(summary_results)
  
  cat("\n=== SIGNIFICANT BEHAVIORS ===\n")
  significant_results <- summary_results[summary_results$Significant, ]
  if (nrow(significant_results) > 0) {
    for (i in 1:nrow(significant_results)) {
      row <- significant_results[i, ]
      direction <- ifelse(row$Percent_Change > 0, "POSITIVE", "NEGATIVE")
      cat(paste0(i, ". ", row$Behavior, ": ", round(row$Percent_Change, 1), 
                 "% change per 1 SD increase (p = ", round(row$P_value, 4), ") - ", direction, " EFFECT\n"))
    }
  } else {
    cat("No significant affiliative behavior effects found\n")
  }
  
} else {
  cat("No successful model fits\n")
}

# Primary significant effect
initiating_contact_model <- glmmTMB(
  adjusted_SSB ~ rate_initiating_contact_std + (1|MotherID),
  weights = proportion_sons,
  family = tweedie(),
  data = analysis_data
)


summary(initiating_contact_model)

# ==============================================================================
# 10. VISUALIZATION FOR SIGNIFICANT EFFECTS
# ==============================================================================

# Create prediction data for plotting
# For initiating contact model
initiating_range <- seq(min(analysis_data$rate_initiating_contact_std, na.rm = TRUE),
                        max(analysis_data$rate_initiating_contact_std, na.rm = TRUE),
                        length.out = 100)

initiating_pred_data <- data.frame(
  rate_initiating_contact_std = initiating_range,
  MotherID = NA  # Set to NA for population-level prediction
)

# Generate predictions for initiating contact model
# Add proportion_sons column to prediction data (can be any value since it's just for prediction structure)
initiating_pred_data$proportion_sons <- 1

# Get predictions on link scale (linear) and response scale
initiating_predictions_link <- predict(initiating_contact_model, 
                                       newdata = initiating_pred_data, 
                                       type = "link", 
                                       se.fit = TRUE,
                                       allow.new.levels = TRUE)

# Calculate confidence intervals on link scale then transform
initiating_lower_link <- initiating_predictions_link$fit - 1.96 * initiating_predictions_link$se.fit
initiating_upper_link <- initiating_predictions_link$fit + 1.96 * initiating_predictions_link$se.fit

# Transform to response scale
initiating_fit <- exp(initiating_predictions_link$fit)
initiating_lower <- exp(initiating_lower_link)
initiating_upper <- exp(initiating_upper_link)

initiating_plot_data <- data.frame(
  x = initiating_range,
  y = initiating_fit,
  lower = initiating_lower,
  upper = initiating_upper
)

# Nature theme
nature_theme <- theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", size = 10),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 10.5, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Plot: Initiating Contact Model
initiating_contact_plot <- ggplot() +
  # Add confidence interval ribbon
  geom_ribbon(data = initiating_plot_data,
              aes(x = x, ymin = lower, ymax = upper),
              fill = "grey90", alpha = 0.8) +
  # Add model line
  geom_line(data = initiating_plot_data,
            aes(x = x, y = y),
            color = "black", linewidth = 1.2) +
  # Add jittered data points - Using standardized version
  geom_jitter(data = analysis_data, 
              aes(x = rate_initiating_contact_std, y = adjusted_SSB),
              size = 1.5, color = "black", alpha = 0.7,
              width = 0.002, height = 0) +
  labs(
    x = "Standardised rate of initiating contact",
    y = "Adjusted SSB"
  ) +
  nature_theme

# Display plot
print(initiating_contact_plot)










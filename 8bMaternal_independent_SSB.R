# H2B: Maternal Care and Sons' SSB - Tweedie Model Analysis
# Hypothesis: Males whose mothers spent more time in maternal care behaviours exhibit higher frequencies of SSB
# Using existing maternal_behaviour_rate.csv dataset

# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(glmmTMB)  # For Tweedie mixed models
library(DHARMa)   # For residual diagnostics
library(ggeffects) # For prediction plots
library(performance) # For model comparison
library(tidyr)
library(tidyverse)


# ==============================================================================
# 1. LOAD EXISTING MATERNAL CARE DATASET AND MALE BEHAVIOR DATA
# ==============================================================================

cat("=== LOADING DATASETS ===\n")

# Load the existing maternal care dataset from the first part of H2B
maternal_data <- read_csv("maternal_behaviour_rate.csv")
male_behavior <- read_csv("male_behaviour_master_sheet2025.csv") 
mother_data <- read_csv("mother_master_sheet2025.csv")

cat("Maternal behavior rates:", nrow(maternal_data), "mothers\n")
cat("Male behavior data:", nrow(male_behavior), "rows\n")
cat("Mother data:", nrow(mother_data), "rows\n")

# Check the structure of maternal data
cat("\nMaternal care dataset structure:\n")
print(names(maternal_data))
cat("\nSample of maternal care data:\n")
print(head(maternal_data, 3))

# ==============================================================================
# 2. IDENTIFY MATERNAL CARE BEHAVIOR COLUMNS
# ==============================================================================

cat("\n=== IDENTIFYING MATERNAL CARE BEHAVIORS ===\n")

# Get all rate columns from the maternal data
rate_columns <- names(maternal_data)[grepl("^rate_", names(maternal_data))]
cat("Found", length(rate_columns), "maternal care rate columns:\n")
cat(paste(rate_columns, collapse = ", "), "\n")

# Check which behaviors have variation and sufficient data
behavior_summary <- maternal_data %>%
  select(all_of(rate_columns)) %>%
  summarise(across(everything(), list(
    mean = ~mean(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE),
    zeros = ~sum(. == 0, na.rm = TRUE),
    nonzeros = ~sum(. > 0, na.rm = TRUE),
    max = ~max(., na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("behavior", "statistic"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(has_variation = sd > 0 & nonzeros >= 3)  # Need at least 3 non-zero observations

cat("\nMaternal care behavior summary:\n")
print(behavior_summary)

# Keep only behaviors with variation
good_behaviors <- behavior_summary$behavior[behavior_summary$has_variation]
cat("\nBehaviors with sufficient variation (", length(good_behaviors), "):\n")
cat(paste(good_behaviors, collapse = ", "), "\n")

if (length(good_behaviors) == 0) {
  stop("No maternal care behaviors have sufficient variation for analysis")
}

# ==============================================================================
# 3. LINK MATERNAL CARE DATA TO SONS' SSB DATA
# ==============================================================================

cat("\n=== LINKING MATERNAL CARE TO SONS' SSB DATA ===\n")

# Get mothers with maternal care data
mothers_with_care <- unique(maternal_data$Actor)
cat("Mothers with maternal care data:", length(mothers_with_care), "\n")

# Create mother-son links and add behavioral data
analysis_data <- mother_data %>%
  select(MotherID, MotherID_unique, SonID, BirthSeason, TotalSons) %>%
  filter(!is.na(SonID)) %>%
  filter(MotherID %in% mothers_with_care) %>%
  left_join(maternal_data, by = c("MotherID" = "Actor")) %>%
  left_join(male_behavior, by = c("SonID" = "ID")) %>%
  filter(!is.na(adjusted_SSB)) 


cat("Dataset before scan filtering:", nrow(analysis_data), "mother-son pairs\n")

# Apply scan threshold
threshold <- 5
analysis_data <- analysis_data %>% filter(Scan_Count >= threshold) %>%
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
# 4. STANDARDIZE BEHAVIORAL VARIABLES
# ==============================================================================

cat("\n=== STANDARDIZING BEHAVIORAL VARIABLES ===\n")

# Get rate columns that exist in the analysis data
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
# 5. EXPLORATORY DATA ANALYSIS
# ==============================================================================

cat("\n=== SSB DISTRIBUTION ANALYSIS ===\n")
cat("SSB Score distribution:\n")
print(summary(analysis_data$rate_comfort_contact))

# Check for zeros
zero_prop <- mean(analysis_data$rate_comfort_contact == 0)
cat("Number of zeros:", sum(analysis_data$rate_comfort_contact == 0), "out of", 
    nrow(analysis_data), "\n")
cat("Proportion of zeros:", round(zero_prop, 3), "\n")

# Check variance-to-mean ratio for overdispersion
var_mean_ratio <- var(analysis_data$rate_comfort_contact) / mean(analysis_data$rate_comfort_contact)
cat("Variance-to-mean ratio:", round(var_mean_ratio, 2), "\n")

# Visualize SSB distribution
hist_plot <- ggplot(analysis_data, aes(x = rate_comfort_contact)) +
  geom_histogram(bins = 20, alpha = 0.7, fill = "lightblue", color = "black") +
  labs(title = "Distribution of Adjusted SSB Scores", 
       subtitle = paste0("Proportion of zeros: ", round(zero_prop, 3), 
                         " | Variance/Mean ratio: ", round(var_mean_ratio, 2)),
       x = "Adjusted SSB Score", y = "Count") +
  theme_minimal()
print(hist_plot)

cat("Tweedie family selected - handles zeros, continuous values, and overdispersion naturally\n")

# Check which standardized behaviors have variation in the linked dataset
good_behaviors_std <- paste0(good_behaviors, "_std")
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
# 6. FUNCTION TO FIT TWEEDIE MODELS FOR EACH STANDARDIZED BEHAVIOR
# ==============================================================================

fit_maternal_care_tweedie <- function(behavior_name, data) {
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
    sim_residuals <- simulateResiduals(tweedie_model)
    
    # Tests
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
    
  }, error = function(e) cat("Diagnostics failed:", e$message, "\n"))
  
  # Extract results
  coef_table <- summary(tweedie_model)$coefficients$cond
  
  # Get maternal care effect
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
# 7. ANALYZE ALL STANDARDIZED MATERNAL CARE BEHAVIORS
# ==============================================================================

cat("\n=== ANALYZING ALL STANDARDIZED MATERNAL CARE BEHAVIORS ===\n")

# Run analysis for each standardized behavior
results_list <- list()
significant_behaviors <- c()

for (behavior in final_good_behaviors) {
  result <- fit_maternal_care_tweedie(behavior, analysis_data)
  if (!is.null(result)) {
    results_list[[behavior]] <- result
    
    # Check if significant
    if (result$effect_p < 0.05) {
      significant_behaviors <- c(significant_behaviors, behavior)
    }
  }
}

# ==============================================================================
# 8. SUMMARY OF ALL RESULTS
# ==============================================================================

cat("\n=== SUMMARY OF ALL STANDARDIZED MATERNAL CARE EFFECTS ===\n")

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
  
  cat("Results for all standardized maternal care behaviors:\n")
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
    cat("No significant maternal care effects found\n")
  }
  
} else {
  cat("No successful model fits\n")
}

# ==============================================================================
# 9. SAVE RESULTS
# ==============================================================================

# Save analysis dataset
write_csv(analysis_data, "maternal_care_SSB_analysis_data_tweedie.csv")

# Save summary results
if (exists("summary_results") && nrow(summary_results) > 0) {
  write_csv(summary_results, "maternal_care_SSB_results_summary_tweedie.csv")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files saved:\n")
cat("- maternal_care_SSB_analysis_data_tweedie.csv (linked analysis dataset)\n")
cat("- maternal_care_SSB_results_summary_tweedie.csv (results summary)\n")

# ==============================================================================
# 10. FINAL HYPOTHESIS CONCLUSION
# ==============================================================================


cat("\n=== TWEEDIE MODEL ADVANTAGES ===\n")

cat("• Handles exact zeros naturally (no artificial adjustments) ✓\n")
cat("• Accommodates continuous positive values ✓\n")
cat("• Built-in overdispersion handling ✓\n")
cat("• Appropriate for right-skewed distributions ✓\n")
cat("• Standardized predictors give interpretable effect sizes ✓\n")
cat("• Random effects for clustering (MotherID) ✓\n")
cat("• Weighted by proportion of sons included ✓\n")

analysis_data <- read.csv("maternal_care_SSB_analysis_data_tweedie.csv")

# Significant Muzzle Contact Model
muzzle_contact_model <- glmmTMB(
  adjusted_SSB ~ rate_muzzle_contact_std + (1|MotherID),
  weights = proportion_sons,
  family = tweedie(),
  data = analysis_data
)
summary(muzzle_contact_model)

retrieve_offspring_model <- glmmTMB(
  adjusted_SSB ~ rate_retrieving_offspring_std + (1|MotherID),
  weights = proportion_sons,
  family = tweedie(),
  data = analysis_data
)
summary(retrieve_offspring_model)
# ==============================================================================
# 11. PLOTTING THESE MODELS 
# ==============================================================================

# Load required libraries
library(ggplot2)
library(glmmTMB)
library(dplyr)

# Assuming your models are already fitted:
# muzzle_contact_model and retrieve_offspring_model

# Create prediction data for plotting
# For muzzle contact model
muzzle_range <- seq(min(analysis_data$rate_muzzle_contact_std, na.rm = TRUE),
                    max(analysis_data$rate_muzzle_contact_std, na.rm = TRUE),
                    length.out = 100)

muzzle_pred_data <- data.frame(
  rate_muzzle_contact_std = muzzle_range,
  MotherID = NA  # Set to NA for population-level prediction
)

# Generate predictions for muzzle contact model
# Add proportion_sons column to prediction data (can be any value since it's just for prediction structure)
muzzle_pred_data$proportion_sons <- 1

# Get predictions on link scale (linear) and response scale
muzzle_predictions_link <- predict(muzzle_contact_model, 
                                   newdata = muzzle_pred_data, 
                                   type = "link", 
                                   se.fit = TRUE,
                                   allow.new.levels = TRUE)

# Calculate confidence intervals on link scale then transform
muzzle_lower_link <- muzzle_predictions_link$fit - 1.96 * muzzle_predictions_link$se.fit
muzzle_upper_link <- muzzle_predictions_link$fit + 1.96 * muzzle_predictions_link$se.fit

# Transform to response scale
muzzle_fit <- exp(muzzle_predictions_link$fit)
muzzle_lower <- exp(muzzle_lower_link)
muzzle_upper <- exp(muzzle_upper_link)

muzzle_plot_data <- data.frame(
  x = muzzle_range,
  y = muzzle_fit,
  lower = muzzle_lower,
  upper = muzzle_upper
)

# For retrieve offspring model - FIXED VERSION
retrieve_range <- seq(min(analysis_data$rate_retrieving_offspring_std, na.rm = TRUE),
                      max(analysis_data$rate_retrieving_offspring_std, na.rm = TRUE),
                      length.out = 100)

retrieve_pred_data <- data.frame(
  rate_retrieving_offspring_std = retrieve_range,  # FIXED: Added _std
  MotherID = NA  # Set to NA for population-level prediction
)

# Generate predictions for retrieve offspring model
# Add proportion_sons column to prediction data
retrieve_pred_data$proportion_sons <- 1

# Get predictions on link scale (linear) and response scale
retrieve_predictions_link <- predict(retrieve_offspring_model, 
                                     newdata = retrieve_pred_data, 
                                     type = "link", 
                                     se.fit = TRUE,
                                     allow.new.levels = TRUE)

# Calculate confidence intervals on link scale then transform
retrieve_lower_link <- retrieve_predictions_link$fit - 1.96 * retrieve_predictions_link$se.fit
retrieve_upper_link <- retrieve_predictions_link$fit + 1.96 * retrieve_predictions_link$se.fit

# Transform to response scale
retrieve_fit <- exp(retrieve_predictions_link$fit)
retrieve_lower <- exp(retrieve_lower_link)
retrieve_upper <- exp(retrieve_upper_link)

retrieve_plot_data <- data.frame(
  x = retrieve_range,
  y = retrieve_fit,
  lower = retrieve_lower,
  upper = retrieve_upper
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

# Plot 1: Muzzle Contact Model
plot1 <- ggplot() +
  # Add confidence interval ribbon
  geom_ribbon(data = muzzle_plot_data,
              aes(x = x, ymin = lower, ymax = upper),
              fill = "grey90", alpha = 0.8) +
  # Add model line
  geom_line(data = muzzle_plot_data,
            aes(x = x, y = y),
            color = "black", linewidth = 1.2) +
  # Add jittered data points - FIXED: Using standardized version
  geom_point(data = analysis_data, 
             aes(x = rate_muzzle_contact_std, y = adjusted_SSB),
             size = 1.5, color = "black", alpha = 0.7,
             position = position_jitter(width = 0.002, height = 0)) +
  labs(
    x = "Standardised muzzle contact rate",
    y = "Adjusted SSB"
  ) +
  nature_theme

# Plot 2: Retrieve Offspring Model
plot2 <- ggplot() +
  # Add confidence interval ribbon
  geom_ribbon(data = retrieve_plot_data,
              aes(x = x, ymin = lower, ymax = upper),
              fill = "grey90", alpha = 0.8) +
  # Add model line
  geom_line(data = retrieve_plot_data,
            aes(x = x, y = y),
            color = "black", linewidth = 1.2) +
  # Add jittered data points - FIXED: Using standardized version
  geom_jitter(data = analysis_data, 
              aes(x = rate_retrieving_offspring_std, y = adjusted_SSB),
              size = 1.5, color = "black", alpha = 0.7,
              width = 0.002, height = 1) +
  labs(
    x = "Standardised rate of retrieving offspring",
    y = "Adjusted SSB"
  ) +
  nature_theme

# Display plots
print(plot1)
print(plot2)

# Optional: Save plots as high-resolution figures
ggsave("muzzle_contact_plot.pdf", plot1, width = 6, height = 4.5, dpi = 300)
ggsave("retrieve_offspring_plot.pdf", plot2, width = 6, height = 4.5, dpi = 300)

# Create Nature-style paneled figure
library(patchwork)
combined_plot <- plot1 / plot2 + 
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = ')'
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 12, hjust = 0),
    plot.margin = margin(15, 15, 15, 15)
  )

print(combined_plot)

# Save combined plot
ggsave("figure1_combined.pdf", combined_plot, 
       width = 89, height = 120, units = "mm", dpi = 300)
ggsave("figurev2_combined.png", combined_plot, 
       width = 89, height = 120, units = "mm", dpi = 300)








# =============================================================================
# MATERNAL FACIAL REDNESS AND SONS' SSB ANALYSIS
# Research Question: Do females with redder faces have sons with higher SSB?
# =============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(corrplot)

cat("=== MATERNAL FACIAL REDNESS → SONS' SSB ANALYSIS ===\n")
cat("Research Question: Do females with redder faces have sons with higher SSB frequencies?\n")
cat("Scan threshold: >= 2\n")
cat("Primary predictor: redness_index (standardized)\n\n")

# =============================================================================
# STEP 1: LOAD AND EXAMINE DATASETS
# =============================================================================

cat("=== STEP 1: LOADING DATASETS ===\n")

# Load facial redness data
facial_redness <- read_csv("redness_data.csv")
cat("Facial redness data:", nrow(facial_redness), "females\n")

# Load male behavior data  
male_behavior <- read_csv("male_behaviour_master_sheet2025.csv")

cat("Male behavior data:", nrow(male_behavior), "males\n")

# Load mother-son linkage data
mother_data <- read_csv("mother_master_sheet2025.csv")
cat("Mother-son linkage data:", nrow(mother_data), "rows\n")

# Examine facial redness dataset structure
cat("\n=== FACIAL REDNESS DATASET STRUCTURE ===\n")
cat("Variables available:\n")
print(names(facial_redness))

cat("\nSample of facial redness data:\n")
print(head(facial_redness, 3))

# Check redness_index distribution
cat("\n=== REDNESS INDEX DISTRIBUTION ===\n")
cat("Redness index summary:\n")
print(summary(facial_redness$redness_index))

# Visualize redness distribution
redness_hist <- ggplot(facial_redness, aes(x = redness_index)) +
  geom_histogram(bins = 20, alpha = 0.7, fill = "red", color = "black") +
  labs(title = "Distribution of Maternal Facial Redness Index",
       x = "Redness Index", y = "Count") +
  theme_minimal()
print(redness_hist)

# Check for outliers or extreme values
cat("Redness index range:", round(min(facial_redness$redness_index), 3), "to", 
    round(max(facial_redness$redness_index), 3), "\n")
cat("Standard deviation:", round(sd(facial_redness$redness_index), 3), "\n")

# =============================================================================
# STEP 2: CREATE MOTHER-SON ANALYSIS DATASET
# =============================================================================

cat("\n=== STEP 2: CREATING ANALYSIS DATASET ===\n")

# Get mothers with facial redness data
mothers_with_redness <- unique(facial_redness$Actor)
cat("Mothers with facial redness data:", length(mothers_with_redness), "\n")

# Debug: Check data before joining
cat("Checking data before joining...\n")
cat("Mother data with sons:", sum(!is.na(mother_data$SonID)), "rows\n")
cat("Mothers in facial redness data:", length(mothers_with_redness), "\n")
cat("Male behavior rows:", nrow(male_behavior), "\n")

# Check if male behavior has the adjusted_SSB column
cat("Columns in male_behavior:", paste(names(male_behavior), collapse = ", "), "\n")
cat("Has adjusted_SSB?", "adjusted_SSB" %in% names(male_behavior), "\n")

# Check sample IDs to ensure matching will work
sample_mother_ids <- head(mother_data$MotherID[!is.na(mother_data$SonID)], 5)
sample_son_ids <- head(mother_data$SonID[!is.na(mother_data$SonID)], 5)
cat("Sample MotherIDs:", paste(sample_mother_ids, collapse = ", "), "\n")
cat("Sample SonIDs:", paste(sample_son_ids, collapse = ", "), "\n")
cat("Sample facial redness Actors:", paste(head(facial_redness$Actor, 5), collapse = ", "), "\n")
cat("Sample male behavior IDs:", paste(head(male_behavior$ID, 5), collapse = ", "), "\n")

# Create the analysis dataset by linking all three datasets
cat("\nStep 1: Starting with mother-son pairs...\n")
step1_data <- mother_data %>%
  filter(!is.na(SonID)) %>%
  filter(MotherID %in% mothers_with_redness)
cat("After filtering for mothers with redness data:", nrow(step1_data), "rows\n")

cat("\nStep 2: Adding facial redness data...\n")
step2_data <- step1_data %>%
  left_join(facial_redness, by = c("MotherID" = "Actor"))
cat("After adding redness data:", nrow(step2_data), "rows\n")
cat("Rows with redness_index:", sum(!is.na(step2_data$redness_index)), "\n")

cat("\nStep 3: Adding sons' behavioral data...\n")

# Debug the join operation
cat("Checking ID matching...\n")
unique_son_ids <- unique(step2_data$SonID[!is.na(step2_data$SonID)])
unique_behavior_ids <- unique(male_behavior$ID[!is.na(male_behavior$ID)])
cat("Unique SonIDs:", length(unique_son_ids), "\n")
cat("Unique behavior IDs:", length(unique_behavior_ids), "\n")

# Check for matches
matches <- sum(unique_son_ids %in% unique_behavior_ids)
cat("SonIDs that match behavior IDs:", matches, "out of", length(unique_son_ids), "\n")

# Show sample IDs for comparison
cat("Sample SonIDs:", paste(head(unique_son_ids, 5), collapse = ", "), "\n")
cat("Sample behavior IDs:", paste(head(unique_behavior_ids, 5), collapse = ", "), "\n")

# Check ID types (character vs numeric)
cat("SonID class:", class(step2_data$SonID), "\n")
cat("Behavior ID class:", class(male_behavior$ID), "\n")

# Convert both to character to ensure matching works
step2_data$SonID <- as.character(step2_data$SonID)
male_behavior$ID <- as.character(male_behavior$ID)

# Now perform the join
step3_data <- step2_data %>%
  left_join(male_behavior, by = c("SonID" = "ID"))

cat("After adding male behavior data:", nrow(step3_data), "rows\n")

# Check if we now have adjusted_SSB
if ("adjusted_SSB" %in% names(step3_data)) {
  cat("Rows with adjusted_SSB:", sum(!is.na(step3_data$adjusted_SSB)), "\n")
} else if ("adjusted_SSB.y" %in% names(step3_data)) {
  cat("✓ Found adjusted_SSB.y (from male behavior data)\n")
  cat("Rows with adjusted_SSB.y:", sum(!is.na(step3_data$adjusted_SSB.y)), "\n")
  # Rename to remove the .y suffix
  step3_data <- step3_data %>%
    rename(adjusted_SSB = adjusted_SSB.y,
           Scan_Count = Scan_Count.y,  # Use male scan count
           DSB_Score = DSB_Score.y,    # Use male DSB score
           SSB_Score = SSB_Score.y,    # Use male SSB score
           BirthSeason = BirthSeason.y) # Use male birth season
  cat("✓ Renamed columns to remove .y suffix\n")
} else {
  cat("ERROR: No adjusted_SSB column found\n")
  cat("Available columns after join:", paste(names(step3_data), collapse = ", "), "\n")
}

# Final filtering and processing
cat("\nStep 4: Final filtering and processing...\n")
analysis_dataset <- step3_data %>%
  # Remove rows without SSB data
  filter(!is.na(adjusted_SSB)) 

# Check final dataset
cat("Final dataset:", nrow(analysis_dataset), "rows\n")
cat("Final columns:", paste(names(analysis_dataset), collapse = ", "), "\n")

# Safely select only columns that exist
available_columns <- names(analysis_dataset)
desired_columns <- c("MotherID", "SonID", "BirthSeason", "TotalSons", "sons_in_study", "proportion_sons",
                     "redness_index", "red_mean", "green_mean", "blue_mean", "red_intensity", "rg_ratio",
                     "adjusted_SSB", "Scan_Count", "DSB_Score", "SSB_Score", "OffspringCount")

# Only select columns that actually exist
final_columns <- desired_columns[desired_columns %in% available_columns]
cat("Selecting columns:", paste(final_columns, collapse = ", "), "\n")

# Update analysis_dataset with safe column selection
analysis_dataset <- analysis_dataset %>%
  select(all_of(final_columns))

# Ensure we have the essential columns
essential_columns <- c("MotherID", "SonID", "redness_index", "adjusted_SSB", "Scan_Count")
missing_essential <- essential_columns[!essential_columns %in% names(analysis_dataset)]

if (length(missing_essential) > 0) {
  cat("ERROR: Missing essential columns:", paste(missing_essential, collapse = ", "), "\n")
  stop("Cannot proceed without essential columns")
} else {
  cat("✓ All essential columns present\n")
}

cat("Dataset before scan filtering:", nrow(analysis_dataset), "mother-son pairs\n")
cat("Mothers:", length(unique(analysis_dataset$MotherID)), "\n")
cat("Sons:", length(unique(analysis_dataset$SonID)), "\n")

# =============================================================================
# STEP 3: APPLY SCAN THRESHOLD AND FINAL FILTERING
# =============================================================================

cat("\n=== STEP 3: APPLYING SCAN THRESHOLD ===\n")

# Apply scan threshold >= 5
SCAN_THRESHOLD <- 5
# AFTER applying scan threshold, recalculate proportion_sons
filtered_dataset <- filtered_dataset %>%
  select(-sons_in_study, -proportion_sons) %>%  # Remove old calculations
  group_by(MotherID) %>%
  mutate(
    sons_in_study = n(),                         # Recalculate on filtered data
    proportion_sons = sons_in_study / first(TotalSons)  # Correct weights
  ) %>%
  ungroup()

cat("After scan threshold >=", SCAN_THRESHOLD, ":\n")
cat("Sample size:", nrow(filtered_dataset), "mother-son pairs\n")
cat("Mothers:", length(unique(filtered_dataset$MotherID)), "\n") 
cat("Sons:", length(unique(filtered_dataset$SonID)), "\n")

# Check if we have sufficient data for analysis
if (length(unique(filtered_dataset$MotherID)) < 5) {
  stop("Too few mothers for mixed effects model. Consider lowering scan threshold.")
}

if (nrow(filtered_dataset) < 20) {
  stop("Too few observations for reliable analysis. Consider lowering scan threshold.")
}

# =============================================================================
# STEP 4: STANDARDIZE REDNESS INDEX AND EXAMINE RELATIONSHIPS
# =============================================================================

cat("\n=== STEP 4: STANDARDIZING REDNESS INDEX ===\n")

# Standardize redness index (mean = 0, SD = 1)
filtered_dataset <- filtered_dataset %>%
  mutate(redness_index_std = as.numeric(scale(redness_index)))

# Check standardization
cat("Original redness index:\n")
cat("Range:", round(min(filtered_dataset$redness_index), 3), "to", 
    round(max(filtered_dataset$redness_index), 3), "\n")
cat("Mean:", round(mean(filtered_dataset$redness_index), 3), 
    "| SD:", round(sd(filtered_dataset$redness_index), 3), "\n")

cat("\nStandardized redness index:\n")
cat("Range:", round(min(filtered_dataset$redness_index_std), 2), "to", 
    round(max(filtered_dataset$redness_index_std), 2), "\n")
cat("Mean:", round(mean(filtered_dataset$redness_index_std), 2), 
    "| SD:", round(sd(filtered_dataset$redness_index_std), 2), "\n")

# =============================================================================
# STEP 5: EXAMINE SSB DISTRIBUTION AND CHOOSE MODEL FAMILY
# =============================================================================

cat("\n=== STEP 5: SSB DISTRIBUTION ANALYSIS ===\n")

cat("Sons' adjusted SSB distribution:\n")
print(summary(filtered_dataset$adjusted_SSB))

# Check for zeros
zero_count <- sum(filtered_dataset$adjusted_SSB == 0)
zero_prop <- mean(filtered_dataset$adjusted_SSB == 0)
cat("Zeros:", zero_count, "out of", nrow(filtered_dataset), 
    "(", round(zero_prop * 100, 1), "%)\n")

# Check variance-to-mean ratio for overdispersion
var_mean_ratio <- var(filtered_dataset$adjusted_SSB) / mean(filtered_dataset$adjusted_SSB)
cat("Variance-to-mean ratio:", round(var_mean_ratio, 2), "\n")

# Visualize SSB distribution
ssb_hist <- ggplot(filtered_dataset, aes(x = adjusted_SSB)) +
  geom_histogram(bins = 20, alpha = 0.7, fill = "lightblue", color = "black") +
  labs(title = "Distribution of Sons' Adjusted SSB Scores",
       subtitle = paste0("Zero proportion: ", round(zero_prop, 3), 
                         " | Variance/Mean ratio: ", round(var_mean_ratio, 2)),
       x = "Adjusted SSB Score", y = "Count") +
  theme_minimal()
print(ssb_hist)

cat("Model family choice: Tweedie (handles zeros, continuous values, overdispersion)\n")



# =============================================================================
# STEP 7: FIT PRIMARY MODEL - MATERNAL REDNESS INDEX → SONS' SSB
# =============================================================================

cat("\n=== STEP 7: FITTING PRIMARY MODEL ===\n")
cat("Model: adjusted_SSB ~ redness_index_std + (1|MotherID)\n")
cat("Family: Tweedie (log link)\n")
cat("Weights: proportion_sons\n\n")

# Fit the main model
redness_model <- NULL

tryCatch({
  redness_model <- glmmTMB(
    adjusted_SSB ~ redness_index_std + (1|MotherID),
    weights = proportion_sons,
    family = tweedie(),
    data = filtered_dataset
  )
  cat("✓ Tweedie model fitted successfully\n")
}, error = function(e) {
  cat("✗ Model fitting failed:", e$message, "\n")
  stop("Model fitting failed")
})

# Display model summary
cat("\n--- MODEL SUMMARY ---\n")
print(summary(redness_model))

# =============================================================================
# STEP 8: MODEL DIAGNOSTICS
# =============================================================================

cat("\n=== STEP 8: MODEL DIAGNOSTICS ===\n")

tryCatch({
  # Generate simulated residuals
  sim_residuals <- simulateResiduals(redness_model)
  
  # Diagnostic tests
  disp_test <- testDispersion(sim_residuals)
  uniform_test <- testUniformity(sim_residuals)
  zero_test <- testZeroInflation(sim_residuals)
  
  cat("Dispersion test p-value:", round(disp_test$p.value, 4), "\n")
  cat("Uniformity test p-value:", round(uniform_test$p.value, 4), "\n")
  cat("Zero-inflation test p-value:", round(zero_test$p.value, 4), "\n")
  
  # Overall diagnostic assessment
  if(disp_test$p.value > 0.05 & uniform_test$p.value > 0.05 & zero_test$p.value > 0.05) {
    cat("✓ Model diagnostics look good!\n")
  } else {
    cat("⚠ Some diagnostic concerns detected\n")
    if(disp_test$p.value <= 0.05) cat("  - Dispersion issues\n")
    if(uniform_test$p.value <= 0.05) cat("  - Non-uniform residuals\n") 
    if(zero_test$p.value <= 0.05) cat("  - Zero-inflation issues\n")
  }
  
  # Plot diagnostics
  plot(sim_residuals, main = "Model Diagnostic Plots")
  
}, error = function(e) {
  cat("Diagnostic tests failed:", e$message, "\n")
  cat("Model may still be valid - continuing with interpretation\n")
})

# =============================================================================
# STEP 9: EXTRACT AND INTERPRET RESULTS
# =============================================================================

cat("\n=== STEP 9: RESULTS INTERPRETATION ===\n")

# Extract coefficients
coef_table <- summary(redness_model)$coefficients$cond

if ("redness_index_std" %in% rownames(coef_table)) {
  
  # Get effect estimates
  effect_estimate <- coef_table["redness_index_std", "Estimate"]
  effect_se <- coef_table["redness_index_std", "Std. Error"]
  effect_z <- coef_table["redness_index_std", "z value"]
  effect_p <- coef_table["redness_index_std", "Pr(>|z|)"]
  
  # Calculate percentage change (log-link transformation)
  percent_change <- (exp(effect_estimate) - 1) * 100
  
  # Calculate 95% confidence interval
  ci_lower <- effect_estimate - 1.96 * effect_se
  ci_upper <- effect_estimate + 1.96 * effect_se
  ci_lower_percent <- (exp(ci_lower) - 1) * 100
  ci_upper_percent <- (exp(ci_upper) - 1) * 100
  
  # Display results
  cat("=== PRIMARY RESULTS ===\n")
  cat("Effect estimate (log scale):", round(effect_estimate, 4), "\n")
  cat("Standard error:", round(effect_se, 4), "\n")
  cat("Z-value:", round(effect_z, 3), "\n")
  cat("P-value:", round(effect_p, 4), "\n\n")
  
  cat("=== INTERPRETED EFFECT SIZE ===\n")
  cat("Percent change in SSB per 1 SD increase in maternal redness:", 
      round(percent_change, 1), "%\n")
  cat("95% CI for percent change: [", round(ci_lower_percent, 1), "%, ", 
      round(ci_upper_percent, 1), "%]\n\n")
  
  # Statistical significance
  if (effect_p < 0.001) {
    significance_level <- "highly significant (p < 0.001)"
  } else if (effect_p < 0.01) {
    significance_level <- "highly significant (p < 0.01)"
  } else if (effect_p < 0.05) {
    significance_level <- "significant (p < 0.05)"
  } else {
    significance_level <- "non-significant (p ≥ 0.05)"
  }
  
  cat("=== STATISTICAL SIGNIFICANCE ===\n")
  cat("Result:", significance_level, "\n")
  
  if (effect_p < 0.05) {
    direction <- ifelse(percent_change > 0, "POSITIVE", "NEGATIVE")
    cat("Direction:", direction, "association\n")
    
    if (direction == "POSITIVE") {
      cat("*** MOTHERS WITH REDDER FACES HAVE SONS WITH HIGHER SSB ***\n")
    } else {
      cat("*** MOTHERS WITH REDDER FACES HAVE SONS WITH LOWER SSB ***\n")
    }
  } else {
    cat("*** NO SIGNIFICANT ASSOCIATION BETWEEN MATERNAL FACIAL REDNESS AND SONS' SSB ***\n")
  }
  
} else {
  cat("Error: Could not extract redness_index_std coefficient\n")
}

# =============================================================================
# STEP 10: MODEL PREDICTIONS AND VISUALIZATION
# =============================================================================

# Clear R Studio Research Code for Redness Model Plot
# Model: adjusted_SSB ~ redness_index_std + (1|MotherID)

# Load required libraries
library(ggplot2)
library(glmmTMB)
library(dplyr)

# Assuming your model is already fitted as:
# redness_model <- glmmTMB(
#   adjusted_SSB ~ redness_index_std + (1|MotherID),
#   weights = proportion_sons,
#   family = tweedie(),
#   data = filtered_dataset
# )

# Generate predictions for straight trend line
pred_data <- data.frame(
  redness_index_std = seq(min(filtered_dataset$redness_index_std, na.rm = TRUE), 
                          max(filtered_dataset$redness_index_std, na.rm = TRUE), 
                          length.out = 100),
  proportion_sons = mean(filtered_dataset$proportion_sons, na.rm = TRUE)
)

# Get predictions with confidence intervals (produces straight line)
pred_results <- predict(redness_model, newdata = pred_data, 
                        type = "response", se.fit = TRUE, re.form = NA)

# Create prediction dataframe for plotting
pred_data <- data.frame(
  x = pred_data$redness_index_std,
  predicted = pred_results$fit,
  conf.low = pred_results$fit - 1.96 * pred_results$se.fit,
  conf.high = pred_results$fit + 1.96 * pred_results$se.fit
)

# Extract model statistics for subtitle
model_summary <- summary(redness_model)
coef_estimate <- round(model_summary$coefficients$cond[2, "Estimate"], 3)
coef_se <- round(model_summary$coefficients$cond[2, "Std. Error"], 3)
coef_p <- round(model_summary$coefficients$cond[2, "Pr(>|z|)"], 3)

# Create the research plot
redness_plot <- ggplot() +
  # Add confidence interval ribbon
  geom_ribbon(data = pred_data, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = "#2c3e50") +
  
  # Add raw data points
  geom_point(data = filtered_dataset, 
             aes(x = redness_index_std, y = adjusted_SSB),
             color = "#34495e", size = 2.5, alpha = 0.7) +
  
  # Add straight trend line
  geom_line(data = pred_data, 
            aes(x = x, y = predicted), 
            color = "#2c3e50", linewidth = 1.2) +
  
  # Research-appropriate labels
  labs(
    title = "Maternal Facial Redness and Son Sexual Behavior",
    x = "Facial Redness Index (standardised)",
    y = "Adjusted Same-Sex Behavior (SSB)",
  ) +
  
  # Clean research theme
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "#2c3e50"),
    plot.subtitle = element_text(size = 11, color = "#7f8c8d", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "#95a5a6", hjust = 0),
    axis.title = element_text(size = 16, color = "#2c3e50"),
    axis.text = element_text(size = 10, color = "#34495e"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # Set appropriate axis breaks
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))

# Display the plot
print(redness_plot)

# Optional: Save the plot for publication
# ggsave("facial_redness_ssb_plot.png", redness_plot, 
#        width = 8, height = 6, dpi = 300, bg = "white")

# Print model summary for reference
cat("=== FACIAL REDNESS MODEL SUMMARY ===\n")
summary(redness_model)
# =============================================================================
# STEP 11: ADDITIONAL REDNESS MEASURES EXPLORATION
# =============================================================================

cat("\n=== STEP 11: EXPLORING OTHER REDNESS MEASURES ===\n")

# Test other redness-related variables
other_redness_vars <- c("red_intensity", "rg_ratio", "normalized_red_excess", 
                        "saturation_weighted_red")

# Check which variables are available in the dataset
available_redness_vars <- other_redness_vars[other_redness_vars %in% names(filtered_dataset)]
cat("Available additional redness measures:", length(available_redness_vars), "\n")
cat("Variables:", paste(available_redness_vars, collapse = ", "), "\n\n")

# Test each additional measure
additional_results <- list()

for (var in available_redness_vars) {
  cat("Testing", var, "...\n")
  
  # Standardize the variable
  var_std_name <- paste0(var, "_std")
  filtered_dataset[[var_std_name]] <- as.numeric(scale(filtered_dataset[[var]]))
  
  # Fit model
  tryCatch({
    formula_str <- paste("adjusted_SSB ~", var_std_name, "+ (1|MotherID)")
    temp_model <- glmmTMB(
      formula = as.formula(formula_str),
      weights = proportion_sons,
      family = tweedie(),
      data = filtered_dataset
    )
    
    # Extract results
    temp_coef <- summary(temp_model)$coefficients$cond
    if (var_std_name %in% rownames(temp_coef)) {
      temp_estimate <- temp_coef[var_std_name, "Estimate"]
      temp_p <- temp_coef[var_std_name, "Pr(>|z|)"]
      temp_percent <- (exp(temp_estimate) - 1) * 100
      
      additional_results[[var]] <- list(
        estimate = temp_estimate,
        p_value = temp_p,
        percent_change = temp_percent
      )
      
      cat("  ", var, ": ", round(temp_percent, 1), "% change, p = ", 
          round(temp_p, 4), "\n")
    }
    
  }, error = function(e) {
    cat("  ", var, ": Model failed\n")
  })
}

# =============================================================================
# STEP 12: FINAL SUMMARY AND CONCLUSIONS
# =============================================================================

cat("\n=== STEP 12: FINAL SUMMARY ===\n")
cat("Research Question: Do females with redder faces have sons with higher SSB frequencies?\n\n")

cat("DATASET SUMMARY:\n")
cat("• Final sample:", nrow(filtered_dataset), "mother-son pairs\n")
cat("• Mothers:", length(unique(filtered_dataset$MotherID)), "\n")
cat("• Sons:", length(unique(filtered_dataset$SonID)), "\n")
cat("• Scan threshold: ≥", SCAN_THRESHOLD, "\n\n")

cat("PRIMARY ANALYSIS RESULTS:\n")
cat("• Predictor: Maternal redness index (standardized)\n")
cat("• Outcome: Sons' adjusted SSB rate\n")
cat("• Model: Tweedie mixed-effects (random intercept for mother)\n")
cat("• Effect size:", round(percent_change, 1), "% change per 1 SD increase\n")
cat("• P-value:", round(effect_p, 4), "\n")
cat("• 95% CI: [", round(ci_lower_percent, 1), "%, ", round(ci_upper_percent, 1), "%]\n\n")

if (effect_p < 0.05) {
  cat("*** HYPOTHESIS SUPPORTED ***\n")
  if (percent_change > 0) {
    cat("Mothers with redder faces DO have sons with significantly higher SSB frequencies\n")
  } else {
    cat("Mothers with redder faces have sons with significantly LOWER SSB frequencies\n")
  }
  cat("This represents a statistically significant association.\n")
} else {
  cat("*** HYPOTHESIS NOT SUPPORTED ***\n")
  cat("No significant association found between maternal facial redness and sons' SSB frequencies\n")
  cat("The relationship may be too weak to detect with this sample size or may not exist.\n")
}

cat("\nBIOLOGICAL INTERPRETATION:\n")
cat("Facial redness can indicate hormonal status, health, or dominance signals.\n")
cat("If significant, this could suggest maternal phenotypic traits influence sons' sexual behaviors.\n")

# Display model summary one final time
cat("\n=== FINAL MODEL SUMMARY ===\n")
print(summary(redness_model))

cat("\n*** ANALYSIS COMPLETE ***\n")
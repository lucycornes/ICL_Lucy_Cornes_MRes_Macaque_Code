# Social Connectedness Hypothesis Analysis - Updated Version
# ============================================================
# Hypothesis: Females with higher social connectedness have higher offspring count

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(lubridate)
library(ggplot2)
library(glmmTMB)

# 1. LOAD AND PREPARE DATA
# ========================
All_Focals_Cleaned <- read_csv("New_All_Focals_Cleaned.csv")
Offspring_Count <- read_csv("Offspring Count.csv")

# Clean data
All_Focals_Cleaned <- All_Focals_Cleaned %>%
  mutate(
    Actor = str_trim(Actor),
    Receiver = str_trim(Receiver),
    Behavior = str_trim(Behavior),
    Duration_s = as.numeric(Duration_s)
  )

# Prepare pedigree data
pedigree_data <- Offspring_Count %>%
  mutate(
    AnimalID = as.character(AnimalID),
    AnimalID = gsub("EE", "E", AnimalID),
    Age = year(Sys.Date()) - BirthSeason
  )

# 2. CALCULATE OBSERVATION TIME PER ACTOR
# =======================================
obs_time <- All_Focals_Cleaned %>%
  distinct(Actor, Log_File_Name) %>%
  group_by(Actor) %>%
  summarise(
    NumFocals = n(),
    TotalObsTime_s = NumFocals * 1800,
    .groups = "drop"
  )

# Add total behavior time per actor
actor_total_time <- All_Focals_Cleaned %>%
  group_by(Actor) %>%
  summarise(Total_seconds = sum(Duration_s, na.rm = TRUE), .groups = "drop") %>%
  mutate(Total_hours = Total_seconds / 3600)

# 3. DEFINE UPDATED SOCIAL CONNECTEDNESS BEHAVIOR CATEGORIES
# ==========================================================
# Updated social connectedness categories based on your specification
# Note: These behaviors are NOT filtered to specific receivers since we want general social connectivity

social_categories <- list(
  "Comfort_Contact" = c("Hugging", "Huddling", "Clasped sleeping"),
  "Proximity" = c("Proximity to offspring", "P+G (give)", "P+G (receive)"), 
  "Grooming" = c("N+G(give)", "N+G(receive)", "Groom (give)", "Groom (receive)", "Soliciting grooming"),
  "Food_Sharing" = c("Food sharing (give)", "Food sharing (receive)"),
  "Muzzle_Contact" = c("Muzzle contact"),
  "Initiating_Contact" = c("Initiating contact"),
  "Passive_Body_Contact" = c("Passive body contact"),
  "Approaching" = c("Approaching"),
  "Stationary_Proximity" = c("Stationary in proximity"),
  "Teeth_Chatter" = c("Teeth chatter"),
  "Reconciliation" = c("Reconciliation"),
  "Lip_Smacking" = c("Lip smacking"),
  "Food_Sharing_Give" = c("Food sharing (give)"),
  "Food_Sharing_Receive" = c("Food sharing (receive)")
)

# 4. CALCULATE SOCIAL CONNECTEDNESS RATES
# =======================================

# Function to calculate rates for social behaviors (no receiver filtering)
calculate_social_rate <- function(behaviors, behavior_name) {
  All_Focals_Cleaned %>%
    filter(Behavior %in% behaviors) %>%
    group_by(Actor) %>%
    summarise(!!paste0("Total_", behavior_name, "_s") := sum(Duration_s, na.rm = TRUE), 
              .groups = "drop")
}

# Calculate rates for each category using dynamic processing
social_data_list <- list()

for(category_name in names(social_categories)) {
  behaviors <- social_categories[[category_name]]
  clean_name <- gsub(" ", "_", tolower(category_name))
  
  # Calculate rates for this category
  social_data_list[[clean_name]] <- calculate_social_rate(behaviors, clean_name)
}

# 5. COMBINE AND CALCULATE RATES PER HOUR
# =======================================
social_rates <- actor_total_time


# Join all behavioral data
for(category_name in names(social_data_list)) {
  social_rates <- social_rates %>%
    left_join(social_data_list[[category_name]], by = "Actor")
}

# Ensure all Total_ columns exist (set to 0 if missing)
all_category_names <- gsub(" ", "_", tolower(names(social_categories)))
for(cat_name in all_category_names) {
  total_col_name <- paste0("Total_", cat_name, "_s")
  if(!total_col_name %in% names(social_rates)) {
    social_rates[[total_col_name]] <- 0
  }
}

# Replace missing values with 0 and calculate rates per hour
social_rates <- social_rates %>%
  mutate(across(starts_with("Total_"), ~replace_na(., 0))) %>%
  mutate(
    rate_comfort_contact = (Total_comfort_contact_s / Total_seconds) * 60,
    rate_proximity = (Total_proximity_s / Total_seconds) * 60,
    rate_grooming = (Total_grooming_s / Total_seconds) * 60,
    rate_food_sharing = (Total_food_sharing_s / Total_seconds) * 60,
    rate_muzzle_contact = (Total_muzzle_contact_s / Total_seconds) * 60,
    rate_initiating_contact = (Total_initiating_contact_s / Total_seconds) * 60,
    rate_passive_body_contact = (Total_passive_body_contact_s / Total_seconds) * 60,
    rate_approaching = (Total_approaching_s / Total_seconds) * 60,
    rate_stationary_proximity = (Total_stationary_proximity_s / Total_seconds) * 60,
    rate_teeth_chatter = (Total_teeth_chatter_s / Total_seconds) * 60,
    rate_reconciliation = (Total_reconciliation_s / Total_seconds) * 60,
    rate_lip_smacking = (Total_lip_smacking_s / Total_seconds) * 60,
    rate_food_sharing_give = (Total_food_sharing_give_s / Total_seconds) * 60,
    rate_food_sharing_receive = (Total_food_sharing_receive_s / Total_seconds) * 60
  )

# Select only the columns we need
rate_columns <- names(social_rates)[grepl("^rate_", names(social_rates))]
social_rates <- social_rates[, c("Actor", "Total_seconds", "Total_hours", rate_columns)]

# 6. CREATE ANALYSIS DATASET
# ==========================
analysis_data2 <- social_rates %>%
  left_join(pedigree_data[, c("AnimalID", "offspring_count", "Age")], 
            by = c("Actor" = "AnimalID")) %>%
  # Remove rows with missing key data
  filter(!is.na(offspring_count) & !is.na(Age)) %>%
  # Replace any remaining NAs with 0 and add log_age
  mutate(
    across(starts_with("rate_"), ~replace_na(., 0)),
    log_age = log(Age)
  )

write.csv(analysis_data2, "social_connectedness.csv")

# 7. EXPLORATORY DATA ANALYSIS
# ============================
cat("Analysis dataset contains", nrow(analysis_data), "females\n")

# Check distribution of offspring count
hist(analysis_data$offspring_count, 
     main = "Distribution of Offspring Count (All Females)", 
     xlab = "Offspring Count")

# Density plot
ggplot(analysis_data, aes(x = offspring_count)) +
  geom_density(fill = "lightgreen", alpha = 0.7) +
  labs(title = "Distribution of Offspring Count", 
       x = "Offspring Count", y = "Density") +
  theme_minimal()

# Test for normality
shapiro_test <- shapiro.test(analysis_data$offspring_count)
cat("Shapiro-Wilk p-value:", round(shapiro_test$p.value, 4), "\n")

# Check for overdispersion with Poisson
poisson_test <- glmmTMB(offspring_count ~ rate_comfort_contact + offset(log_age), 
                        family = poisson, 
                        data = analysis_data2)

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

# 8. CHOOSE APPROPRIATE MODEL FAMILY
# ==================================
# Based on dispersion test results, choose appropriate family
if (dispersion_result < 1.5) {
  model_family <- poisson()
  cat("Using Poisson family\n")
} else if (dispersion_result < 2) {
  model_family <- poisson()  # or negative binomial
  cat("Using Poisson family (could consider negative binomial)\n")
} else {
  model_family <- gaussian()
  cat("Using Gaussian family due to overdispersion\n")
}


# 9. SOCIAL CONNECTEDNESS MODEL ANALYSIS
# ======================================

# Full model with all 14 social connectedness categories
full_social_model <- glmmTMB(
  offspring_count ~ rate_comfort_contact + rate_proximity + rate_grooming + 
    rate_food_sharing + rate_muzzle_contact + rate_initiating_contact + 
    rate_passive_body_contact + rate_approaching + rate_stationary_proximity + 
    rate_teeth_chatter + rate_reconciliation + rate_lip_smacking + 
    rate_food_sharing_give + rate_food_sharing_receive + offset(log_age),
  family = poisson(),
  data = analysis_data2
)

summary(full_social_model)

# 10. STEPWISE MODEL SELECTION
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

# Perform stepwise elimination with all 14 social connectedness categories
social_vars <- c("rate_comfort_contact", "rate_proximity", "rate_grooming", 
                 "rate_food_sharing", "rate_muzzle_contact", "rate_initiating_contact", 
                 "rate_passive_body_contact", "rate_approaching", "rate_stationary_proximity", 
                 "rate_teeth_chatter", "rate_reconciliation", "rate_lip_smacking", 
                 "rate_food_sharing_give", "rate_food_sharing_receive")

stepwise_result <- stepwise_elimination(analysis_data, social_vars, model_family)
final_social_model <- stepwise_result$final_model

cat("\n=== FINAL SOCIAL CONNECTEDNESS MODEL ===\n")
summary(final_social_model)

social_connect_model <- glmmTMB(
  offspring_count ~ rate_food_sharing_receive + offset(log_age),
  family = poisson(),
  data = analysis_data
)

# Display model summary
summary(social_connect_model)

# 11. MODEL VALIDATION
# ====================

# Check model with DHARMa
if(require(DHARMa, quietly = TRUE)) {
  residuals_dharma <- simulateResiduals(final_social_model)
  plot(residuals_dharma, main = "DHARMa Residual Diagnostics")
} else {
  cat("Install DHARMa package for residual diagnostics: install.packages('DHARMa')\n")
}

# 12. MODEL VISUALIZATION
# ======================

# Get significant predictors
significant_vars <- stepwise_result$final_vars

# Create plots for significant predictors
for(var in significant_vars) {
  # Create prediction data
  pred_data <- data.frame(
    x = seq(min(analysis_data[[var]]), max(analysis_data[[var]]), length.out = 100)
  )
  names(pred_data)[1] <- var
  
  # Add other variables at their means
  for(other_var in setdiff(significant_vars, var)) {
    pred_data[[other_var]] <- mean(analysis_data[[other_var]])
  }
  pred_data$log_age <- mean(analysis_data$log_age)
  
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
    labs(x = paste("Rate of", behavior_name, "(min/hour)"),
         y = "Offspring Count",
         title = paste("Effect of", behavior_name, "on Offspring Count")) +
    theme_minimal()
  
  print(p)
}

# 13. SAVE RESULTS
# ================
write.csv(analysis_data, "social_connectedness.csv", row.names = FALSE)

# Create summary table of results
if(length(significant_vars) > 0) {
  model_summary <- summary(final_social_model)
  coef_table <- model_summary$coefficients$cond
  results_table <- data.frame(
    Variable = rownames(coef_table),
    Estimate = coef_table[,1],
    SE = coef_table[,2],
    Z_value = coef_table[,3],
    P_value = coef_table[,4]
  )
  write.csv(results_table, "social_connectedness_model_results_updated.csv", row.names = FALSE)
  
  # Print results summary
  cat("\n=== MODEL RESULTS SUMMARY ===\n")
  print(results_table)
}

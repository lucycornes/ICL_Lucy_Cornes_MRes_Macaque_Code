# Maternal Care Hypothesis Analysis - Updated Version
# ===================================================
# Hypothesis: Mothers who spend more time engaging in maternal care behaviours—
# including time spent nursing, grooming their offspring, retrieving, and 
# protective interventions—have a greater number of offspring.

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
# Read as all characters first, then convert
All_Focals_Cleaned <- read_csv("New_All_Focals_Cleaned.csv",
                               col_types = "ccccccccccc",  # All as character
                               show_col_types = FALSE)

# Then convert the numeric columns manually
All_Focals_Cleaned <- All_Focals_Cleaned %>%
  mutate(
    Time_Relative_s = as.numeric(Time_Relative_s),
    Duration_s = as.numeric(Duration_s)
  )
Offspring_Count <- read_csv("Offspring Count.csv")
Infant_count <- read_csv("Infant_count.csv")

# Check what values are causing the problem
problematic_time <- All_Focals_Cleaned$Time_Relative_s[is.na(as.numeric(All_Focals_Cleaned$Time_Relative_s))]
print("Problematic Time_Relative_s values:")
print(unique(problematic_time))

# Check how many NAs are introduced
cat("Original non-NA values:", sum(!is.na(All_Focals_Cleaned$Time_Relative_s)), "\n")
test_conversion <- as.numeric(All_Focals_Cleaned$Time_Relative_s)
cat("After conversion non-NA values:", sum(!is.na(test_conversion)), "\n")
cat("NAs introduced:", sum(is.na(test_conversion)) - sum(is.na(All_Focals_Cleaned$Time_Relative_s)), "\n")



# Check what parsing issues occurred
problems(New_All_Focals_Cleaned)

# Clean data
All_Focals_Cleaned <- New_All_Focals_Cleaned %>%
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

# 3. DEFINE UPDATED MATERNAL CARE BEHAVIOR CATEGORIES
# ===================================================

# Updated maternal care categories based on your specification
maternal_categories <- list(
  "Nursing Behaviors" = c("Nursing", "N+G(give)", "N+G(receive)", "N+LA (G)", "N+LA(R)", "N+HA(g)"),
  "Carrying Behaviors" = c("Carrying offspring", "C+HA(G)", "C+HA(R)", "C+LA(G)", "C+LA(R)"),
  "Proximity Behaviors" = c("Proximity to offspring", "P+G (give)", "P+G (receive)"),
  "Active Grooming" = c("N+G(give)", "Groom (give)"),
  "Protective Interventions" = c("Protective interventions"),
  "Retrieving Offspring" = c("Retrieving offspring"),
  "Soliciting Grooming" = c("Soliciting grooming"),
  "Comfort Contact" = c("Hugging", "Huddling", "Clasped sleeping"),
  "Food Sharing" = c("Food sharing (give)"),
  "Stationary Touch" = c("Stationary touch"),
  "Muzzle Contact" = c("Muzzle contact"),
  "Initiating Contact" = c("Initiating contact"),
  "Reconciliation" = c("Reconciliation"),
  "Teeth Chatter" = c("Teeth chatter")
)

# Categories that are already offspring-specific (don't need filtering)
offspring_specific_categories <- c("Nursing Behaviors", "Carrying Behaviors", "Proximity Behaviors", 
                                   "Protective Interventions", "Retrieving Offspring")

# Categories that need filtering to "Own Offspring"
need_filtering_categories <- c("Active Grooming", "Soliciting Grooming", "Comfort Contact", 
                               "Food Sharing", "Stationary Touch", "Muzzle Contact", 
                               "Initiating Contact", "Reconciliation", "Teeth Chatter")

# 4. CALCULATE MATERNAL CARE RATES
# ================================

# Function for behaviors that DON'T need filtering (already offspring-specific)
calculate_maternal_rate_no_filter <- function(behaviors, behavior_name) {
  All_Focals_Cleaned %>%
    filter(Behavior %in% behaviors) %>%
    group_by(Actor) %>%
    summarise(!!paste0("Total_", behavior_name, "_s") := sum(Duration_s, na.rm = TRUE), 
              .groups = "drop")
}

# Function for behaviors that DO need filtering to "Own Offspring"
calculate_maternal_rate_filtered <- function(behaviors, behavior_name) {
  All_Focals_Cleaned %>%
    filter(Behavior %in% behaviors & Receiver == "Own Offspring") %>%
    group_by(Actor) %>%
    summarise(!!paste0("Total_", behavior_name, "_s") := sum(Duration_s, na.rm = TRUE), 
              .groups = "drop")
}

# Calculate rates for each category
maternal_data_list <- list()

for(category_name in names(maternal_categories)) {
  behaviors <- maternal_categories[[category_name]]
  clean_name <- gsub(" ", "_", tolower(category_name))
  
  if(category_name %in% offspring_specific_categories) {
    # Don't filter - already offspring-specific
    maternal_data_list[[clean_name]] <- calculate_maternal_rate_no_filter(behaviors, clean_name)
  } else {
    # Filter to "Own Offspring"
    maternal_data_list[[clean_name]] <- calculate_maternal_rate_filtered(behaviors, clean_name)
  }
}

# 5. COMBINE AND CALCULATE RATES PER HOUR
# =======================================
maternal_rates <- actor_total_time

# Join all behavioral data
for(category_name in names(maternal_data_list)) {
  maternal_rates <- maternal_rates %>%
    left_join(maternal_data_list[[category_name]], by = "Actor")
}

# Ensure all Total_ columns exist (set to 0 if missing)
all_category_names <- gsub(" ", "_", tolower(names(maternal_categories)))
for(cat_name in all_category_names) {
  total_col_name <- paste0("Total_", cat_name, "_s")
  if(!total_col_name %in% names(maternal_rates)) {
    maternal_rates[[total_col_name]] <- 0
  }
}

# Replace missing values with 0 and calculate rates per hour
maternal_rates <- maternal_rates %>%
  mutate(across(starts_with("Total_"), ~replace_na(., 0))) %>%
  mutate(
    rate_nursing_behaviors = (Total_nursing_behaviors_s / Total_seconds) * 60,
    rate_carrying_behaviors = (Total_carrying_behaviors_s / Total_seconds) * 60,
    rate_proximity_behaviors = (Total_proximity_behaviors_s / Total_seconds) * 60,
    rate_active_grooming = (Total_active_grooming_s / Total_seconds) * 60,
    rate_protective_interventions = (Total_protective_interventions_s / Total_seconds) * 60,
    rate_retrieving_offspring = (Total_retrieving_offspring_s / Total_seconds) * 60,
    rate_soliciting_grooming = (Total_soliciting_grooming_s / Total_seconds) * 60,
    rate_comfort_contact = (Total_comfort_contact_s / Total_seconds) * 60,
    rate_food_sharing = (Total_food_sharing_s / Total_seconds) * 60,
    rate_stationary_touch = (Total_stationary_touch_s / Total_seconds) * 60,
    rate_muzzle_contact = (Total_muzzle_contact_s / Total_seconds) * 60,
    rate_initiating_contact = (Total_initiating_contact_s / Total_seconds) * 60,
    rate_reconciliation = (Total_reconciliation_s / Total_seconds) * 60,
    rate_teeth_chatter = (Total_teeth_chatter_s / Total_seconds) * 60
  )

# Select only the columns we need
rate_columns <- names(maternal_rates)[grepl("^rate_", names(maternal_rates))]
maternal_rates <- maternal_rates[, c("Actor", "Total_seconds", "Total_hours", rate_columns)]

# 6. CREATE ANALYSIS DATASET
# ==========================
analysis_data <- maternal_rates %>%
  left_join(pedigree_data[, c("AnimalID", "offspring_count", "Age")], 
            by = c("Actor" = "AnimalID")) %>%
  left_join(Infant_count[, c("AnimalID", "Infant Age")],
            by = c("Actor" = "AnimalID")) %>%
  # Filter to mothers with infants under 1 year
  filter(`Infant Age` < 1) %>%
  # Remove duplicates (mothers with multiple infants under 1)
  distinct(Actor, .keep_all = TRUE) %>%
  # Replace any remaining NAs with 0 and add log_age
  mutate(
    across(starts_with("rate_"), ~replace_na(., 0)),
    log_age = log(Age)
  )

write.csv(analysis_data, "maternal_behaviour_rate.csv")

# 7. EXPLORATORY DATA ANALYSIS
# ============================
cat("Analysis dataset contains", nrow(analysis_data), "mothers with infants < 1 year\n")

# Check distribution of offspring count
hist(analysis_data$offspring_count, 
     main = "Distribution of Offspring Count (Mothers with Infants < 1yr)", 
     xlab = "Offspring Count")

# Density plot
ggplot(analysis_data, aes(x = offspring_count)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  labs(title = "Distribution of Offspring Count", 
       x = "Offspring Count", y = "Density") +
  theme_minimal()

# Test for normality
shapiro_test <- shapiro.test(analysis_data$offspring_count)
cat("Shapiro-Wilk p-value:", round(shapiro_test$p.value, 4), "\n")

# Check for overdispersion with Poisson
poisson_test <- glmmTMB(offspring_count ~ rate_nursing_behaviors + offset(log_age), 
                        family = poisson, 
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

# 9. MATERNAL CARE MODEL ANALYSIS
# ===============================

# Full model with all 14 maternal care categories
full_maternal_model <- glmmTMB(
  offspring_count ~ rate_nursing_behaviors + rate_carrying_behaviors + rate_proximity_behaviors + 
    rate_active_grooming + rate_protective_interventions + rate_retrieving_offspring + 
    rate_soliciting_grooming + rate_comfort_contact + rate_food_sharing + 
    rate_stationary_touch + rate_muzzle_contact + rate_initiating_contact + 
    rate_reconciliation + rate_teeth_chatter + offset(log_age),
  family = poisson(),
  data = analysis_data
)

summary(full_maternal_model)

# 10. STEPWISE MODEL SELECTION
# ============================

# Function to perform stepwise elimination for glmmTMB
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

# Perform stepwise elimination with all 14 behavioral categories
maternal_vars <- c("rate_nursing_behaviors", "rate_carrying_behaviors", "rate_proximity_behaviors", 
                   "rate_active_grooming", "rate_protective_interventions", "rate_retrieving_offspring", 
                   "rate_soliciting_grooming", "rate_comfort_contact", "rate_food_sharing", 
                   "rate_stationary_touch", "rate_muzzle_contact", "rate_initiating_contact", 
                   "rate_reconciliation", "rate_teeth_chatter")

stepwise_result <- stepwise_elimination(analysis_data, maternal_vars, model_family)
final_maternal_model <- stepwise_result$final_model

cat("\n=== FINAL MATERNAL CARE MODEL ===\n")
summary(final_maternal_model)


# 13. SAVE RESULTS
# ================
write.csv(analysis_data, "maternal_care_analysis_data_updated.csv", row.names = FALSE)


print(behavior_summary)
write.csv(behavior_summary, "behavioral_category_summary.csv", row.names = FALSE)

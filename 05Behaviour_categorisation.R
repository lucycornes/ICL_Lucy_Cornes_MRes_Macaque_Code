
## ==============================================================================
# BEHAVIORAL CATEGORIZATION AND RATE CALCULATION
# Rhesus Macaque Maternal Effects on Same-Sex Sexual Behavior Study
# ==============================================================================

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(lubridate)

print("=== BEHAVIORAL CATEGORIZATION ANALYSIS ===")

# ==============================================================================
# 1. DATA IMPORT
# ==============================================================================

print("Loading cleaned focal data and demographic files...")

# Load cleaned focal data and demographic data
focal_data_scaled <- read_csv("focal_data_scaled.csv")
Offspring_Count <- read_csv("Offspring Count.csv")
Infant_count <- read_csv("Infant_count.csv")

print("Data loaded successfully")

# ==============================================================================
# 2. BASIC DATA PREPARATION
# ==============================================================================

print("=== PREPARING BEHAVIORAL DATA ===")

# Ensure Duration is numeric and clean whitespace
focal_data_scaled$Duration_s <- as.numeric(focal_data_scaled$Duration_s)

# Clean actor, receiver, and behavior names
focal_data_scaled <- focal_data_scaled %>%
  mutate(
    Actor = str_trim(Actor),
    Receiver = str_trim(Receiver),
    Behavior = str_trim(Behavior)
  )

print(paste("Total behaviors found:", length(unique(focal_data_scaled$Behavior))))

# ==============================================================================
# 3. CALCULATE BASIC BEHAVIORAL RATES (ALL BEHAVIORS)
# ==============================================================================

print("=== CALCULATING BASIC BEHAVIORAL RATES ===")

# Calculate duration of each behavior per Actor
behavior_totals <- focal_data_scaled %>%
  group_by(Actor, Behavior) %>%
  summarise(TotalBehaviorDuration_s = sum(Duration_s, na.rm = TRUE), .groups = "drop")

# Calculate total observation time per Actor
obs_time <- focal_data_scaled %>%
  distinct(Actor, Log_File_Name) %>%
  group_by(Actor) %>%
  summarise(NumFocals = n(), .groups = "drop") %>%
  mutate(TotalObsTime_s = NumFocals * 1800)  # 1800 seconds = 30 minutes per focal

# Join and calculate rate per hour
behavior_rates <- behavior_totals %>%
  left_join(obs_time, by = "Actor") %>%
  mutate(Rate_per_hour = (TotalBehaviorDuration_s / TotalObsTime_s) * 60)

# Reshape to wide format - one row per Actor, all behaviors as columns
behavior_wide <- behavior_rates %>%
  select(Actor, Behavior, Rate_per_hour) %>%
  pivot_wider(
    names_from = Behavior,
    values_from = Rate_per_hour,
    values_fill = 0  
  )

print(paste("Created behavioral rate matrix for", nrow(behavior_wide), "individuals"))

# ==============================================================================
# 4. ADD DEMOGRAPHIC DATA
# ==============================================================================

print("=== ADDING DEMOGRAPHIC DATA ===")

# Calculate total focal time per individual
actor_focal_time <- focal_data_scaled %>%
  group_by(Actor) %>%
  summarise(Total_seconds = sum(Duration_s, na.rm = TRUE), .groups = "drop") %>%
  mutate(Total_hours = Total_seconds / 3600)

# Join focal time data
Female_rate_behaviour <- behavior_wide %>%
  left_join(actor_focal_time, by = "Actor")

# Prepare demographic data from offspring count
Offspring_Count$AnimalID <- as.character(Offspring_Count$AnimalID)
Offspring_Count$AnimalID <- gsub("EE", "E", Offspring_Count$AnimalID)

today_year <- year(Sys.Date())
pedigree_data <- Offspring_Count %>%
  mutate(Age = today_year - BirthSeason)

# Join demographic data
Female_rate_behaviour <- Female_rate_behaviour %>%
  left_join(pedigree_data %>% select(AnimalID, Age, offspring_count), 
            by = c("Actor" = "AnimalID"))

print(paste("Added demographic data for", sum(!is.na(Female_rate_behaviour$Age)), "individuals"))

# ==============================================================================
# 5. DEFINE MATERNAL CARE BEHAVIORAL CATEGORIES
# ==============================================================================

print("=== DEFINING MATERNAL CARE CATEGORIES ===")

# Define behavioral categories using the original format from focal data
maternal_care_categories <- list(
  "Nursing" = c("Nursing", "N+G(give)", "N+G(receive)", "N+LA (G)", "N+LA(R)", "N+HA(g)"),
  "Carrying" = c("Carrying offspring", "C+HA(G)", "C+HA(R)", "C+LA(G)", "C+LA(R)"), 
  "Proximity" = c("Proximity to offspring", "P+G (give)", "P+G (receive)"),
  "Active_Grooming" = c("N+G(give)", "Groom (give)"),
  "Protective" = c("Protective interventions"),
  "Retrieving" = c("Retrieving offspring"),
  "Soliciting_Grooming" = c("Soliciting grooming"),
  "Comfort_Contact" = c("Hugging", "Huddling", "Clasped sleeping"),
  "Food_Sharing" = c("Food sharing (give)"),
  "Stationary_Touch" = c("Stationary touch"),
  "Muzzle_Contact" = c("Muzzle contact"),
  "Initiating_Contact" = c("Initiating contact"),
  "Reconciliation" = c("Reconciliation"),
  "Teeth_Chatter" = c("Teeth chatter")
)

# Define which categories need to be filtered to "Own Offspring" only
offspring_specific_categories <- c("Active_Grooming", "Comfort_Contact", "Food_Sharing", 
                                   "Stationary_Touch", "Muzzle_Contact")

print("Maternal care categories defined:")
for(category in names(maternal_care_categories)) {
  filter_type <- ifelse(category %in% offspring_specific_categories, 
                        "(filtered to own offspring)", 
                        "(all maternal behaviors)")
  print(paste("-", category, ":", length(maternal_care_categories[[category]]), 
              "behaviors", filter_type))
}

# ==============================================================================
# 6. CALCULATE MATERNAL CARE CATEGORY RATES (WITH PROPER FILTERING)
# ==============================================================================

print("=== CALCULATING MATERNAL CARE CATEGORY RATES ===")

# Create behavior durations dataframe for maternal care calculations
behavior_durations <- focal_data_scaled %>%
  select(Actor, Receiver, Behavior, Duration_s) %>%
  mutate(Duration_s = as.numeric(Duration_s))

# Calculate rates for your 14 specified maternal care categories
maternal_category_rates <- list()

for(category_name in names(maternal_care_categories)) {
  behaviors <- maternal_care_categories[[category_name]]
  
  # Determine if this category needs offspring filtering
  if(category_name %in% offspring_specific_categories) {
    # Filter to Own Offspring only for specific categories
    category_data <- focal_data_scaled %>%
      filter(Behavior %in% behaviors & Receiver == "Own Offspring") %>%
      group_by(Actor) %>%
      summarise(Total_seconds_category = sum(Duration_s, na.rm = TRUE), .groups = "drop")
    
    print(paste("Category", category_name, "filtered to own offspring"))
  } else {
    # No offspring filter for general maternal behaviors
    category_data <- focal_data_scaled %>%
      filter(Behavior %in% behaviors) %>%
      group_by(Actor) %>%
      summarise(Total_seconds_category = sum(Duration_s, na.rm = TRUE), .groups = "drop")
    
    print(paste("Category", category_name, "includes all maternal behaviors"))
  }
  
  # Calculate rates
  category_data <- category_data %>%
    left_join(Female_rate_behaviour %>% select(Actor, Total_seconds), by = "Actor") %>%
    mutate(
      rate_category = (Total_seconds_category / Total_seconds) * 60  # minutes per hour
    ) %>%
    select(Actor, rate_category)
  
  # Rename the rate column to match category
  category_rate_name <- paste0("rate_", tolower(gsub("_", "", category_name)))
  names(category_data)[names(category_data) == "rate_category"] <- category_rate_name
  
  maternal_category_rates[[category_name]] <- category_data
}

# Join all 14 maternal care category rates
maternal_rates <- Female_rate_behaviour %>%
  select(Actor, Total_seconds)

for(category_name in names(maternal_category_rates)) {
  maternal_rates <- maternal_rates %>%
    left_join(maternal_category_rates[[category_name]], by = "Actor")
}

# Replace NAs with 0 for mothers who didn't perform these behaviors
rate_columns <- grep("^rate_", names(maternal_rates), value = TRUE)
maternal_rates <- maternal_rates %>%
  mutate(across(all_of(rate_columns), ~replace_na(., 0)))

print(paste("Calculated rates for", length(rate_columns), "maternal care categories"))

# Summarize rates (in case of duplicates)
maternal_rates <- maternal_rates %>%
  group_by(Actor) %>%
  summarise(across(starts_with("rate_"), ~mean(.x, na.rm = TRUE)), .groups = "drop")

print(paste("Calculated maternal care rates for", nrow(maternal_rates), "individuals"))

# ==============================================================================
# 7. CREATE ANALYSIS DATASET WITH INFANT AGE FILTER
# ==============================================================================

print("=== CREATING FILTERED ANALYSIS DATASET ===")

# Create initial analysis dataset
analysis_data <- Female_rate_behaviour %>%
  select(Actor, offspring_count, Age) %>%
  left_join(maternal_rates, by = "Actor") %>%
  mutate(
    across(starts_with("rate_"), ~replace_na(., 0)),
    log_age = log(Age)
  )

# Add infant age information
analysis_data <- analysis_data %>%
  left_join(Infant_count %>% select(AnimalID, "Infant Age"),
            by = c("Actor" = "AnimalID"))

print(paste("Before infant age filter:", nrow(analysis_data), "individuals"))

# Filter for females with infants under 1 year old
analysis_data <- analysis_data %>%
  filter(`Infant Age` < 1)

print(paste("After infant age filter (<1 year):", nrow(analysis_data), "individuals"))

# Remove duplicates (in case of multiple infants under 1 year)
analysis_data <- analysis_data %>%
  distinct(Actor, .keep_all = TRUE)

print(paste("After removing duplicates:", nrow(analysis_data), "individuals"))

# ==============================================================================
# 8. SOCIAL CONNECTEDNESS CATEGORIES
# ==============================================================================

print("=== CALCULATING SOCIAL CONNECTEDNESS CATEGORIES ===")

# Define social connectedness categories using original format from focal data
social_connectedness_categories <- list(
  "Comfort_Contact" = c("Hugging", "Huddling", "Clasped sleeping"),
  "Proximity_Behaviors" = c("Stationary in proximity", "Approaching", "P+G (give)", 
                            "P+G (receive)", "Proximity to offspring", "Sit together"),
  "Grooming_Behaviors" = c("N+G(give)", "N+G(receive)", "Groom (give)", "Groom (receive)", 
                           "Soliciting grooming"),
  "Food_Sharing" = c("Food sharing (give)", "Food sharing (receive)"),
  "Muzzle_Contact" = c("Muzzle contact"),
  "Initiating_Contact" = c("Initiating contact"),
  "Passive_Body_Contact" = c("Passive body contact"),
  "Approaching" = c("Approaching"),
  "Stationary_Proximity" = c("Stationary in proximity"),
  "Teeth_Chatter" = c("Teeth chatter"),
  "Reconciliation" = c("Reconciliation"),
  "Lip_Smacking" = c("Lip smacking")
)

# Calculate social connectedness category rates
for(category_name in names(social_connectedness_categories)) {
  behaviors <- social_connectedness_categories[[category_name]]
  
  # Convert behavior names to column names (how they appear after pivot_wider)
  behavior_columns <- gsub("[^A-Za-z0-9]", ".", behaviors)
  behavior_columns <- gsub("\\.+", ".", behavior_columns)  # Replace multiple dots with single dot
  behavior_columns <- gsub("\\.$", "", behavior_columns)   # Remove trailing dot
  
  # Calculate rates for this social category
  column_name <- paste0("rate_social_", tolower(gsub("_", "", category_name)))
  
  Female_rate_behaviour[[column_name]] <- 0
  
  for(behavior_col in behavior_columns) {
    if(behavior_col %in% colnames(Female_rate_behaviour)) {
      Female_rate_behaviour[[column_name]] <- Female_rate_behaviour[[column_name]] + 
        replace_na(Female_rate_behaviour[[behavior_col]], 0)
    }
  }
}

print("Social connectedness categories calculated")

# ==============================================================================
# 10. FINAL OUTPUT AND SUMMARY
# ==============================================================================

print("=== SAVING RESULTS ===")

# Save the complete behavioral dataset
write.csv(Female_rate_behaviour, "Updated_female_beh_rate.csv", row.names = FALSE)

# Save the filtered analysis dataset
write.csv(analysis_data, "maternal_care_analysis_data.csv", row.names = FALSE)

print("=== BEHAVIORAL CATEGORIZATION COMPLETE ===")
print("Files created:")
print("- Updated_female_beh_rate.csv (complete behavioral dataset)")
print("- maternal_care_analysis_data.csv (filtered dataset for maternal care analysis)")

# Print summary statistics
print("\n=== SUMMARY STATISTICS ===")
print(paste("Total individuals in complete dataset:", nrow(Female_rate_behaviour)))
print(paste("Individuals with demographic data:", sum(!is.na(Female_rate_behaviour$Age))))
print(paste("Individuals in filtered analysis dataset:", nrow(analysis_data)))

print("\n=== FILTERING SUMMARY ===")
print("Categories filtered to 'Own Offspring' only:")
for(cat in offspring_specific_categories) {
  print(paste("-", cat))
}
print("Categories including ALL maternal behaviors:")
for(cat in names(maternal_care_categories)[!names(maternal_care_categories) %in% offspring_specific_categories]) {
  print(paste("-", cat))
}


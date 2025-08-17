# ==============================================================================
# BEHAVIORAL SUMMARY AND DEMOGRAPHIC DATA INTEGRATION
# Rhesus Macaque Maternal Effects on Same-Sex Sexual Behavior Study
# ==============================================================================

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

# ==============================================================================
# 1. DATA IMPORT
# ==============================================================================

print("=== IMPORTING DATA ===")

# Read in the scaled female focal data and anonymized demographic files
focal_data_scaled <- read_csv("New_All_Focals_Cleaned.csv")
Infant_count <- read_csv("Infant_count.csv")
Offspring_Count <- read_csv("Offspring Count.csv")

# Check column names to identify the correct infant age column
print("Infant count column names:")
print(colnames(Infant_count))
print("Offspring count column names:")
print(colnames(Offspring_Count))

print("Data import complete")

# ==============================================================================
# 2. BEHAVIORAL TIME SUMMARY
# ==============================================================================

print("=== CREATING BEHAVIORAL SUMMARY ===")

# Summarise data by FemaleID and Behaviour, summing the Duration_s
behavior_summary <- focal_data_scaled %>%
  group_by(Actor, Behavior) %>%
  summarize(Total_Time = sum(Duration_s, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Behavior, values_from = Total_Time, values_fill = list(Total_Time = 0))  # Pivot behaviors into columns

# This will give a table with one row per female and each behavior as a column 
print("Behavioral summary created:")
head(behavior_summary)

# Save initial behavioral summary
write.csv(behavior_summary, "behaviour_summary.csv", row.names = FALSE)

# ==============================================================================
# 3. DATA QUALITY CHECK - Ensuring all behaviours add up to 9000 seconds 
# ==============================================================================

print("=== PERFORMING DATA QUALITY CHECKS ===")

# Calculate total time for each female before pivoting (i.e., total time across all behaviors)
total_time_before <- focal_data_scaled %>%
  group_by(Actor) %>%
  summarize(Total_Time_All_Behaviors = sum(Duration_s, na.rm = TRUE), .groups = 'drop')

# Recreate behavior summary for checking (ensuring consistency)
behavior_summary_check <- focal_data_scaled %>%
  group_by(Actor, Behavior) %>%
  summarize(Total_Time = sum(Duration_s, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Behavior, values_from = Total_Time, values_fill = list(Total_Time = 0))

# Calculate total time for each female across all behaviors after pivoting (sum across columns)
# Note: Using all numeric columns instead of starts_with("Behavior") since behavior names vary
total_time_after <- behavior_summary_check %>%
  rowwise() %>%
  mutate(Total_Time_All_Behaviors_Pivoted = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Actor, Total_Time_All_Behaviors_Pivoted)

# Merge the two results to check if they match
check_total_time <- total_time_before %>%
  left_join(total_time_after, by = "Actor") %>%
  mutate(Totals_Match = abs(Total_Time_All_Behaviors - Total_Time_All_Behaviors_Pivoted) < 0.01)  # Allow for small rounding differences

# View the results
print("Data quality check results:")
print(paste("Total females:", nrow(check_total_time)))
print(paste("Females with matching totals:", sum(check_total_time$Totals_Match, na.rm = TRUE)))

if(any(!check_total_time$Totals_Match, na.rm = TRUE)) {
  print("WARNING: Some females have mismatched totals!")
  print(check_total_time[!check_total_time$Totals_Match, ])
} else {
  print("✓ All totals match - data quality check passed")
}

# ==============================================================================
# 4. ADDING DEMOGRAPHIC DATA
# ==============================================================================

print("=== INTEGRATING DEMOGRAPHIC DATA ===")

# Adding the age of each female to this table
# Ensure that the 'AnimalID' column is treated as a character type
Offspring_Count$AnimalID <- as.character(Offspring_Count$AnimalID)

# Note: The 'EE' to 'E' replacement may not be needed with anonymized IDs, 
# but keeping for compatibility
Offspring_Count$AnimalID <- gsub("EE", "E", Offspring_Count$AnimalID)

# Calculate current year for age calculation
today_year <- year(Sys.Date())  # Gets the current year, e.g., 2025

# Create an age column
pedigree_data <- Offspring_Count %>%
  mutate(Age = today_year - BirthSeason)

print("Pedigree data processed:")
print(paste("Total individuals in pedigree data:", nrow(pedigree_data)))

# Join Age and offspring count information to behavioral summary
female_behaviour <- behavior_summary %>%
  left_join(pedigree_data %>% select(AnimalID, Age, offspring_count), 
            by = c("Actor" = "AnimalID"))

print("After joining pedigree data:")
print(paste("Females with age data:", sum(!is.na(female_behaviour$Age))))

# Adding infant age to this table
# First check if infant age column exists and identify correct name
infant_age_col <- NULL
possible_names <- c("Infant Age", "Infant.Age", "infant_age", "InfantAge")

for(name in possible_names) {
  if(name %in% colnames(Infant_count)) {
    infant_age_col <- name
    break
  }
}

if(!is.null(infant_age_col)) {
  print(paste("Found infant age column:", infant_age_col))
  female_behaviour <- female_behaviour %>%
    left_join(Infant_count %>% select(AnimalID, all_of(infant_age_col)),
              by = c("Actor" = "AnimalID"))
  print("After joining infant data:")
  print(paste("Females with infant age data:", sum(!is.na(female_behaviour[[infant_age_col]]))))
} else {
  print("Warning: No infant age column found. Available columns in Infant_count:")
  print(colnames(Infant_count))
  print("Skipping infant age join...")
}

# ==============================================================================
# 5. FINAL OUTPUT
# ==============================================================================

# Save the complete behavioral summary with demographics
write.csv(female_behaviour, "female_behaviour.csv", row.names = FALSE)

print("=== BEHAVIORAL SUMMARY COMPLETE ===")
print("Files created:")
print("- behaviour_summary.csv (time spent on each behavior)")
print("- female_behaviour.csv (behavioral summary with demographics)")

# Final summary statistics
print("\nFinal dataset summary:")
print(paste("Total females in final dataset:", nrow(female_behaviour)))
print(paste("Females with complete demographic data:", sum(!is.na(female_behaviour$Age) & !is.na(female_behaviour$offspring_count))))
print(paste("Number of behavioral variables:", sum(sapply(female_behaviour, is.numeric)) - 3))  # Subtract Age, offspring_count, Infant Age

# ==============================================================================
# DATA CLEANING AND PREPARATION PIPELINE
# Rhesus Macaque Maternal Effects on Same-Sex Sexual Behavior Study
# ==============================================================================

# Load required libraries
library(dplyr)

# ==============================================================================
# 1. DATA IMPORT AND INITIAL CLEANING
# ==============================================================================

All_Focals_Combined <- read.csv("All_Focals_Combined.csv", 
                                header = TRUE, 
                                stringsAsFactors = FALSE,
                                check.names = FALSE)  # Prevents R from changing column names

# Check column names (should already be correct)
print("Column names:")
print(colnames(All_Focals_Combined))

# Remove X columns and empty columns first
print("Original column names:")
print(names(All_Focals_Combined))

# Remove X columns and columns that are entirely NA
All_Focals_Combined <- All_Focals_Combined %>%
  select(-starts_with("X")) %>%  # Remove all X, X.1, X.2, etc. columns
  select_if(~ !all(is.na(.)))     # Remove columns that are completely NA

# Ensure column names are strings
colnames(All_Focals_Combined) <- as.character(colnames(All_Focals_Combined))

# Remove any spaces before column names
colnames(All_Focals_Combined) <- trimws(colnames(All_Focals_Combined))

print("Cleaned column names:")
print(names(All_Focals_Combined))

# Remove columns that are not needed for analysis
New_All_Focals_Cleaned <- All_Focals_Combined %>%
  select(-'Time_Lag_s', -'Observation_Name', -'Observer_ID', -'Coding_Scheme')

# Check and remove quotation marks (") from ACTOR and RECIPIENT columns
New_All_Focals_Cleaned <- New_All_Focals_Cleaned %>%
  mutate(
    Actor = gsub('"', '', Actor),
    Receiver = gsub('"', '', Receiver),
    Behavior = gsub('"', '', Behavior)
  )

# Inspect the cleaned data
head(New_All_Focals_Cleaned)

# Save the initial cleaned data
write.csv(New_All_Focals_Cleaned, "New_All_Focals_Cleaned.csv", row.names = FALSE)

# ==============================================================================
# 2. FOCAL DURATION CALCULATIONS AND SCALING
# ==============================================================================

# Convert multiple columns to numeric for duration calculations
New_All_Focals_Cleaned <- New_All_Focals_Cleaned %>%
  mutate(
    Time_Relative_s = as.numeric(Time_Relative_s),
    Duration_s = as.numeric(Duration_s)
  )

# Calculate focal duration for each individual to check for scaling needs
focal_duration_check <- New_All_Focals_Cleaned %>%
  mutate(Duration_s = as.numeric(Duration_s)) %>%  # Ensure Duration_s is numeric
  group_by(Actor) %>%
  summarise(
    Total_Duration = round(sum(Duration_s, na.rm = TRUE)),
    Max_Time_Relative_s = max(Time_Relative_s, na.rm = TRUE)
  )

# Check the duration totals
print("Focal duration check:")
print(focal_duration_check)

# Save duration check results
write.csv(focal_duration_check, "focal_duration_check.csv", row.names = FALSE)

# Get actors with 7200 total duration (these need scaling to 9000 seconds)
actors_to_scale <- focal_duration_check$Actor[focal_duration_check$Total_Duration == 7200]

# Verify actors being scaled
original_totals <- New_All_Focals_Cleaned %>%
  filter(Actor %in% actors_to_scale) %>%
  group_by(Actor) %>%
  summarise(Total = sum(Duration_s))

print("Actors to be scaled (7200s to 9000s):")
print(original_totals)

# ==============================================================================
# 3. DURATION SCALING APPLICATION
# ==============================================================================

# Create a copy of focal data for scaling
focal_data_scaled <- New_All_Focals_Cleaned

# Scale Duration_s by 1.25 for matching actors only (7200 * 1.25 = 9000)
focal_data_scaled$Duration_s <- ifelse(
  focal_data_scaled$Actor %in% actors_to_scale,
  focal_data_scaled$Duration_s * 1.25,
  focal_data_scaled$Duration_s
)

# Verify scaling worked correctly
scaled_totals <- focal_data_scaled %>%
  filter(Actor %in% actors_to_scale) %>%
  group_by(Actor) %>%
  summarise(Total = sum(Duration_s))

print("After scaling - should now be 9000s:")
print(scaled_totals)

# ==============================================================================
# 4. FINAL OUTPUT
# ==============================================================================

# Save the final scaled dataset
write.csv(focal_data_scaled, "focal_data_scaled.csv", row.names = FALSE)

print("Data cleaning and scaling complete!")
print("Files created:")
print("- New_All_Focals_Cleaned.csv (initial cleaning)")
print("- focal_duration_check.csv (duration verification)")
print("- focal_data_scaled.csv (final scaled dataset)")

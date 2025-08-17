
library(glmmTMB)

# Load your dataset with unique MotherID values
mother_data <- read.csv("mother_master_sheet2025.csv", stringsAsFactors = FALSE)

# Step 1: Count how many sons are included per mother (based on rows in the dataset)
sons_included <- mother_data %>%
  group_by(MotherID) %>%
  summarise(n_sons_included = n())

# Step 2: Join this count to the main dataset (which already has TotalSons)
mother_data <- mother_data %>%
  left_join(sons_included, by = "MotherID") %>%
  mutate(
    prop_included = n_sons_included / TotalSons
  )

# Clean and prepare
mother_data <- mother_data %>%
  filter(
    !is.na(Total_Offspring_Count),
    !is.na(BirthSeason),
    !is.na(adjusted_SSB), !is.na(adjusted_DSB),
  ) %>%
  mutate(
    log_SSB = log(adjusted_SSB + 1),
    log_DSB = log(adjusted_DSB + 1),
    offset_birth = log(BirthSeason),
  )

# Fit the GLMM
glmm_model <- glmmTMB(
  Total_Offspring_Count ~ 
    log_SSB + 
    log_DSB + 
    offset(offset_birth) + 
    (1 | MotherID),
  weights = mother_data$prop_included,
  family = nbinom2,
  data = mother_data
)

library(ggplot2)

# View model summary
summary(glmm_model)

plot_SSB <- ggplot(mother_data, aes(x = log_SSB, y = Total_Offspring_Count)) +
  geom_point(color = "blue", alpha = 0.6) +  # Scatter plot for log adjusted SSB
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Regression line with confidence interval
  labs(title = "Significant Relationship: Total Offspring Count vs log Adjusted SSB",
       x = "Log Adjusted SSB", y = "Total Offspring Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


# Show the two plots side by side
plot_SSB 


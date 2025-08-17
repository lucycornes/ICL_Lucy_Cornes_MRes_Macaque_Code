# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(gridExtra)

# Read the CSV file (convert your Excel file to CSV first)
data <- read_csv("focal_count.csv")

# Function to extract hour from time slot strings
extract_start_hour <- function(time_slot) {
  if (is.na(time_slot)) return(NA)
  # Extract the start time from formats like "7:37 - 8:07" or "10:57-11:27"
  start_time <- gsub("\\s*-.*", "", time_slot)  # Remove everything after the dash
  start_time <- gsub("^\\s+|\\s+$", "", start_time)  # Trim whitespace
  
  # Extract hour
  hour <- as.numeric(gsub(":.*", "", start_time))
  return(hour)
}

# Create a long format dataset with all observation times
focal_times <- data %>%
  select(`Individual ID`, `Time Slot`, `Time slot2`, `Time slot3`, `Time slot4`, `Time slot5`) %>%
  pivot_longer(cols = starts_with("Time"), 
               names_to = "session", 
               values_to = "time_slot") %>%
  filter(!is.na(time_slot)) %>%
  mutate(start_hour = sapply(time_slot, extract_start_hour)) %>%
  filter(!is.na(start_hour)) %>%
  rename(Individual_ID = `Individual ID`)

# Get unique individuals and sort alphabetically
individuals <- sort(unique(focal_times$Individual_ID))

# Calculate y-axis limits to make them consistent across all plots
max_count <- focal_times %>%
  group_by(Individual_ID, start_hour) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pull(count) %>%
  max()

# Create individual histograms
plot_list <- list()

for (i in seq_along(individuals)) {
  individual_data <- focal_times %>% filter(Individual_ID == individuals[i])
  
  p <- ggplot(individual_data, aes(x = start_hour)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.7) +
    scale_x_continuous(breaks = 7:14, limits = c(6.5, 14.5)) +
    scale_y_continuous(limits = c(0, max_count)) +
    labs(title = individuals[i]) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text = element_text(size = 8),
          axis.title = element_blank())  # Remove individual axis titles
  
  plot_list[[i]] <- p
}

# Arrange plots in a grid (max 5 columns for Word document formatting)
n_individuals <- length(individuals)
ncol_grid <- 5  # Number of columns in the grid (max 5 for Word document)
nrow_grid <- ceiling(n_individuals / ncol_grid)

# Create the combined plot with shared axis labels
combined_plot <- grid.arrange(grobs = plot_list, ncol = ncol_grid,
                              bottom = "Start Hour",
                              left = "Number of Observations")

# Save the plot
ggsave("focal_follow_histograms.png", combined_plot, 
       width = 18, height = 3 * nrow_grid, dpi = 300)


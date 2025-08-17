
#####
# VIF Analysis for Behavioral Data
# =================================
# This script calculates Variance Inflation Factor (VIF) scores to assess
# multicollinearity in behavioral categories for research analysis.
#
# Author: [Your Name]
# Date: [Date]
# Dependencies: car, knitr

library(car)
library(knitr)

# ============================================================================
# CALCULATE VIF SCORES FOR BEHAVIORAL CATEGORIES
# ============================================================================

calculate_category_vif <- function(data, behavior_list, category_name) {
  # Calculate VIF scores for a category of behaviors
  #
  # Args:
  #   data: dataframe containing behavioral data
  #   behavior_list: vector of behavior column names
  #   category_name: string name for the category
  #
  # Returns:
  #   list with VIF analysis results
  
  # Get available behaviors in the dataset
  available <- intersect(behavior_list, names(data))
  
  # Handle single behavior case
  if(length(available) < 2) {
    return(list(
      category = category_name,
      behaviors = paste(available, collapse = ", "),
      vif_range = "N/A (single behavior)",
      cor_range = "N/A",
      n_behaviors = length(available),
      assessment = "Standalone",
      decision = "Use individually"
    ))
  }
  
  # Get clean data (complete cases only)
  clean_data <- data[available]
  clean_data <- clean_data[complete.cases(clean_data), ]
  
  # Check for sufficient data
  if(nrow(clean_data) < 10) {
    return(list(
      category = category_name,
      behaviors = paste(available, collapse = ", "),
      vif_range = "N/A (insufficient data)",
      cor_range = "N/A",
      n_behaviors = length(available),
      assessment = "Insufficient data",
      decision = "Cannot assess"
    ))
  }
  
  # Calculate correlation matrix and range
  cors <- cor(clean_data)
  cor_values <- cors[upper.tri(cors)]
  cor_range <- paste0(round(min(cor_values), 2), " to ", round(max(cor_values), 2))
  
  # Calculate VIF scores (requires 3+ behaviors)
  if(length(available) >= 3) {
    tryCatch({
      # Create linear model formula
      formula_str <- paste(available[1], "~", paste(available[-1], collapse = " + "))
      model <- lm(as.formula(formula_str), data = clean_data)
      vif_values <- vif(model)
      vif_range <- paste0(round(min(vif_values), 1), " - ", round(max(vif_values), 1))
      
      # Assess multicollinearity severity
      max_vif <- max(vif_values)
      if(max_vif > 10) {
        assessment <- "❌ High multicollinearity"
        decision <- "Remove redundant behaviors"
      } else if(max_vif > 5) {
        assessment <- "⚠️ Moderate concern"
        decision <- "Use with caution"
      } else {
        assessment <- "✅ Good"
        decision <- "Keep as composite"
      }
      
    }, error = function(e) {
      vif_range <- "Error calculating VIF"
      assessment <- "❌ Calculation failed"
      decision <- "Review behaviors"
    })
  } else {
    # Handle 2-behavior case using correlation
    if(abs(cor_values[1]) > 0.9) {
      vif_range <- "N/A (2 behaviors)"
      assessment <- "❌ Too highly correlated"
      decision <- "Use one behavior only"
    } else if(abs(cor_values[1]) > 0.3) {
      vif_range <- "N/A (2 behaviors)"
      assessment <- "✅ Good"
      decision <- "Keep as composite"
    } else {
      vif_range <- "N/A (2 behaviors)"
      assessment <- "⚠️ Weak correlation"
      decision <- "Consider using separately"
    }
  }
  
  return(list(
    category = category_name,
    behaviors = paste(available, collapse = ", "),
    vif_range = vif_range,
    cor_range = cor_range,
    n_behaviors = length(available),
    assessment = assessment,
    decision = decision
  ))
}

# ============================================================================
# BEHAVIORAL CATEGORY DEFINITIONS
# ============================================================================

# Define social connectedness categories using exact column names from dataset
social_connectedness_categories <- list(
  "Comfort_Contact" = c("Hugging", "Huddling", "Clasped.sleeping"),
  "Proximity_Behaviors" = c("Stationary.in.proximity", "Approaching", "P.G..give.", 
                            "P.G..receive.", "Proximity.to.offspring", "Sit.together"),
  "Grooming_Behaviors" = c("N.G.give.", "N.G.receive.", "Groom..give.", "Groom..receive.", 
                           "Soliciting.grooming"),
  "Food_Sharing" = c("Food.sharing..give.", "Food.sharing..receive."),
  "Muzzle_Contact" = c("Muzzle.contact"),
  "Initiating_Contact" = c("Initiating.contact"),
  "Passive_Body_Contact" = c("Passive.body.contact"),
  "Approaching" = c("Approaching"),
  "Stationary_Proximity" = c("Stationary.in.proximity"),
  "Teeth_Chatter" = c("Teeth.chatter"),
  "Reconciliation" = c("Reconciliation"),
  "Lip_Smacking" = c("Lip.smacking")
)

# Define maternal care categories using exact column names from dataset
maternal_care_categories <- list(
  "Nursing" = c("Nursing", "N.G.give.", "N.G.receive.", "N.LA..G.", "N.LA.R.", "N.HA.g."),
  "Carrying" = c("Carrying.offspring", "C.HA.G.", "C.HA.R.", "C.LA.G.", "C.LA.R."),
  "Proximity" = c("Proximity.to.offspring", "P.G..give.", "P.G..receive."),
  "Active_Grooming" = c("N.G.give.", "Groom..give."),
  "Protective_Interventions" = c("Protective.interventions"),
  "Retrieving_Offspring" = c("Retrieving.offspring"),
  "Comfort_Contact" = c("Hugging", "Huddling", "Clasped.sleeping"),
  "Food_Sharing" = c("Food.sharing..give."),
  "Stationary_Touch" = c("Stationary.touch"),
  "Muzzle_Contact" = c("Muzzle.contact")
)

# ============================================================================
# VIF TABLE GENERATION FUNCTION
# ============================================================================

generate_vif_table <- function(data, categories, table_title) {
  # Generate formatted VIF analysis table
  #
  # Args:
  #   data: behavioral dataset
  #   categories: list of behavioral categories
  #   table_title: title for the output table
  #
  # Returns:
  #   dataframe with VIF analysis results
  
  # Calculate VIF for all categories
  results <- list()
  for(name in names(categories)) {
    results[[name]] <- calculate_category_vif(data, categories[[name]], name)
  }
  
  # Convert results to data frame
  table_data <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      Category = x$category,
      Behaviors = x$behaviors,
      VIF_Range = x$vif_range,
      Correlation_Range = x$cor_range,
      Assessment = x$assessment,
      Decision = x$decision,
      stringsAsFactors = FALSE
    )
  }))
  
  # Print formatted markdown table
  cat("\n# ", table_title, "\n\n")
  cat("| Category | Behaviors Included | VIF Range | Correlation Range | Assessment | Decision |\n")
  cat("|----------|-------------------|-----------|-------------------|------------|----------|\n")
  
  for(i in 1:nrow(table_data)) {
    cat(sprintf("| %s | %s | %s | %s | %s | %s |\n",
                table_data$Category[i],
                table_data$Behaviors[i],
                table_data$VIF_Range[i],
                table_data$Correlation_Range[i],
                table_data$Assessment[i],
                table_data$Decision[i]))
  }
  
  # Add interpretation notes
  cat("\n**Notes:**\n")
  cat("- VIF Interpretation: <5 = acceptable, 5-10 = moderate concern, >10 = problematic\n")
  cat("- ✅ = Good for analysis, ⚠️ = Use with caution, ❌ = Problematic\n")
  cat("- Sample size: n =", nrow(data[complete.cases(data)]), "\n\n")
  
  return(table_data)
}

# ============================================================================
# FORMATTED TABLE EXPORT FUNCTION
# ============================================================================

create_formatted_table <- function(table_data, title, filename) {
  # Export formatted table to text file
  #
  # Args:
  #   table_data: dataframe with VIF results
  #   title: table title
  #   filename: output filename
  
  # Write formatted table to file
  sink(filename)
  
  cat(title, "\n")
  cat(rep("=", nchar(title)), "\n\n")
  
  # Create markdown table
  cat("| Category | Behaviors Included | VIF Range | Correlation Range | Assessment | Decision |\n")
  cat("|----------|-------------------|-----------|-------------------|------------|----------|\n")
  
  for(i in 1:nrow(table_data)) {
    cat(sprintf("| %s | %s | %s | %s | %s | %s |\n",
                table_data$Category[i],
                table_data$Behaviors[i],
                table_data$VIF_Range[i],
                table_data$Correlation_Range[i],
                table_data$Assessment[i],
                table_data$Decision[i]))
  }
  
  cat("\n**Notes:**\n")
  cat("- VIF Interpretation: <5 = acceptable, 5-10 = moderate concern, >10 = problematic\n")
  cat("- ✅ = Good for analysis, ⚠️ = Use with caution, ❌ = Problematic\n")
  
  # Add summary statistics
  composites <- sum(grepl("composite|Keep as composite", table_data$Decision))
  individuals <- sum(grepl("individually", table_data$Decision))
  cat("\n**Final Variables:**\n")
  cat("- Composite variables:", composites, "\n")
  cat("- Individual variables:", individuals, "\n")
  cat("- Total variables:", nrow(table_data), "\n")
  
  sink()
  cat("✅ Formatted table saved to:", filename, "\n")
}

# ============================================================================
# MAIN ANALYSIS EXECUTION
# ============================================================================

# Load behavioral data
cat("Loading behavioral data...\n")
df <- read.csv("Updated_female_beh_rate.csv")  # Update filename as needed

cat("Generating VIF tables with calculated values...\n\n")

# Generate VIF analysis tables
social_table <- generate_vif_table(df, social_connectedness_categories, 
                                   "Table 1: Multicollinearity Assessment of Social Connectedness Variables")

maternal_table <- generate_vif_table(df, maternal_care_categories, 
                                     "Table 2: Multicollinearity Assessment of Maternal Care Variables")

# Print summary statistics
cat("## Analysis Summary\n")
cat("**Social Connectedness Variables:**\n")
social_composites <- sum(grepl("composite|Keep as composite", social_table$Decision))
social_individual <- sum(grepl("individually", social_table$Decision))
cat("- Composite variables:", social_composites, "\n")
cat("- Individual variables:", social_individual, "\n")

cat("\n**Maternal Care Variables:**\n") 
maternal_composites <- sum(grepl("composite|Keep as composite", maternal_table$Decision))
maternal_individual <- sum(grepl("individually", maternal_table$Decision))
cat("- Composite variables:", maternal_composites, "\n")
cat("- Individual variables:", maternal_individual, "\n")

# ============================================================================
# EXPORT RESULTS
# ============================================================================

# Save formatted tables to text files
create_formatted_table(social_table, 
                       "Table 1: Multicollinearity Assessment of Social Connectedness Variables",
                       "social_connectedness_vif_table.txt")

create_formatted_table(maternal_table,
                       "Table 2: Multicollinearity Assessment of Maternal Care Variables", 
                       "maternal_care_vif_table.txt")

# Export to CSV for further analysis
write.csv(social_table, "social_connectedness_vif_results.csv", row.names = FALSE)
write.csv(maternal_table, "maternal_care_vif_results.csv", row.names = FALSE)

# Create comprehensive summary report
sink("vif_analysis_summary_report.txt")
cat("VIF ANALYSIS SUMMARY REPORT\n")
cat("============================\n\n")
cat("Analysis Date:", Sys.Date(), "\n")
cat("Dataset: Updated_female_beh_rate.csv\n")
cat("Sample Size:", nrow(df[complete.cases(df)]), "complete cases\n\n")

cat("METHODOLOGY:\n")
cat("- VIF scores calculated using linear regression models\n")
cat("- Categories with <2 behaviors assessed individually\n")
cat("- Correlation analysis for 2-behavior categories\n")
cat("- Multicollinearity thresholds: VIF >5 (moderate), VIF >10 (severe)\n\n")

cat("RESULTS SUMMARY:\n")
cat("Social Connectedness: ", nrow(social_table), " categories analyzed\n")
cat("Maternal Care: ", nrow(maternal_table), " categories analyzed\n\n")

cat("For detailed results, see individual table files.\n")
sink()

cat("\n✅ Analysis complete! Files generated:\n")
cat("- social_connectedness_vif_table.txt\n")
cat("- maternal_care_vif_table.txt\n")
cat("- social_connectedness_vif_results.csv\n")
cat("- maternal_care_vif_results.csv\n")
cat("- vif_analysis_summary_report.txt\n")
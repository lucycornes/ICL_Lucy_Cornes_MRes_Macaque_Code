# Load required libraries
library(tidyverse)
library(corrplot)
library(factoextra)

# Read the data
social_data <- read.csv("social_connectedness.csv")

# View the structure of the data
str(social_data)
head(social_data)

# =============================================================================
# IDENTIFY SOCIAL BEHAVIOR COLUMNS AND CHECK VARIATION
# =============================================================================

# Get all rate columns from the social data
rate_columns <- names(social_data)[grepl("^rate_", names(social_data))]
cat("Found", length(rate_columns), "social behavior rate columns:\n")
cat(paste(rate_columns, collapse = ", "), "\n")

# Check which behaviors have variation and sufficient data
behavior_summary <- social_data[, rate_columns] %>%
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

cat("\nSocial behavior summary:\n")
print(behavior_summary)

# Keep only behaviors with variation
good_behaviors <- behavior_summary$behavior[behavior_summary$has_variation]
cat("\nBehaviors with sufficient variation (", length(good_behaviors), "):\n")
cat(paste(good_behaviors, collapse = ", "), "\n")

if (length(good_behaviors) == 0) {
  stop("No social behaviors have sufficient variation for analysis")
}

# =============================================================================
# SELECT BEHAVIORS FOR PCA
# =============================================================================

# Select only the social behavior rate variables with sufficient variation for PCA
behavioral_vars <- social_data[, good_behaviors]

# Check for missing values
sum(is.na(behavioral_vars))

# Remove any rows with missing values if present
behavioral_vars_clean <- na.omit(behavioral_vars)

# Check correlations between variables (PCA works best when variables are correlated)
correlation_matrix <- cor(behavioral_vars_clean)
print(correlation_matrix)

# Visualize correlation matrix
corrplot(correlation_matrix, method = "color", type = "upper", 
         order = "hclust", tl.cex = 0.8, tl.col = "black")

# Perform PCA (scale = FALSE since all rates are on same scale, center = TRUE)
pca_result <- prcomp(behavioral_vars_clean, scale = FALSE, center = TRUE)

# Summary of PCA results
summary(pca_result)

# Scree plot to visualize explained variance
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

# Biplot showing variables and individuals
fviz_pca_biplot(pca_result, 
                col.var = "contrib", # Color variables by their contribution
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE, # Avoid text overlapping
                title = "Social Connectedness PCA Biplot")

# Variable contributions to PC1 and PC2
fviz_contrib(pca_result, choice = "var", axes = 1, top = 10,
             title = "Variable Contributions to PC1")

fviz_contrib(pca_result, choice = "var", axes = 2, top = 10,
             title = "Variable Contributions to PC2")

# Extract PC1 scores as the "social connectedness" score
# PC1 typically captures the most variation and often represents overall behavior intensity
social_connectedness_scores <- pca_result$x[, 1]

# Create a dataframe with Actor IDs and their social connectedness scores
# Note: Make sure the row order matches between original data and PCA input
social_scores_df <- data.frame(
  Actor = social_data$Actor[complete.cases(behavioral_vars)], # Only keep actors with complete data
  Social_Connectedness_Score = social_connectedness_scores,
  PC2_Score = pca_result$x[, 2]  # Include PC2 in case it's also meaningful
)

# Add back other variables you might need for your models
final_data <- social_data[complete.cases(behavioral_vars), ] %>%
  mutate(Social_Connectedness_Score = social_connectedness_scores,
         PC2_Score = pca_result$x[, 2])

# View the results
print(social_scores_df)

# Summary statistics of the social connectedness scores
summary(social_connectedness_scores)

# Plot distribution of social connectedness scores
ggplot(social_scores_df, aes(x = Social_Connectedness_Score)) +
  geom_histogram(bins = 8, fill = "lightgreen", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(density) * nrow(social_scores_df) * 0.5), 
               color = "red", size = 1) +
  labs(title = "Distribution of Social Connectedness Scores (PC1)",
       x = "Social Connectedness Score", y = "Frequency") +
  theme_minimal()

# Check which behaviors load most strongly on PC1 (your social connectedness score)
pc1_loadings <- pca_result$rotation[, 1]
pc1_loadings_sorted <- sort(abs(pc1_loadings), decreasing = TRUE)
print("Variables with strongest loadings on PC1 (Social Connectedness Score):")
print(pc1_loadings_sorted)

# Save the results
write.csv(final_data, "social_data_with_PCA_scores.csv", row.names = FALSE)

# Load required libraries
library(tidyverse)
library(MASS)      # for glm.nb() negative binomial - should be pre-installed

# Read the data
social_data <- read.csv("social_data_with_PCA_scores.csv")

# Examine the data structure
str(social_data)
summary(social_data)

# Examine offspring count distribution
table(social_data$offspring_count)
hist(social_data$offspring_count, 
     main = "Distribution of Offspring Count", 
     xlab = "Offspring Count", 
     col = "lightblue")

# Basic descriptive statistics for offspring count
mean_offspring <- mean(social_data$offspring_count)
var_offspring <- var(social_data$offspring_count)
cat("Mean:", mean_offspring, "\n")
cat("Variance:", var_offspring, "\n")
cat("Variance/Mean ratio:", var_offspring/mean_offspring, "\n")

# If variance/mean ratio >> 1, we likely have overdispersion
if(var_offspring/mean_offspring > 2) {
  cat("Variance much larger than mean - suggests overdispersion\n")
} else {
  cat("Variance similar to mean - Poisson might be appropriate\n")
}

# Examine the relationship visually
plot(social_data$Social_Connectedness_Score, social_data$offspring_count,
     xlab = "Social Connectedness Score", ylab = "Offspring Count",
     main = "Social Connectedness Score vs Offspring Count")

# Check for missing values
sum(is.na(social_data))


# Load required libraries
library(tidyverse)

# Read the datasets
social_data <- read.csv("social_data_with_PCA_scores.csv")
davidscore_data <- read.csv("Female_DS_byGroup.csv")

# Check the structure of both datasets
cat("Social data structure:\n")
str(social_data)
cat("\nDavidScore data structure:\n")
str(davidscore_data)

# Check unique actors/animals to understand the linking
cat("\nUnique actors in social data:", length(unique(social_data$Actor)), "\n")
cat("Unique animals in DavidScore data:", length(unique(davidscore_data$AnimalID)), "\n")

# Preview the linking variables
cat("\nSample Actor IDs from social data:\n")
print(head(social_data$Actor, 10))
cat("\nSample AnimalIDs from DavidScore data:\n")
print(head(davidscore_data$AnimalID, 10))

# Check for any mismatches between Actor and AnimalID
actors_not_in_davidscore <- setdiff(social_data$Actor, davidscore_data$AnimalID)
davidscore_not_in_social <- setdiff(davidscore_data$AnimalID, social_data$Actor)

cat("\nActors in social data but not in DavidScore data:", length(actors_not_in_davidscore), "\n")
if(length(actors_not_in_davidscore) > 0) {
  cat("Missing actors:", paste(actors_not_in_davidscore, collapse = ", "), "\n")
}

cat("\nAnimalIDs in DavidScore data but not in social data:", length(davidscore_not_in_social), "\n")
if(length(davidscore_not_in_social) > 0) {
  cat("Extra animals:", paste(davidscore_not_in_social, collapse = ", "), "\n")
}

# Standardize social data IDs to EE format
social_data <- social_data %>%
  mutate(Actor = as.character(Actor),
         Actor = gsub("^(\\d+)E(\\d+)$", "\\1EE\\2", Actor))

# Standardize DavidScore data IDs to EE format  
davidscore_data <- davidscore_data %>%
  mutate(AnimalID = as.character(AnimalID),
         AnimalID = gsub("^(\\d+)E(\\d+)$", "\\1EE\\2", AnimalID))



# Merge the datasets
social_data_with_davidscore <- social_data %>%
  left_join(davidscore_data[, c("AnimalID", "DavidScore", "Group", "Rank")], 
            by = c("Actor" = "AnimalID"))

# Check the merge results
cat("\nMerge results:\n")
cat("Original social data rows:", nrow(social_data), "\n")
cat("Final merged data rows:", nrow(social_data_with_davidscore), "\n")
cat("Rows with DavidScore:", sum(!is.na(social_data_with_davidscore$DavidScore)), "\n")
cat("Rows missing DavidScore:", sum(is.na(social_data_with_davidscore$DavidScore)), "\n")

# Check the new variables
cat("\nSummary of added variables:\n")
cat("DavidScore summary:\n")
print(summary(social_data_with_davidscore$DavidScore))
cat("\nGroup distribution:\n")
print(table(social_data_with_davidscore$Group, useNA = "ifany"))
cat("\nRank summary:\n")
print(summary(social_data_with_davidscore$Rank))

# Check for any issues with the merge
cat("\nChecking for duplicate actors after merge:\n")
duplicate_actors <- social_data_with_davidscore$Actor[duplicated(social_data_with_davidscore$Actor)]
if(length(duplicate_actors) > 0) {
  cat("WARNING: Duplicate actors found:", paste(duplicate_actors, collapse = ", "), "\n")
} else {
  cat("No duplicate actors - merge looks good!\n")
}

# Display the structure of the final dataset
cat("\nFinal dataset structure:\n")
str(social_data_with_davidscore)

# Preview the merged data
cat("\nFirst few rows of merged data:\n")
print(head(social_data_with_davidscore %>% 
             select(Actor, Social_Connectedness_Score, DavidScore, Group, Rank), 10))

# Save the merged dataset
write.csv(social_data_with_davidscore, "social_data_with_PCA_and_DavidScore.csv", row.names = FALSE)
cat("\nMerged dataset saved as 'social_data_with_PCA_and_DavidScore.csv'\n")

# =============================================================================
# MODEL FITTING AND FAMILY SELECTION
# =============================================================================

social_DS_PCA <- read.csv("social_data_with_PCA_and_DavidScore.csv")

# Start with Poisson regression (most restrictive assumption)
model_poisson <- glm(offspring_count ~ Social_Connectedness_Score + DavidScore + offset(log_age), 
                     data = social_DS_PCA, 
                     family = poisson(link = "log"))

summary(model_poisson)

# Create the negative binomial model to match your Poisson model
library(MASS)  # Make sure MASS is loaded for glm.nb

model_nb <- glm.nb(offspring_count ~ Social_Connectedness_Score + DavidScore + offset(log_age), 
                   data = social_DS_PCA)

# Then you can compare them fairly
summary(model_nb)
AIC(model_poisson)
AIC(model_nb)

# Manual overdispersion test for Poisson
# Calculate Pearson chi-square statistic
pearson_resid <- residuals(model_poisson, type = "pearson")
pearson_chisq <- sum(pearson_resid^2)
df_residual <- df.residual(model_poisson)
dispersion_ratio <- pearson_chisq / df_residual

cat("\n=== OVERDISPERSION TEST ===\n")
cat("Pearson chi-square:", round(pearson_chisq, 2), "\n")
cat("Degrees of freedom:", df_residual, "\n")
cat("Dispersion ratio:", round(dispersion_ratio, 3), "\n")

if(dispersion_ratio > 1.5) {
  cat("Dispersion ratio > 1.5 suggests overdispersion - try negative binomial\n")
} else {
  cat("Dispersion ratio close to 1 - Poisson may be appropriate\n")
}

# Basic residual plots for Poisson
par(mfrow = c(2, 2))
plot(model_poisson)
par(mfrow = c(1, 1))


# =============================================================================
# MODEL COMPARISON
# =============================================================================

# Compare Poisson vs Negative Binomial using AIC
aic_poisson <- AIC(model_poisson)
aic_nb <- AIC(model_nb)

cat("\n=== MODEL COMPARISON ===\n")
cat("Poisson AIC:", round(aic_poisson, 2), "\n")
cat("Negative Binomial AIC:", round(aic_nb, 2), "\n")
cat("AIC difference:", round(aic_poisson - aic_nb, 2), "\n")

if(aic_nb < aic_poisson - 2) {
  cat("Negative binomial is substantially better (ΔAIC > 2)\n")
  best_model <- "negative_binomial"
} else if(aic_poisson < aic_nb - 2) {
  cat("Poisson is substantially better (ΔAIC > 2)\n")
  best_model <- "poisson"
} else {
  cat("Models are similar - choose based on biological reasoning\n")
  best_model <- "similar"
}

# Likelihood ratio test (manual)
loglik_poisson <- logLik(model_poisson)
loglik_nb <- logLik(model_nb)
lr_statistic <- 2 * (loglik_nb - loglik_poisson)
p_value_lr <- 1 - pchisq(lr_statistic, df = 1)

cat("\nLikelihood Ratio Test:\n")
cat("LR statistic:", round(lr_statistic, 3), "\n")
cat("P-value:", round(p_value_lr, 4), "\n")

if(p_value_lr < 0.05) {
  cat("Negative binomial significantly better than Poisson (p < 0.05)\n")
}

# =============================================================================
# CHOOSE FINAL MODEL
# =============================================================================

# Select final model based on tests
if(best_model == "negative_binomial" || p_value_lr < 0.05 || dispersion_ratio > 1.5) {
  final_model <- model_nb
  model_type <- "Negative Binomial"
} else {
  final_model <- model_poisson
  model_type <- "Poisson"
}

cat("\n=== FINAL MODEL SELECTED ===\n")
cat("Using:", model_type, "regression\n")

# =============================================================================
# FINAL MODEL RESULTS AND INTERPRETATION
# =============================================================================

cat("\n=== FINAL MODEL SUMMARY ===\n")
summary(final_model)

# Calculate confidence intervals (manual method)
coef_estimates <- coef(final_model)
se_estimates <- summary(final_model)$coefficients[, "Std. Error"]

# 95% confidence intervals
ci_lower <- coef_estimates - 1.96 * se_estimates
ci_upper <- coef_estimates + 1.96 * se_estimates

cat("\n=== CONFIDENCE INTERVALS ===\n")
for(i in 1:length(coef_estimates)) {
  cat(names(coef_estimates)[i], ":", 
      round(ci_lower[i], 4), "to", round(ci_upper[i], 4), "\n")
}

# Effect size interpretation
coef_social <- coef(final_model)["Social_Connectedness_Score"]
effect_size <- exp(coef_social)

cat("\n=== BIOLOGICAL INTERPRETATION ===\n")
cat("For each 1-unit increase in Social Connectedness Score:\n")
cat("Offspring count multiplied by", round(effect_size, 3), "\n")
cat("This represents a", round((effect_size - 1) * 100, 1), "% change in offspring count\n")





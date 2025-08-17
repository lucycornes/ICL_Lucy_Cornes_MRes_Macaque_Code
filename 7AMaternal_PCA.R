# Load required libraries
install.packages(c("tidyverse", "corrplot", "factoextra"))
library(tidyverse)
library(corrplot)
library(factoextra)
library(MASS)

# Read the data

maternal_data <- read.csv("maternal_behaviour_rate.csv")

# View the structure of the data
str(maternal_data)
head(maternal_data)

# Select only the maternal behavior rate variables for PCA
behavioral_vars <- maternal_data %>%
  dplyr::select(rate_nursing_behaviors, rate_carrying_behaviors, rate_proximity_behaviors,
                rate_active_grooming, rate_protective_interventions, rate_retrieving_offspring,
                rate_comfort_contact, rate_food_sharing, rate_stationary_touch, rate_muzzle_contact)

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

# Perform PCA (scale = TRUE standardizes the variables)
pca_result <- prcomp(behavioral_vars_clean, scale = TRUE, center = TRUE)

# Summary of PCA results
summary(pca_result)

# Scree plot to visualize explained variance
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

# Biplot showing variables and individuals
fviz_pca_biplot(pca_result, 
                col.var = "contrib", # Color variables by their contribution
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE, # Avoid text overlapping
                title = "Maternal Behavior PCA Biplot")

# Variable contributions to PC1 and PC2
fviz_contrib(pca_result, choice = "var", axes = 1, top = 10,
             title = "Variable Contributions to PC1")

fviz_contrib(pca_result, choice = "var", axes = 2, top = 10,
             title = "Variable Contributions to PC2")

# Extract PC1 scores as the "maternality" score
# PC1 typically captures the most variation and often represents overall behavior intensity
maternality_scores <- pca_result$x[, 1]

# Create a dataframe with Actor IDs and their maternality scores
# Note: Make sure the row order matches between original data and PCA input
maternal_scores_df <- data.frame(
  Actor = maternal_data$Actor[complete.cases(behavioral_vars)], # Only keep actors with complete data
  Maternality_Score = maternality_scores,
  PC2_Score = pca_result$x[, 2]  # Include PC2 in case it's also meaningful
)

# Add back other variables you might need for your models
final_data <- maternal_data[complete.cases(behavioral_vars), ] %>%
  mutate(Maternality_Score = maternality_scores,
         PC2_Score = pca_result$x[, 2])

# View the results
print(maternal_scores_df)

# Summary statistics of the maternality scores
summary(maternality_scores)

# Plot distribution of maternality scores
ggplot(maternal_scores_df, aes(x = Maternality_Score)) +
  geom_histogram(bins = 8, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(density) * nrow(maternal_scores_df) * 0.5), 
               color = "red", size = 1) +
  labs(title = "Distribution of Maternality Scores (PC1)",
       x = "Maternality Score", y = "Frequency") +
  theme_minimal()

# Check which behaviors load most strongly on PC1 (your maternality score)
pc1_loadings <- pca_result$rotation[, 1]
pc1_loadings_sorted <- sort(abs(pc1_loadings), decreasing = TRUE)
print("Variables with strongest loadings on PC1 (Maternality Score):")
print(pc1_loadings_sorted)

# Save the results
write.csv(final_data, "new_maternal_data_with_PCA_scores.csv", row.names = FALSE)



######
# Load required libraries
library(tidyverse)

# Read the datasets
maternal_data <- read.csv("new_maternal_data_with_PCA_scores.csv")
davidscore_data <- read.csv("Female_DS_byGroup.csv")

# STANDARDIZE ID FORMATS TO MATCH (EE format is Excel-safe)
cat("=== STANDARDIZING ID FORMATS ===\n")

# Check original formats
cat("Original Actor IDs sample:\n")
print(head(maternal_data$Actor, 10))
cat("Original DavidScore AnimalIDs sample:\n")
print(head(davidscore_data$AnimalID, 10))

# Standardize maternal data IDs to EE format
maternal_data <- maternal_data %>%
  mutate(Actor = as.character(Actor),
         Actor = gsub("^(\\d+)E(\\d+)$", "\\1EE\\2", Actor))

# Standardize DavidScore data IDs to EE format  
davidscore_data <- davidscore_data %>%
  mutate(AnimalID = as.character(AnimalID),
         AnimalID = gsub("^(\\d+)E(\\d+)$", "\\1EE\\2", AnimalID))

# Check standardized formats
cat("\nStandardized Actor IDs sample:\n")
print(head(maternal_data$Actor, 10))
cat("Standardized DavidScore AnimalIDs sample:\n")
print(head(davidscore_data$AnimalID, 10))

# Check the structure of both datasets
cat("\nMaternal data structure:\n")
str(maternal_data)
cat("\nDavidScore data structure:\n")
str(davidscore_data)

# Check unique actors/animals to understand the linking
cat("\nUnique actors in maternal data:", length(unique(maternal_data$Actor)), "\n")
cat("Unique animals in DavidScore data:", length(unique(davidscore_data$AnimalID)), "\n")

# Check for any mismatches between Actor and AnimalID
actors_not_in_davidscore <- setdiff(maternal_data$Actor, davidscore_data$AnimalID)
davidscore_not_in_maternal <- setdiff(davidscore_data$AnimalID, maternal_data$Actor)

cat("\nActors in maternal data but not in DavidScore data:", length(actors_not_in_davidscore), "\n")
if(length(actors_not_in_davidscore) > 0) {
  cat("Missing actors:", paste(actors_not_in_davidscore, collapse = ", "), "\n")
}

cat("\nAnimalIDs in DavidScore data but not in maternal data:", length(davidscore_not_in_maternal), "\n")
if(length(davidscore_not_in_maternal) > 0) {
  cat("Extra animals:", paste(davidscore_not_in_maternal, collapse = ", "), "\n")
}

# Merge the datasets
# Use left_join to keep all maternal data and add DavidScore where available
# Use base R subsetting instead of select()
maternal_data_with_davidscore <- maternal_data %>%
  left_join(davidscore_data[, c("AnimalID", "DavidScore", "Group", "Rank")], 
            by = c("Actor" = "AnimalID"))

# Check the merge results
cat("\nMerge results:\n")
cat("Original maternal data rows:", nrow(maternal_data), "\n")
cat("Final merged data rows:", nrow(maternal_data_with_davidscore), "\n")
cat("Rows with DavidScore:", sum(!is.na(maternal_data_with_davidscore$DavidScore)), "\n")
cat("Rows missing DavidScore:", sum(is.na(maternal_data_with_davidscore$DavidScore)), "\n")

# Check the new variables
cat("\nSummary of added variables:\n")
cat("DavidScore summary:\n")
print(summary(maternal_data_with_davidscore$DavidScore))
cat("\nGroup distribution:\n")
print(table(maternal_data_with_davidscore$Group, useNA = "ifany"))
cat("\nRank summary:\n")
print(summary(maternal_data_with_davidscore$Rank))

# Check for any issues with the merge
cat("\nChecking for duplicate actors after merge:\n")
duplicate_actors <- maternal_data_with_davidscore$Actor[duplicated(maternal_data_with_davidscore$Actor)]
if(length(duplicate_actors) > 0) {
  cat("WARNING: Duplicate actors found:", paste(duplicate_actors, collapse = ", "), "\n")
} else {
  cat("No duplicate actors - merge looks good!\n")
}

# Display the structure of the final dataset
cat("\nFinal dataset structure:\n")
str(maternal_data_with_davidscore)

# Preview the merged data
cat("\nFirst few rows of merged data:\n")
print(head(maternal_data_with_davidscore %>% 
             dplyr::select(Actor, Maternality_Score, DavidScore, Group, Rank), 10))

# Save the merged dataset
write.csv(maternal_data_with_davidscore, "new_maternal_data_with_PCA_and_DavidScore.csv", row.names = FALSE)
cat("\nMerged dataset saved as 'new_maternal_data_with_PCA_and_DavidScore.csv'\n")



# Load required libraries (using more common packages)
library(tidyverse)
library(MASS)      # for glm.nb() negative binomial - should be pre-installed

# Read the data
maternal_data <- read.csv("maternal_data_with_PCA_scores.csv")

# Examine the data structure
str(maternal_data)
summary(maternal_data)

# Examine offspring count distribution
table(maternal_data$offspring_count)
hist(maternal_data$offspring_count, 
     main = "Distribution of Offspring Count", 
     xlab = "Offspring Count", 
     col = "lightblue")

# Basic descriptive statistics for offspring count
mean_offspring <- mean(maternal_data$offspring_count)
var_offspring <- var(maternal_data$offspring_count)
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
plot(maternal_data$Maternality_Score, maternal_data$offspring_count,
     xlab = "Maternality Score", ylab = "Offspring Count",
     main = "Offspring Count vs Maternality Score")

# Check for missing values
sum(is.na(maternal_data))
view(maternal_data)
# =============================================================================
# MODEL FITTING AND FAMILY SELECTION
# =============================================================================

new_maternal_data <- read.csv("new_maternal_data_with_PCA_and_DavidScore.csv")


# Start with Poisson regression (most restrictive assumption)
model_poisson <- glm(offspring_count ~ Maternality_Score + DavidScore +offset(log_age), 
                     data = new_maternal_data, 
                     family = poisson(link = "log"))

summary(model_poisson)

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
# NEGATIVE BINOMIAL MODEL (for overdispersion)
# =============================================================================

# Negative binomial model (allows for overdispersion)
model_nb <- glm.nb(offspring_count ~ Maternality_Score + offset(log_age), 
                   data = new_maternal_data)

summary(model_nb)

# Basic residual plots for negative binomial
par(mfrow = c(2, 2))
plot(model_nb)
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
coef_maternality <- coef(final_model)["Maternality_Score"]
effect_size <- exp(coef_maternality)

cat("\n=== BIOLOGICAL INTERPRETATION ===\n")
cat("For each 1-unit increase in Maternality Score:\n")
cat("Offspring count multiplied by", round(effect_size, 3), "\n")
cat("This represents a", round((effect_size - 1) * 100, 1), "% change in offspring count\n")

# Significance test
p_value <- summary(final_model)$coefficients["Maternality_Score", "Pr(>|z|)"]
cat("P-value for Maternality Score effect:", round(p_value, 4), "\n")

if(p_value < 0.05) {
  cat("*** Maternality Score has a SIGNIFICANT effect on offspring count ***\n")
} else {
  cat("Maternality Score does NOT have a significant effect on offspring count\n")
}

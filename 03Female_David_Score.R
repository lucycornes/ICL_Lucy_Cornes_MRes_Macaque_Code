# Female David's Score Analysis - Simple Clean Version
# ====================================================

# Load required libraries
library(readr)
library(dplyr)
library(glmmTMB)

# 1. LOAD DATA
# ============
focal_data_scaled <- read_csv("focal_data_scaled.csv")

# 2. DEFINE AGGRESSIVE BEHAVIORS 
# ===============================
aggressive_codes <- c(
  "High aggression (give)", "High aggression (receive)", 
  "Low aggression (give)", "Low aggression (receive)",
  "N+LA (G)", "N+HA (g)", "C+HA(G)", "C+LA(G)",
  "N+LA(R)", "C+HA(R)", "C+LA(R)"
)

# 3. PREPARE AGGRESSION DATA
# ==========================
aggression_df <- focal_data_scaled %>%
  filter(Behavior %in% aggressive_codes) %>%
  mutate(
    Winner = ifelse(grepl("give|\\(G\\)|\\(g\\)", Behavior), Actor, Receiver),
    Loser = ifelse(grepl("give|\\(G\\)|\\(g\\)", Behavior), Receiver, Actor)
  )

cat("Total number of aggressive interactions analyzed:", nrow(aggression_df), "\n")
# 4. FUNCTION TO COMPUTE DAVID'S SCORE BY GROUP 
# =============================================
compute_davids_score <- function(group_data, group_label) {
  all_ids <- unique(c(group_data$Winner, group_data$Loser))
  win_matrix <- matrix(0, nrow = length(all_ids), ncol = length(all_ids),
                       dimnames = list(all_ids, all_ids))
  
  for (i in 1:nrow(group_data)) {
    w <- group_data$Winner[i]
    l <- group_data$Loser[i]
    win_matrix[w, l] <- win_matrix[w, l] + 1
  }
  
  total_interactions <- win_matrix + t(win_matrix)
  Pij <- ifelse(total_interactions == 0, 0, win_matrix / total_interactions)
  
  W <- rowSums(Pij)
  L <- colSums(Pij)
  
  W2 <- Pij %*% W
  L2 <- t(Pij) %*% L
  
  DS <- W + W2 - L - L2
  names(DS) <- rownames(win_matrix)
  
  data.frame(
    Individual = names(DS),
    DavidScore = round(as.numeric(DS), 3),
    Group = group_label
  )
}

# 5. CALCULATE DAVID'S SCORES BY GROUP
# ====================================
groups <- unique(aggression_df$Group)
ds_list <- list()

for (grp in groups) {
  group_data <- aggression_df %>% filter(Group == grp)
  
  if (nrow(group_data) == 0) {
    next
  }
  
  ids <- unique(c(group_data$Winner, group_data$Loser))
  if (length(ids) < 2) {
    next
  }
  
  ds_list[[grp]] <- compute_davids_score(group_data, grp)
}

# 6. COMBINE RESULTS AND FILTER TO FOCAL FEMALES
# ===============================================
focal_females <- unique(aggression_df$Actor)
all_aggression_participants <- unique(c(aggression_df$Winner, aggression_df$Loser))

cat("Total focal females (actors):", length(focal_females), "\n")
cat("Focal females who participated in aggression:", 
    length(intersect(focal_females, all_aggression_participants)), "\n")

# Check which focal females (if any) didn't participate in aggression
non_participants <- setdiff(focal_females, all_aggression_participants)
if(length(non_participants) > 0) {
  cat("Focal females who did NOT participate in aggression:", 
      paste(non_participants, collapse = ", "), "\n")
} else {
  cat("All focal females participated in aggressive interactions\n")
}

DS_df <- bind_rows(ds_list) %>%
  filter(Individual %in% focal_females) %>%
  group_by(Group) %>%
  mutate(Rank = rank(-DavidScore, ties.method = "average")) %>%
  ungroup()

cat("Final dataset contains:", nrow(DS_df), "individuals with David's scores\n")

# 7. ADD DEMOGRAPHIC DATA
# =======================
davids_score_final <- DS_df %>%
  rename(AnimalID = Individual) %>%
  left_join(pedigree_data[, c("AnimalID", "offspring_count", "Age")], by = "AnimalID")

# 8. SAVE DATASET
# ===============
write.csv(davids_score_final, "Female_DS_byGroup.csv", row.names = FALSE)


# Print final summary
cat("Final dataset contains", nrow(davids_score_final), "individuals with David's scores\n")
cat("These are the individuals who participated in aggressive interactions\n")



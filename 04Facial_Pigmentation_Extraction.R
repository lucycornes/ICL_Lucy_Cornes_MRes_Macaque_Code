# ==============================================================================
# FACIAL PIGMENTATION EXTRACTION AND ANALYSIS
# Rhesus Macaque Maternal Effects on Same-Sex Sexual Behavior Study
# ==============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(lubridate)

# Try to load image packages (install if needed)
if (!require(jpeg, quietly = TRUE)) {
  install.packages("jpeg")
  library(jpeg)
}

if (!require(png, quietly = TRUE)) {
  install.packages("png")
  library(png)
}

print("=== FACIAL PIGMENTATION ANALYSIS ===")

# ==============================================================================
# 1. CORE REDNESS EXTRACTION FUNCTION
# ==============================================================================

extract_face_redness <- function(image_path, face_region = NULL) {
  
  cat("Processing:", basename(image_path), "\n")
  
  # Determine file type and load accordingly
  file_ext <- tolower(tools::file_ext(image_path))
  
  if (file_ext %in% c("jpg", "jpeg")) {
    img_array <- readJPEG(image_path)
  } else if (file_ext == "png") {
    img_array <- readPNG(image_path)
  } else {
    stop("Unsupported file format. Use JPEG or PNG files.")
  }
  
  # Check if image is grayscale or color
  if (length(dim(img_array)) == 2) {
    stop("Image appears to be grayscale. Need color images for redness analysis.")
  }
  
  # img_array is now height x width x channels (R, G, B)
  height <- dim(img_array)[1]
  width <- dim(img_array)[2]
  
  # Crop to face region if specified (central 60% of image by default)
  if (!is.null(face_region)) {
    y_min <- round(face_region$y_min * height)
    y_max <- round(face_region$y_max * height)
    x_min <- round(face_region$x_min * width)
    x_max <- round(face_region$x_max * width)
    
    # Ensure bounds are valid
    y_min <- max(1, y_min)
    y_max <- min(height, y_max)
    x_min <- max(1, x_min)
    x_max <- min(width, x_max)
    
    img_array <- img_array[y_min:y_max, x_min:x_max, ]
  }
  
  # Extract RGB channels
  r_channel <- img_array[, , 1]
  g_channel <- img_array[, , 2]
  b_channel <- img_array[, , 3]
  
  # Calculate color statistics
  mean_r <- mean(r_channel, na.rm = TRUE)
  mean_g <- mean(g_channel, na.rm = TRUE)
  mean_b <- mean(b_channel, na.rm = TRUE)
  
  # Calculate redness metrics
  redness_metrics <- list(
    # File info
    image_file = basename(image_path),
    
    # Raw RGB values
    raw_red = mean_r,
    raw_green = mean_g,
    raw_blue = mean_b,
    
    # Main redness measures (lighting-robust)
    red_ratio = mean_r / (mean_r + mean_g + mean_b),
    red_dominance = mean_r / mean_g,
    red_minus_green = mean_r - mean_g,
    red_chroma = mean_r / sqrt(mean_r^2 + mean_g^2 + mean_b^2),
    
    # Lighting indicators
    overall_brightness = (mean_r + mean_g + mean_b) / 3,
    contrast = sd(c(r_channel, g_channel, b_channel), na.rm = TRUE),
    overexposed_pct = mean(r_channel > 0.95 | g_channel > 0.95 | b_channel > 0.95),
    underexposed_pct = mean(r_channel < 0.05 & g_channel < 0.05 & b_channel < 0.05),
    
    # Color quality measures
    red_variance = var(as.vector(r_channel), na.rm = TRUE),
    saturation = max(mean_r, mean_g, mean_b) - min(mean_r, mean_g, mean_b)
  )
  
  # Print diagnostics
  cat("  Brightness:", round(redness_metrics$overall_brightness, 3),
      " | Red ratio:", round(redness_metrics$red_ratio, 3),
      " | Overexposed:", round(redness_metrics$overexposed_pct * 100, 1), "%\n")
  
  return(redness_metrics)
}

# ==============================================================================
# 2. BATCH PROCESSING FUNCTION
# ==============================================================================

process_photo_directory <- function(photo_directory, face_crop = TRUE) {
  
  print("=== PROCESSING PHOTO DIRECTORY ===")
  
  # Check if directory exists
  if(!dir.exists(photo_directory)) {
    stop("Directory does not exist: ", photo_directory)
  }
  
  # Find image files
  jpeg_files <- list.files(photo_directory, pattern = "\\.(jpg|jpeg)$", 
                           ignore.case = TRUE, full.names = TRUE)
  png_files <- list.files(photo_directory, pattern = "\\.png$", 
                          ignore.case = TRUE, full.names = TRUE)
  image_files <- c(jpeg_files, png_files)
  
  cat("Found", length(image_files), "image files in", photo_directory, "\n")
  cat("JPEG files:", length(jpeg_files), "| PNG files:", length(png_files), "\n\n")
  
  if(length(image_files) == 0) {
    stop("No JPEG or PNG files found! Check directory path.")
  }
  
  # Define face region (crop to central 60% of image)
  if(face_crop) {
    face_region <- list(x_min = 0.2, x_max = 0.8, y_min = 0.2, y_max = 0.8)
  } else {
    face_region <- NULL
  }
  
  # Process each photo
  redness_data_list <- list()
  successful_extractions <- 0
  
  for(i in 1:length(image_files)) {
    cat("Image", i, "of", length(image_files), ": ")
    
    tryCatch({
      metrics <- extract_face_redness(image_files[i], face_region)
      redness_data_list[[i]] <- metrics
      successful_extractions <- successful_extractions + 1
      
    }, error = function(e) {
      cat("ERROR -", e$message, "\n")
      redness_data_list[[i]] <- NULL
    })
  }
  
  # Convert to dataframe
  valid_data <- redness_data_list[!sapply(redness_data_list, is.null)]
  
  if(length(valid_data) == 0) {
    stop("No photos were successfully processed!")
  }
  
  redness_df <- do.call(rbind, lapply(valid_data, data.frame, stringsAsFactors = FALSE))
  
  cat("\n=== EXTRACTION SUMMARY ===\n")
  cat("Successfully processed:", successful_extractions, "/", length(image_files), "photos\n")
  cat("Final dataframe dimensions:", nrow(redness_df), "rows x", ncol(redness_df), "columns\n")
  
  return(redness_df)
}

# ==============================================================================
# 3. LIGHTING CORRECTION
# ==============================================================================

create_lighting_corrected_redness <- function(redness_data) {
  
  print("=== APPLYING LIGHTING CORRECTION ===")
  
  # Check required columns
  required_cols <- c("red_ratio", "red_dominance", "overall_brightness", "contrast")
  missing_cols <- required_cols[!required_cols %in% names(redness_data)]
  
  if(length(missing_cols) > 0) {
    cat("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    return(redness_data)
  }
  
  # Check data sufficiency
  if(nrow(redness_data) < 3) {
    cat("WARNING: Too few photos (", nrow(redness_data), ") for lighting correction\n")
    return(redness_data)
  }
  
  tryCatch({
    # Statistical residual correction (removes lighting effects)
    cat("Creating residual corrections...\n")
    redness_data$red_ratio_residual <- residuals(lm(red_ratio ~ overall_brightness + contrast, 
                                                    data = redness_data))
    
    redness_data$red_dominance_residual <- residuals(lm(red_dominance ~ overall_brightness + contrast, 
                                                        data = redness_data))
    
    # Brightness-bin standardization
    cat("Creating z-score corrections...\n")
    redness_data <- redness_data %>%
      mutate(
        brightness_category = cut(overall_brightness, 
                                  breaks = 3, 
                                  labels = c("Dark", "Medium", "Bright"))
      ) %>%
      group_by(brightness_category) %>%
      mutate(
        red_ratio_zscore = as.numeric(scale(red_ratio)),
        red_dominance_zscore = as.numeric(scale(red_dominance))
      ) %>%
      ungroup()
    
    # Effectiveness assessment
    original_cor <- cor(redness_data$red_ratio, redness_data$overall_brightness, use = "complete.obs")
    residual_cor <- cor(redness_data$red_ratio_residual, redness_data$overall_brightness, use = "complete.obs")
    
    cat("Lighting correction effectiveness:\n")
    cat("  Original correlation (red ratio ~ brightness):", round(original_cor, 3), "\n")
    cat("  Corrected correlation (residual ~ brightness):", round(residual_cor, 3), "\n")
    
  }, error = function(e) {
    cat("ERROR in lighting correction:", e$message, "\n")
    cat("Continuing without lighting correction...\n")
  })
  
  return(redness_data)
}

# ==============================================================================
# 4. QUALITY ASSESSMENT
# ==============================================================================

assess_photo_quality <- function(redness_data, create_plots = TRUE) {
  
  print("=== PHOTO QUALITY ASSESSMENT ===")
  
  # Basic statistics
  cat("Brightness range:", round(min(redness_data$overall_brightness), 3), 
      "to", round(max(redness_data$overall_brightness), 3), "\n")
  
  brightness_cv <- sd(redness_data$overall_brightness) / mean(redness_data$overall_brightness)
  cat("Brightness variation (CV):", round(brightness_cv, 3), "\n")
  
  # Flag problematic photos
  n_overexposed <- sum(redness_data$overexposed_pct > 0.1)
  n_too_bright <- sum(redness_data$overall_brightness > 0.8)
  n_too_dark <- sum(redness_data$overall_brightness < 0.2)
  
  cat("Potentially problematic photos:\n")
  cat("  Too bright (>0.8):", n_too_bright, "\n")
  cat("  Too dark (<0.2):", n_too_dark, "\n")
  cat("  Overexposed (>10%):", n_overexposed, "\n")
  
  # Create diagnostic plots
  if(create_plots && nrow(redness_data) > 1) {
    # Brightness distribution
    p1 <- ggplot(redness_data, aes(x = overall_brightness)) +
      geom_histogram(bins = 15, fill = "lightblue", alpha = 0.7, color = "black") +
      geom_vline(xintercept = c(0.2, 0.8), linetype = "dashed", color = "red") +
      labs(title = "Photo Brightness Distribution", 
           x = "Overall Brightness", y = "Count") +
      theme_classic()
    
    # Red ratio vs brightness relationship
    p2 <- ggplot(redness_data, aes(x = overall_brightness, y = red_ratio)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(title = "Red Ratio vs Brightness", 
           x = "Overall Brightness", y = "Red Ratio") +
      theme_classic()
    
    print(p1)
    print(p2)
  }
  
  return(redness_data)
}

# ==============================================================================
# 5. MAIN EXTRACTION WORKFLOW
# ==============================================================================

extract_facial_pigmentation <- function(photo_directory, output_file = "facial_pigmentation_data.csv") {
  
  print("=== FACIAL PIGMENTATION EXTRACTION WORKFLOW ===")
  
  # Step 1: Extract redness metrics from photos
  print("Step 1: Extracting redness metrics...")
  redness_data <- process_photo_directory(photo_directory, face_crop = TRUE)
  
  # Step 2: Apply lighting correction
  print("Step 2: Applying lighting correction...")
  redness_data <- create_lighting_corrected_redness(redness_data)
  
  # Step 3: Assess photo quality
  print("Step 3: Assessing photo quality...")
  redness_data <- assess_photo_quality(redness_data, create_plots = TRUE)
  
  # Step 4: Prepare data for export
  print("Step 4: Preparing final dataset...")
  
  # Extract individual IDs from filenames
  redness_data$individual_id <- gsub("\\.(jpg|jpeg|png)$", "", redness_data$image_file, ignore.case = TRUE)
  
  # Select key columns for analysis
  key_columns <- c("individual_id", "image_file", "red_ratio", "red_dominance", 
                   "red_minus_green", "red_ratio_residual", "red_dominance_residual", 
                   "overall_brightness", "contrast", "overexposed_pct",
                   "raw_red", "raw_green", "raw_blue")
  
  export_data <- redness_data[, key_columns[key_columns %in% names(redness_data)]]
  
  # Add photo date column (placeholder for manual entry)
  export_data$photo_date <- "YYYY-MM-DD"
  
  # Save results
  write.csv(export_data, output_file, row.names = FALSE)
  
  print("=== EXTRACTION WORKFLOW COMPLETE ===")
  cat("Results saved to:", output_file, "\n")
  cat("Dataset contains:", nrow(export_data), "photos from", 
      length(unique(export_data$individual_id)), "individuals\n")
  
  # Display summary of individuals
  cat("Individual IDs found:\n")
  individual_summary <- table(export_data$individual_id)
  print(individual_summary)
  
  cat("\n⚠️  IMPORTANT: Add real photo dates to 'photo_date' column before analysis!\n")
  cat("   Replace 'YYYY-MM-DD' with actual dates in format: YYYY-MM-DD\n")
  
  return(redness_data)
}

# ==============================================================================
# 6. INTEGRATE WITH BEHAVIORAL DATA
# ==============================================================================

analyze_pigmentation_behavior_association <- function(pigmentation_file = "facial_pigmentation_data.csv",
                                                      behavioral_file = "female_behaviour.csv",
                                                      offspring_file = "Offspring Count_anonymised.csv") {
  
  print("=== PIGMENTATION-BEHAVIOR ASSOCIATION ANALYSIS ===")
  
  # Load pigmentation data
  print("Loading pigmentation data...")
  pigmentation_data <- read.csv(pigmentation_file, stringsAsFactors = FALSE)
  
  # Check for date completion
  placeholder_dates <- sum(pigmentation_data$photo_date %in% c("YYYY-MM-DD", "", NA), na.rm = TRUE)
  if(placeholder_dates > 0) {
    cat("WARNING:", placeholder_dates, "photos still have placeholder dates\n")
    cat("Add real dates before proceeding with temporal analysis\n")
  }
  
  # Load behavioral data
  print("Loading behavioral data...")
  behavioral_data <- read.csv(behavioral_file, stringsAsFactors = FALSE)
  
  # Load offspring data
  print("Loading offspring data...")
  offspring_data <- read.csv(offspring_file, stringsAsFactors = FALSE)
  
  # Merge datasets
  print("Merging datasets...")
  
  # First merge pigmentation with behavioral data
  merged_data <- merge(pigmentation_data, behavioral_data, 
                       by.x = "individual_id", by.y = "Actor", all = FALSE)
  
  # Then merge with offspring data
  merged_data <- merge(merged_data, offspring_data %>% select(AnimalID, offspring_count, BirthSeason), 
                       by.x = "individual_id", by.y = "AnimalID", all.x = TRUE)
  
  cat("Final merged dataset:", nrow(merged_data), "individuals\n")
  
  if(nrow(merged_data) == 0) {
    cat("ERROR: No individuals found in all datasets!\n")
    return(NULL)
  }
  
  # Prepare analysis variables
  merged_data <- merged_data %>%
    mutate(
      Age = 2025 - BirthSeason,
      log_age = log(Age)
    )
  
  # Statistical analysis - test association with offspring count
  print("=== STATISTICAL ANALYSIS ===")
  
  if(!is.null(merged_data$offspring_count) && sum(!is.na(merged_data$offspring_count)) >= 5) {
    
    cat("Testing association between facial pigmentation and reproductive success...\n")
    
    # Fit models for different redness measures
    models <- list()
    
    # Model 1: Red ratio (lighting corrected)
    if("red_ratio_residual" %in% names(merged_data)) {
      models$red_ratio <- glm(offspring_count ~ red_ratio_residual + offset(log(Age)), 
                              family = poisson(), data = merged_data)
    }
    
    # Model 2: Red dominance (lighting corrected) 
    if("red_dominance_residual" %in% names(merged_data)) {
      models$red_dominance <- glm(offspring_count ~ red_dominance_residual + offset(log(Age)), 
                                  family = poisson(), data = merged_data)
    }
    
    # Print results
    for(model_name in names(models)) {
      cat("\n--- Model:", model_name, "---\n")
      model_summary <- summary(models[[model_name]])
      print(model_summary)
      
      # Extract effect size
      coef_table <- model_summary$coefficients
      predictor_name <- paste0(model_name, "_residual")
      if(predictor_name %in% rownames(coef_table)) {
        effect <- coef_table[predictor_name, "Estimate"]
        se <- coef_table[predictor_name, "Std. Error"]
        p_val <- coef_table[predictor_name, "Pr(>|z|)"]
        
        cat("Effect size:", round(effect, 4), "±", round(se, 4), "\n")
        cat("P-value:", round(p_val, 4), "\n")
        cat("95% CI: [", round(effect - 1.96*se, 4), ",", round(effect + 1.96*se, 4), "]\n")
      }
    }
    
    # Create visualization
    if("red_ratio_residual" %in% names(merged_data)) {
      plot_data <- merged_data[!is.na(merged_data$red_ratio_residual) & 
                                 !is.na(merged_data$offspring_count), ]
      
      p <- ggplot(plot_data, aes(x = red_ratio_residual, y = offspring_count)) +
        geom_point(alpha = 0.7, size = 2.5) +
        geom_smooth(method = "glm", method.args = list(family = "poisson"), 
                    se = TRUE, color = "red") +
        labs(x = "Facial Redness (Lighting-Corrected)",
             y = "Offspring Count",
             title = "Facial Pigmentation and Reproductive Success",
             subtitle = paste("n =", nrow(plot_data), "females")) +
        theme_classic()
      
      print(p)
      
      # Save plot
      ggsave("facial_pigmentation_offspring_plot.png", p, width = 8, height = 6, dpi = 300)
    }
    
  } else {
    cat("Insufficient offspring data for statistical analysis\n")
  }
  
  # Save merged dataset
  write.csv(merged_data, "pigmentation_behavioral_merged.csv", row.names = FALSE)
  
  print("=== ANALYSIS COMPLETE ===")
  print("Files created:")
  print("- pigmentation_behavioral_merged.csv (merged dataset)")
  if(exists("p")) {
    print("- facial_pigmentation_offspring_plot.png (visualization)")
  }
  
  return(list(
    data = merged_data,
    models = if(exists("models")) models else NULL,
    sample_size = nrow(merged_data)
  ))
}

# ==============================================================================
# 7. EXAMPLE USAGE
# ==============================================================================

# Set your photo directory path
photo_directory <- "/Users/lucycornes/Desktop/New_female_faces"

# Run the complete extraction workflow
print("Starting facial pigmentation extraction...")
print("Note: Update the photo_directory path above to match your system")

# This would extract pigmentation data from photos
# pigmentation_results <- extract_facial_pigmentation(photo_directory, "facial_pigmentation_data.csv")

# After manually adding dates to the CSV file, run the behavioral analysis
# analysis_results <- analyze_pigmentation_behavior_association()

print("=== FACIAL PIGMENTATION ANALYSIS SCRIPT READY ===")
print("To run:")
print("1. Update photo_directory path")
print("2. Run: pigmentation_results <- extract_facial_pigmentation(photo_directory)")
print("3. Manually add photo dates to the output CSV")
print("4. Run: analysis_results <- analyze_pigmentation_behavior_association()")
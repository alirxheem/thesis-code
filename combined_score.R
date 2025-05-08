# run after substrate.list.generation
# create the combined score list
combined_kinase_scores <- list()

# Iterate over each kinase in kinase_scores
for (kinase in names(kinase_scores)) {
  
  # Extract the Score array for the kinase from kinase_scores
  score_array <- kinase_scores[[kinase]]$input_control$Score
  
  # Ensure the kinase exists in the profileScoreMatrix column names
  if (kinase %in% colnames(profileScoreMatrix)) {
    
    # Access the kinase column properly using dataframe syntax
    profile_scores <- profileScoreMatrix[, kinase]  # <-- CORRECTED ACCESS
    
    # Ensure the dimensions match
    if (length(score_array) == nrow(profileScoreMatrix)) {
      
      # Add the two score arrays element-wise
      combined_score <- score_array + profile_scores
      
      # Store the combined scores in the list
      combined_kinase_scores[[kinase]] <- combined_score
      
    } else {
      cat(paste("Warning: Dimension mismatch for kinase:", kinase, "\n"))
    }
  } else {
    cat(paste("Warning: Kinase", kinase, "not found in profileScoreMatrix\n"))
  }
}

print(combined_kinase_scores)

# Create a matrix with 20583 rows (one for each row in our_data) 
# and columns for each kinase in combined_kinase_scores
kinase_names_tf <- names(combined_kinase_scores)
num_rows_tf <- nrow(our_data)

# create matrix with row names as substrates and columns as kinases
result_matrix <- matrix(NA, nrow = num_rows_tf, ncol = length(kinase_names_tf),
                        dimnames = list(our_data$New_Gene_Residue, kinase_names_tf))

# Loop through each kinase and populate the matrix with scores
for (kinase in kinase_names_tf) {
  kinase_scores_tf <- combined_kinase_scores[[kinase]]  # Get the scores for the kinase
  result_matrix[, kinase] <- kinase_scores_tf  # Place scores in the matrix by row order
}

# Convert the matrix to a data frame for easier handling
result_df <- as.data.frame(result_matrix)


head(result_df)




###combined score heat map


library(ggplot2)
library(pheatmap)

#  Preprocessing: Handle NAs and Inf values
# Replace NAs with 0 for heatmap visualization
result_df[is.na(result_df)] <- 0
result_df[!is.finite(as.matrix(result_df))] <- 0

pheatmap(
  as.matrix(result_df),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("green", "white", "purple"))(100),
  main = "Kinase-Phosphosite Heatmap",
  fontsize_row = 8,
  fontsize_col = 10,
  na_col = "white",
  filename = file.path("~/Documents/THESIS/Plots/transcription_factors", paste0(kinase, "combinedscore_heatmap.png"))
)







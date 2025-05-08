#run after combined_score 




library(readxl)


#tf_table <- read_excel("~/Documents/THESIS/directories/Mouse_TFs_Kinases_webpage-3-30-2017.xlsx", sheet = 1)

# Extract transcription factors (capitalized)
#transcription_factors <- toupper(tf_table$'Gene Symbol')

# create an empty list to store kinase-TF predictions
#kinase_tf_list <- list()

# Loop through the substrate list
#for (kinase in names(substrate.list)) {
  
  # Extract the substrates for the current kinase
 # substrates_pro <- substrate.list[[kinase]]
  
  # Extract the gene names (everything before the first semicolon)
 # substrate_genes <- sapply(substrates_pro, function(x) strsplit(x, ";")[[1]][1])
  
  # Identify transcription factors in the substrates
 # tf_matches <- intersect(substrate_genes, transcription_factors)
  
  # Add to the list only if there are matches
 # if (length(tf_matches) > 0) {
 #   kinase_tf_list[[kinase]] <- tf_matches
 # }
#}


#print(kinase_tf_list)




######sanity
# Sanity check: verify kinase-TF pairs exist in substrate.list

# create lists to store matched and mismatched pairs
#matched_pairs <- list()
#mismatched_pairs <- list()

# Loop through each kinase in kinase_tf_list
#for (kinase in names(kinase_tf_list)) {
  # TFs predicted by the kinase
#  predicted_tfs <- kinase_tf_list[[kinase]]
  
  # Extract the substrates from substrate.list
#  kinase_substrates <- substrate.list[[kinase]]
  
  # Extract the gene names from the substrate phosphosites
#  substrate_genes <- unique(sapply(kinase_substrates, function(x) strsplit(x, ";")[[1]][1]))
  
  # Find matches and mismatches
#  matched <- intersect(predicted_tfs, substrate_genes)
#  mismatched <- setdiff(predicted_tfs, substrate_genes)
  
  # Store the results
#  matched_pairs[[kinase]] <- matched
#  mismatched_pairs[[kinase]] <- mismatched
#}

# Display the results
#cat("\n Matched Pairs (Present in both lists):\n")
#print(matched_pairs)

#cat("\n Mismatched Pairs (Only in kinase_tf_list):\n")
#print(mismatched_pairs)







########phosphosite tf

tf_table_phos <- read_excel("~/Documents/THESIS/directories/Mouse_TFs_Kinases_webpage-3-30-2017.xlsx", sheet = 1)

# Extract transcription factors (capitalized)
transcription_factors <- toupper(tf_table_phos$'Gene Symbol')

# create an empty list to store kinase-TF phosphosite predictions
kinase_tf_phosphosite_list <- list()

# Loop through the substrate list
for (kinase in names(substrate.list)) {
  
  # Extract the full phosphosites for the current kinase
  substrates_phos <- substrate.list[[kinase]]
  
  # Extract the gene names (everything before the first semicolon)
  substrate_genes_phos <- sapply(substrates_phos, function(x) strsplit(x, ";")[[1]][1])
  
  # Identify matching TFs
  tf_matches_phos <- intersect(substrate_genes_phos, transcription_factors)
  
  # If there are TF matches, include all associated phosphosites
  if (length(tf_matches_phos > 0)) {
    # Filter all phosphosites for the matching TFs
    phosphosites_phos <- substrates_phos[substrate_genes_phos %in% tf_matches_phos]
    
    # Store the kinase with its TF phosphosites
    kinase_tf_phosphosite_list[[kinase]] <- phosphosites_phos
  }
}

# Remove trailing semicolon from each substrate
kinase_tf_phosphosite_list <- lapply(kinase_tf_phosphosite_list, function(substrates) {
  gsub(";$", "", substrates)
})


print(kinase_tf_phosphosite_list)



library(ggplot2)
library(pheatmap)

# create an empty list to store the combined scores for each kinase-substrate pair
master_table <- data.frame(Kinase = character(),
                           Substrate = character(),
                           Combined_Score = numeric(),
                           stringsAsFactors = FALSE)

#  Create a heatmap for each kinase with its respective substrates
output_dir <- "~/Documents/THESIS/Plots/transcription_factor"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Iterate over each kinase in kinase_tf_phosphosite_list
for (kinase in names(kinase_tf_phosphosite_list)) {
  
  # Extract the substrates for the current kinase
  substrates_tf <- kinase_tf_phosphosite_list[[kinase]]
  
  # Filter result_df to include only the substrates for the current kinase
  valid_substrates <- intersect(substrates_tf, rownames(result_df))
  
  # Get the number of substrates
  num_substrates <- length(valid_substrates)
  
  if (num_substrates > 0) {
    # Create a submatrix with only the kinase's substrates
    kinase_matrix <- result_df[valid_substrates, kinase, drop = FALSE]
    
    # Store the kinase-substrate pairs and their scores in the master table
    temp_table <- data.frame(
      Kinase = rep(kinase, num_substrates),
      Substrate = valid_substrates,
      Combined_Score = as.vector(kinase_matrix[, 1]),
      stringsAsFactors = FALSE
    )
    
    # add the temporary table to the master table
    master_table <- rbind(master_table, temp_table)
    
    attr(kinase_matrix, "label") <- "Kinase Score"
    library(grid)

    pheatmap(
      as.matrix(kinase_matrix),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = colorRampPalette(c("green", "white", "purple"))(100),
      main = paste(kinase, "- Number of Substrates (Transcription Factors):", num_substrates),
      fontsize_row = 8,
      fontsize_col = 10,
      na_col = "white",
      filename = file.path(output_dir, paste0(kinase, "_heatmap.png"))
    )

    cat("Generated heatmap for:", kinase, "with", num_substrates, "substrates\n")
  } else {
    cat("No valid substrates for:", kinase, "\n")
  }
}

# Save the master table as a CSV file
write.csv(master_table, file = "~/Documents/THESIS/Plots/kinase_substrate_scores.csv", row.names = FALSE)

cat("\nAll heatmaps saved in '~/Documents/THESIS/Plots/transcription_factors' folder!\n")
cat("\nMaster table saved as 'kinase_substrate_scores.csv'!\n")
















#library(ggplot2)
#library(reshape2)

#Get the list of unique phosphosites from kinase_tf_phosphosite_list
#all_substrates <- unique(unlist(lapply(kinase_tf_phosphosite_list, function(x) gsub("'", "", x))))

# Filter result_df to include only substrates in kinase_tf_phosphosite_list
#filtered_df <- result_df[rownames(result_df) %in% all_substrates, ]

# Create a dataframe for heatmap by melting the filtered dataframe
#heatmap_data <- melt(as.matrix(filtered_df))

# Rename the columns for readability
#colnames(heatmap_data) <- c("Phosphosite", "Kinase", "Score")

#Mask the non-interacting kinase-phosphosite pairs
# Create a new dataframe to mark kinase-phosphosite pairs
#heatmap_data$Masked_Score <- NA

# Iterate through each kinase in kinase_tf_phosphosite_list and mask irrelevant pairs
#for (kinase in names(kinase_tf_phosphosite_list)) {
#  substrates <- gsub("'", "", kinase_tf_phosphosite_list[[kinase]])  # Clean substrate names
#  heatmap_data$Masked_Score[heatmap_data$Kinase == kinase & 
#                              heatmap_data$Phosphosite %in% substrates] <- heatmap_data$Score
#}


#ggplot(heatmap_data, aes(x = Kinase, y = Phosphosite, fill = Masked_Score)) +
 # geom_tile() +
#  scale_fill_gradient2(low = "green", mid = "white", high = "purple", midpoint = median(heatmap_data$Score, na.rm = TRUE)) +
 # theme_minimal() +
 # labs(title = "Kinase-Phosphosite Interaction Heatmap",
#       x = "Kinase",
#       y = "Phosphosite",
#       fill = "Combined Score") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1),
#        axis.text.y = element_text(size = 8))






# run after predictions made
# Add a new column 'row_number' to our_data that represents the row index
our_data$row_number <- 1:nrow(our_data)

# create an empty list to store the substrates for each kinase
substrate.list <- list()

# Loop through each kinase in finalPredictedSub
for (kinase in names(finalPredictedSub)) {
  # Get the row indices for the current kinase from finalPredictedSub
  substrate_rows <- finalPredictedSub[[kinase]]
  
  # create an empty vector to store substrates for this kinase
  substrates <- character(0)
  
  # Loop through each row index in substrate_rows
  for (row_num in substrate_rows) {
    # Ensure that the row_num exists in the 'row_number' column
    matched_row <- our_data[our_data$row_number == row_num, ]
    
    # If we found a matching row, extract the substrate (New_Gene_Residue)
    if (nrow(matched_row) > 0) {
      substrates <- c(substrates, toupper(paste0(matched_row$New_Gene_Residue,";")))
    } else {
      # If no match is found, append NA (or handle this differently if needed)
      substrates <- c(substrates, NA)
    }
  }
  
  # Store the substrates for the current kinase in the substrate.list
  substrate.list[[kinase]] <- substrates
}

# run after combined_score




# List of target genes to check
genes_of_interest <- c("NANOG", "SOX2", "POU5F1", "MYC;", "STAT3", "LIN28A")

# create an empty list to store the results
results <- list()

# For loop to check for target genes in each kinase's substrates
for (kinase in names(substrate.list)) {
  substrates <- substrate.list[[kinase]]  # Get the substrates for the current kinase
  
  # Check if any of the substrates contain any of the genes of interest
  matching_substrates <- substrates[sapply(substrates, function(sub) {
    any(sapply(genes_of_interest, function(gene) grepl(gene, sub)))
  })]
  
  # Store the matching substrates in the results list
  if (length(matching_substrates) > 0) {
    results[[kinase]] <- matching_substrates
  }
}


results

# Create an empty data frame to store results
kinase_data <- data.frame(
  kinase = character(),
  substrate = character(),
  score = numeric(),
  stringsAsFactors = FALSE
)

# Named list of kinase-to-substrate mappings cell cycle regulated
kinase_map <- list(
  CDK1     = c("STAT3;S727;"),
  CDK2     = c("STAT3;S727;"),
  CDK5     = c("MYC;S62;", "STAT3;S727;"),
  CSNK2A1  = c("POU5F1;S102;"),
  GSK3B    = c("MYC;S62;", "STAT3;S727;"),
  MAPK1    = c("LIN28A;S200;"),
  MAPK14   = c("MYC;S62;", "STAT3;S727;"),
  MAPK3    = c("LIN28A;S200;"),
  MAPK9    = c("LIN28A;S200;"),
  MTOR     = c("MYC;S62;", "STAT3;S727;"),
  PRKAA1   = c("POU5F1;S102;"),
  RPS6KB1  = c("POU5F1;S102;")
)

# Named list of kinase-to-substrate mappings non-cell cycle regulated
#kinase_map <- list(
#  MAPK1    = c("STAT3;S727;"),
#  MAPK14   = c("STAT3;S727;"),
#  MAPK3    = c("POU5F1;S12;"),
#  MTOR     = c("SOX2;S253;", "STAT3;S727;")
#)

# Loop over the map and extract scores from result_df
for (kinase in names(kinase_map)) {
  for (substrate in kinase_map[[kinase]]) {
    # Get score from result_df
    score <- result_df[substrate, kinase]
    
    # Add to final data frame
    kinase_data <- rbind(kinase_data, data.frame(
      kinase = kinase,
      substrate = substrate,
      score = score,
      stringsAsFactors = FALSE
    ))
  }
}








library(ggplot2)


ggplot(kinase_data, aes(x = substrate, y = kinase, color = score)) +
  geom_point(size = 4) +  # Adjust point size as necessary
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(kinase_data$score)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(
    title = "Cell-Cycle Regulated Pluripotent-Associated\n Kinase Predictions", 
    x = "Substrate", 
    y = "Kinase", 
    color = "Kinase Score"
  ) +
  theme(legend.title = element_text(size = 12))  # Adjust the legend title size if necessary


#SANTIY: kinase_scores[["RPS6KB1"]]$input_control$Score[3079]


##sanity check
# create an empty data frame to store results
#substrate_table <- data.frame(Kinase = character(), Substrate = character(), stringsAsFactors = FALSE)

# Loop through the substrate list and convert it to a table
#for (kinase in names(substrate.list)) {
# substrates <- substrate.list[[kinase]]
# for (substrate in substrates) {
# Add each substrate as a row in the data frame
#   substrate_table <- rbind(substrate_table, data.frame(Kinase = kinase, Substrate = substrate, stringsAsFactors = FALSE))
#  }
#}


#print(substrate_table)

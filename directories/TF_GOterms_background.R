# Extract phosphosites from kinase_tf_phosphosite_list
all_phosphosites_tf <- unique(unlist(kinase_tf_phosphosite_list))
# Filter our_data for matching phosphosites
tf_datainput <- our_data %>% filter(New_Gene_Residue %in% all_phosphosites_tf) %>%  
  dplyr::select(New_Gene_Residue, starts_with("phase_"))         # Keep phosphosites + phase_ columns

tf_datainput <- as.data.frame(tf_datainput)

head(tf_datainput)



# Directory to save the cluster files
output_dir <- "~/Documents/THESIS/clusters"

# Loop over each kinase in kinase_tf_phosphosite_list
for (i in 1:length(kinase_tf_phosphosite_list)) {
  # Get the current kinase and its substrates (phosphosites)
  kinase <- names(kinase_tf_phosphosite_list)[i]
  substrates <- kinase_tf_phosphosite_list[[i]]
  
  # Filter the substrates that match the New_Gene_Residue in tf_datainput
  matched_phosphosites <- tf_datainput %>%
    filter(New_Gene_Residue %in% substrates)  # Filter based on the substrates
  
  # Create the data frame with relevant columns
  cluster_df <- matched_phosphosites %>%
    dplyr::select(New_Gene_Residue, starts_with("phase_")) %>%
    rename(RowNames = New_Gene_Residue) %>%
    mutate(cluster = i)  # Add the Cluster column (set to the current cluster number)
  
  # Write the data frame to an Excel file for this cluster
  write.xlsx(cluster_df, file.path(output_dir, paste0("cluster_", i, ".xlsx")), rowNames = FALSE)
  
  print(paste("Saved cluster file for kinase", kinase, "as cluster_", i, ".xlsx"))
}




# Load necessary libraries
library(readxl)
library(clusterProfiler)
library(org.Mm.eg.db)  # Use org.Hs.eg.db for human
library(openxlsx)
library(ggplot2)
library(ComplexHeatmap)
library(stringr)

# Define the directory where your cluster files are stored
cluster_dir <- "~/Documents/THESIS/clusters_background"

# Get the list of cluster files
cluster_files <- list.files(cluster_dir, pattern = "cluster_.*\\.xlsx", full.names = TRUE)

# Read each cluster and extract genes
clusters_list <- lapply(cluster_files, function(file) {
  data <- read_excel(file)  # Read each cluster file
  genes <- unique(data$RowNames)  # Extract unique gene names
  return(genes)
})
names(clusters_list) <- paste0("cluster_", seq_along(clusters_list))  # Assign names

# Cleaning function to remove extra characters
clean_gene_names <- function(genes) {
  gsub(";.*", "", genes)  # Remove everything after the first semicolon
}

# Apply cleaning to each cluster
clusters_list <- lapply(clusters_list, clean_gene_names)

# Convert uppercase gene names to CamelCase (if needed)
convert_to_camel_case <- function(gene) {
  str_to_title(tolower(gene))  # Convert to lowercase first, then title case
}

clusters_list <- lapply(clusters_list, function(genes) {
  sapply(genes, convert_to_camel_case)  # Convert each gene name
})

# --- SET BACKGROUND GENES (All transcription factors in kinase_tf_phosphosite_list) ---
background_genes <- unique(unlist(kinase_tf_phosphosite_list))  # Extract all unique TFs
background_genes <- clean_gene_names(background_genes)  # Clean names
background_genes <- sapply(background_genes, convert_to_camel_case)  # Convert names

# Verify background gene count
print(paste("Total background genes:", length(background_genes)))

# --- Perform GO enrichment analysis with BACKGROUND ---
enrichment_results <- lapply(clusters_list, function(gene_set) {
  enrichGO(
    gene = gene_set,
    OrgDb = org.Mm.eg.db,  # Mouse database
    keyType = "SYMBOL",
    ont = "BP",  # "BP" = Biological Process (useful for lineage differentiation)
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    universe = background_genes  # Setting the background here
  )
})

# Convert results to a dataframe for easy viewing
enrichment_df <- lapply(enrichment_results, function(res) {
  if (is.null(res)) return(NULL)  # Skip empty results
  as.data.frame(res)  # Convert enrichment output to data frame
})

# Save enrichment results
write.xlsx(enrichment_df, "Lineage_Specific_Enrichment.xlsx")

# --- PLOTTING ---
output_dir <- "~/Documents/THESIS/clusters_background"

for (i in 1:length(enrichment_results)) {
  kinase_name <- names(kinase_tf_phosphosite_list)[i]
  num_phosphosites <- length(kinase_tf_phosphosite_list[[i]])
  
  dotplot_file <- file.path(output_dir, paste0("cluster_", i, "_", kinase_name, "_dotplot.png"))
  barplot_file <- file.path(output_dir, paste0("cluster_", i, "_", kinase_name, "_barplot.png"))
  
  # Check if enrichment_results[[i]] is not NULL and contains valid data
  if (!is.null(enrichment_results[[i]]) && nrow(as.data.frame(enrichment_results[[i]])) > 0) {
    
    dotplot_obj <- dotplot(enrichment_results[[i]], showCategory = 10, 
                           title = paste("cluster", i, "-", kinase_name, "(", num_phosphosites, "Phosphosites)"))
    ggsave(dotplot_file, plot = dotplot_obj, width = 8, height = 6, dpi = 300)
    
    barplot_obj <- barplot(enrichment_results[[i]], showCategory = 10, 
                           title = paste("Top 10 Lineage-Specific Processes (cluster", i, "-", kinase_name, ")", 
                                         "(", num_phosphosites, "Phosphosites)"))
    ggsave(barplot_file, plot = barplot_obj, width = 8, height = 6, dpi = 300)
    
    print(paste("Saved plots for cluster", i, "-", kinase_name, "with", num_phosphosites, "phosphosites"))
  } else {
    print(paste("No significant enrichment for cluster", i, "-", kinase_name))
  }
}

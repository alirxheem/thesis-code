#run after predictions are made
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

# Load location data 
cell_loc <- read_tsv("~/Documents/THESIS/directories/subcellular_location.tsv", show_col_types = FALSE)

# Define broader categories once
nuclear_structures <- c("Nucleoplasm", "Nuclear speckles", "Nucleoli", 
                        "Nuclear bodies", "Nuclear membrane", "Nucleoli rim", 
                        "Nucleoli fibrillar center")

cytoplasmic_structures <- c("Cytosol", "Plasma membrane", "Vesicles", 
                            "Golgi apparatus", "Cell Junctions", "Mitochondria", 
                            "Centrosome", "Endoplasmic reticulum", "Actin filaments", 
                            "Microtubules", "Centriolar satellite", "Focal adhesion sites", 
                            "Intermediate filaments", "Mitotic spindle", "Mitotic chromosome", 
                            "Microtubule ends", "Cytoplasmic bodies", "Cytokinetic bridge", 
                            "Midbody ring", "Kinetochore", "Lysosomes", "Endosomes", 
                            "Peroxisomes", "Midbody", "Rods & Rings")

process_kinase <- function(kinase_name, data, sub_list, cell_loc) {
  kinase_data <- data[sub_list[[kinase_name]], ] %>%
    mutate(Gene = toupper(Gene)) %>%
    mutate(Main_location = cell_loc$`Main location`[match(Gene, cell_loc$`Gene name`)],
           Additional_location = cell_loc$`Additional location`[match(Gene, cell_loc$`Gene name`)],
           Extracellular_location = cell_loc$`Extracellular location`[match(Gene, cell_loc$`Gene name`)]) %>%
    filter(!is.na(Main_location)) %>%
    separate_rows(Main_location, sep = ";") %>%
    mutate(Broader_Location = case_when(
      Main_location %in% nuclear_structures ~ "Nuclear",
      Main_location %in% cytoplasmic_structures ~ "Cytoplasmic",
      TRUE ~ "Other"
    )) %>%
    group_by(Gene) %>%
    mutate(Location_Type = ifelse(n_distinct(Broader_Location) > 1, "Multiple Locations", "Unique Location")) %>%
    ungroup() %>%
    filter(Location_Type == "Unique Location") %>%
    group_by(Broader_Location) %>%
    summarise(phosphosite_count = n(), .groups = 'drop') %>%
    filter(Broader_Location %in% c("Nuclear", "Cytoplasmic")) %>%  # Only keep the two main types
    arrange(Broader_Location) %>%
    mutate(Total = sum(phosphosite_count),
           Unique_Location_Percentage = (phosphosite_count / Total) * 100) %>%
    summarise(
      Kinase = kinase_name,
      Difference = abs(diff(phosphosite_count)),
      enrichment_difference = diff(Unique_Location_Percentage),
      Total = sum(phosphosite_count),
      Enriched_phosphosite_count = max(phosphosite_count)
    )
  
  return(kinase_data)
}



# Apply function to all kinases
phosphosite_enrichment_all <- map_dfr(names(finalPredictedSub), process_kinase,
                                      data = our_data,
                                      sub_list = finalPredictedSub,
                                      cell_loc = cell_loc)

phosphosite_enrichment_all <- map_dfr(names(finalPredictedSub), process_kinase,
                                      data = our_data,
                                      sub_list = finalPredictedSub,
                                      cell_loc = cell_loc)

phosphosite_enrichment_all <- phosphosite_enrichment_all %>%
  arrange(enrichment_difference) %>%
  mutate(Kinase = factor(Kinase, levels = Kinase))

phosphosite_enrichment_all <- phosphosite_enrichment_all %>%
  mutate(
    Enriched_location = ifelse(enrichment_difference > 0, "Nuclear", "Cytoplasmic"),
    Label = paste0(
      Enriched_phosphosite_count, " (", Enriched_location, ")\n",
      "\u2015\u2015\u2015\u2015\u2015\u2015\n",  # Longer solid dash line (Unicode)
      Total, " (Total)"
    )
  )








ggplot(phosphosite_enrichment_all, aes(x = Kinase, y = enrichment_difference)) +
  geom_col(fill = "purple", color = "purple", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +
  geom_text(
    aes(label = Label, vjust = ifelse(enrichment_difference > 0, -0.25, 1.4)),
    color = "black", size = 1.9
  ) +
  labs(
    title = "Subcellular Enrichment of Kinase Substrates",
    x = "Kinase", y = "Nuclear - Cytoplasmic (% Difference)"
  ) +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    limits = c(-110, 40),
    breaks = seq(-100, 30, by = 25)
  )

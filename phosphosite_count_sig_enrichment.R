library(pheatmap)
library(ggplot2)
library(PhosR)
library(limma)
library(GGally)
library(readxl)
library(openxlsx)
library(dplyr)
library(stringr)
library(tidyr)
library(gplots)
library(calibrate)
library(directPA)
library(org.Rn.eg.db)
library(reactome.db)
library(annotate)
library(SummarizedExperiment)
library(enrichR)
library(tibble)
library(robustHD)
library(NbClust)
library(reshape2)
library(tidyverse)

# Load datasets LG1
CCReg_G1 <- read.table("~/Documents/THESIS/directories/signalome_construction/LateG1_Significant_Collapsed.txt",  header=T, sep="\t")
NCCReg_G1 <- read.table("~/Documents/THESIS/directories/signalome_construction/LateG1_Collapsed_NonCCRegulated.txt",  header=T, sep="\t")
cell_loc <- read_tsv("~/Documents/THESIS/directories/subcellular_location.tsv", show_col_types = FALSE)

# Match subcellular location for CCReg and NCCReg based on gene names
CCReg_G1$Gene <- toupper(CCReg_G1$Gene)
NCCReg_G1$Gene <- toupper(NCCReg_G1$Gene)

# Add location columns by matching gene names with cell_loc data
CCReg_G1 <- CCReg_G1 %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

NCCReg_G1 <- NCCReg_G1 %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

# Remove rows with NA in Main_location and separate multiple locations
CCReg_G1 <- CCReg_G1 %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

NCCReg_G1 <- NCCReg_G1 %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

# Calculate total counts for each Main_location in both datasets
total_counts_G1 <- CCReg_G1 %>%
  group_by(Main_location) %>%
  summarise(CC_count = n()) %>%
  left_join(NCCReg_G1 %>% group_by(Main_location) %>% summarise(NCC_count = n()), by = "Main_location") %>%
  mutate(
    total_count = CC_count + NCC_count,
    CC_percentage = (CC_count / total_count),
    NCC_percentage = (NCC_count / total_count),
    enrichment_difference = (CC_percentage - NCC_percentage) * 100
  )

# Perform Fisher's exact test for each subcellular location
total_counts_G1 <- total_counts_G1 %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(CC_count, NCC_count, total_count - CC_count, total_count - NCC_count), nrow = 2))$p.value,
    significance = ifelse(p_value < 0.05, TRUE, FALSE)
  )

# Filter to only significant subcellular locations
significant_counts_G1 <- total_counts_G1 %>%
  filter(significance == TRUE)

# Load datasets S
CCReg_S <- read.table("~/Documents/THESIS/directories/signalome_construction/S_Significant_Collapsed.txt",  header=T, sep="\t")
NCCReg_S <- read.table("~/Documents/THESIS/directories/signalome_construction/S_Collapsed_NonCCRegulated.txt",  header=T, sep="\t")
cell_loc <- read_tsv("~/Documents/THESIS/directories/subcellular_location.tsv", show_col_types = FALSE)

# Match subcellular location for CCReg and NCCReg based on gene names
CCReg_S$Gene <- toupper(CCReg_S$Gene)
NCCReg_S$Gene <- toupper(NCCReg_S$Gene)

# Add location columns by matching gene names with cell_loc data
CCReg_S <- CCReg_S %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

NCCReg_S <- NCCReg_S %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

# Remove rows with NA in Main_location and separate multiple locations
CCReg_S <- CCReg_S %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

NCCReg_S <- NCCReg_S %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

# Calculate total counts for each Main_location in both datasets
total_counts_S <- CCReg_S %>%
  group_by(Main_location) %>%
  summarise(CC_count = n()) %>%
  left_join(NCCReg_S %>% group_by(Main_location) %>% summarise(NCC_count = n()), by = "Main_location") %>%
  mutate(
    total_count = CC_count + NCC_count,
    CC_percentage = (CC_count / total_count),
    NCC_percentage = (NCC_count / total_count),
    enrichment_difference = (CC_percentage - NCC_percentage) * 100
  )

# Perform Fisher's exact test for each subcellular location
total_counts_S <- total_counts_S %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(CC_count, NCC_count, total_count - CC_count, total_count - NCC_count), nrow = 2))$p.value,
    significance = ifelse(p_value < 0.05, TRUE, FALSE)
  )

# Filter to only significant subcellular locations
significant_counts_S <- total_counts_S %>%
  filter(significance == TRUE)

# Load datasets G2
CCReg_G2 <- read.table("~/Documents/THESIS/directories/signalome_construction/G2_Significant_Collapsed.txt",  header=T, sep="\t")
NCCReg_G2 <- read.table("~/Documents/THESIS/directories/signalome_construction/G2_Collapsed_NonCCRegulated.txt",  header=T, sep="\t")
cell_loc <- read_tsv("~/Documents/THESIS/directories/subcellular_location.tsv", show_col_types = FALSE)

# Match subcellular location for CCReg and NCCReg based on gene names
CCReg_G2$Gene <- toupper(CCReg_G2$Gene)
NCCReg_G2$Gene <- toupper(NCCReg_G2$Gene)

# Add location columns by matching gene names with cell_loc data
CCReg_G2 <- CCReg_G2 %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

NCCReg_G2 <- NCCReg_G2 %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

# Remove rows with NA in Main_location and separate multiple locations
CCReg_G2 <- CCReg_G2 %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

NCCReg_G2 <- NCCReg_G2 %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

# Calculate total counts for each Main_location in both datasets
total_counts_G2 <- CCReg_G2 %>%
  group_by(Main_location) %>%
  summarise(CC_count = n()) %>%
  left_join(NCCReg_G2 %>% group_by(Main_location) %>% summarise(NCC_count = n()), by = "Main_location") %>%
  mutate(
    total_count = CC_count + NCC_count,
    CC_percentage = (CC_count / total_count),
    NCC_percentage = (NCC_count / total_count),
    enrichment_difference = (CC_percentage - NCC_percentage) * 100
  )

# Perform Fisher's exact test for each subcellular location
total_counts_G2 <- total_counts_G2 %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(CC_count, NCC_count, total_count - CC_count, total_count - NCC_count), nrow = 2))$p.value,
    significance = ifelse(p_value < 0.05, TRUE, FALSE)
  )

# Filter to only significant subcellular locations
significant_counts_G2 <- total_counts_G2 %>%
  filter(significance == TRUE)

# Load datasets M
CCReg_M<- read.table("~/Documents/THESIS/directories/signalome_construction/M_Significant_Collapsed.txt",  header=T, sep="\t")
NCCReg_M <- read.table("~/Documents/THESIS/directories/signalome_construction/M_Collapsed_NonCCRegulated.txt",  header=T, sep="\t")
cell_loc <- read_tsv("~/Documents/THESIS/directories/subcellular_location.tsv", show_col_types = FALSE)

# Match subcellular location for CCReg and NCCReg based on gene names
CCReg_M$Gene <- toupper(CCReg_M$Gene)
NCCReg_M$Gene <- toupper(NCCReg_M$Gene)

# Add location columns by matching gene names with cell_loc data
CCReg_M <- CCReg_M %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

NCCReg_M <- NCCReg_M %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')])

# Remove rows with NA in Main_location and separate multiple locations
CCReg_M <- CCReg_M %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

NCCReg_M <- NCCReg_M %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")

# Calculate total counts for each Main_location in both datasets
total_counts_M <- CCReg_M %>%
  group_by(Main_location) %>%
  summarise(CC_count = n()) %>%
  left_join(NCCReg_M %>% group_by(Main_location) %>% summarise(NCC_count = n()), by = "Main_location") %>%
  mutate(
    total_count = CC_count + NCC_count,
    CC_percentage = (CC_count / total_count),
    NCC_percentage = (NCC_count / total_count),
    enrichment_difference = (CC_percentage - NCC_percentage) * 100
  )

# Perform Fisher's exact test for each subcellular location
total_counts_M <- total_counts_M %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(CC_count, NCC_count, total_count - CC_count, total_count - NCC_count), nrow = 2))$p.value,
    significance = ifelse(p_value < 0.05, TRUE, FALSE)
  )

# Filter to only significant subcellular locations
significant_counts_M <- total_counts_M %>%
  filter(significance == TRUE)

# Rename the enrichment_difference column for each dataset
significant_counts_G1 <- significant_counts_G1 %>%
  dplyr::rename(EG1vLG1 = enrichment_difference)

significant_counts_S <- significant_counts_S %>%
  dplyr::rename(EG1vS = enrichment_difference)

significant_counts_G2 <- significant_counts_G2 %>%
  dplyr::rename(EG1vG2 = enrichment_difference)

significant_counts_M <- significant_counts_M %>%
  dplyr::rename(EG1vM = enrichment_difference)

# Merge datasets by Main_location
merged_data <- significant_counts_G1 %>%
  dplyr::select(Main_location, EG1vLG1) %>%
  full_join(significant_counts_S %>% dplyr::select(Main_location, EG1vS), by = "Main_location") %>%
  full_join(significant_counts_G2 %>% dplyr::select(Main_location, EG1vG2), by = "Main_location") %>%
  full_join(significant_counts_M %>% dplyr::select(Main_location, EG1vM), by = "Main_location")

# Reshape data to long format for plotting
long_data <- merged_data %>%
  pivot_longer(cols = c(EG1vLG1, EG1vS, EG1vG2, EG1vM),
               names_to = "Phase",
               values_to = "Enrichment_Difference") %>%
  drop_na(Enrichment_Difference)  # Remove NA values for plotting

long_data <- long_data %>%
  mutate(Phase = factor(Phase, levels = c("EG1vLG1", "EG1vS", "EG1vG2", "EG1vM")))

library(ggplot2)
ggplot(long_data, aes(x = Main_location, y = Enrichment_Difference, color = Phase)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "Significant Percentage Difference of Phosphosite Enrichment by Subcellular Location",
    x = "Subcellular Location",
    y = "Cell Cycle Regulated - Non-Cell Cycle Regulated (% Difference)",
    color = "Phase Comparison"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



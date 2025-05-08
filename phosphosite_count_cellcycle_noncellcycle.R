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

# Load data. replace NCCREG with collapsed cell cycle or non cell cycle regulated data
NCCReg <- read_excel("~/Documents/THESIS/directories/signalome_construction/cell cycle collapsed data.xlsx")
cell_loc <- read_tsv("~/Documents/THESIS/directories/subcellular_location.tsv", show_col_types = FALSE)



# Match locations for NCCReg
NCCReg <- NCCReg %>%
  mutate(Gene = toupper(Gene)) %>%
  mutate(Main_location = cell_loc$'Main location'[match(Gene, cell_loc$'Gene name')],
         Additional_location = cell_loc$'Additional location'[match(Gene, cell_loc$'Gene name')],
         Extracellular_location = cell_loc$'Extracellular location'[match(Gene, cell_loc$'Gene name')]) %>%
  filter(!is.na(Main_location)) %>%
  separate_rows(Main_location, sep = ";")



# Reshape data for plotting - Non-cell cycle regulated (NCCReg)
data_long1 <- NCCReg %>%
  dplyr::select(Gene, Main_location, starts_with("phase")) %>%
  gather(key = "phase", value = "phospho_abundance", starts_with("phase")) %>%
  separate(phase, into = c("ignore", "Phase", "Sample"), sep = "_") %>%
  mutate(Phase = recode(Phase, "earlyG1" = "earlyG1", "lateG1" = "lateG1", "S" = "S", "G2" = "G2", "M" = "M")) %>%
  dplyr::select(-ignore)

data_long1$Phase <- factor(data_long1$Phase, levels = c("earlyG1", "lateG1", "S", "G2", "M"))



# Summarize total phosphosites by Main_location
phosphosite_counts_NCCR <- NCCReg %>%
  group_by(Main_location) %>%
  summarise(phosphosite_count = n(), .groups = "drop")

# Sort locations by phosphosite count for non-cell cycle regulated
phosphosite_counts_NCCR <- phosphosite_counts_NCCR %>%
  mutate(Main_location = factor(Main_location, levels = Main_location[order(phosphosite_count)]))

# Plot for non-cell cycle regulated phosphosites (ranked) chage title name based on cell cycle/ non cell cycle data
ggplot(phosphosite_counts_NCCR, aes(x = Main_location, y = phosphosite_count, fill = "Collapsed")) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("purple")) +
  scale_y_continuous(labels = scales::comma_format()) +  # Use regular numeric labels
  labs(
    title = "Phosphosite Counts by Subcellular Location",
    x = "Subcellular Location",
    y = "Phosphosite Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Removing current data
rm(list = ls())
# Importing clinical info
library(tidyverse)
TCGA_OV_CLINICAL_DF <- read.csv(file = "data/clinical/clinical.cart.2021-12-03/clinical.tsv", sep = "\t", header = T)
# Loading sample_sheet
TCGA_OV_SAMPLE_INFO <- read.csv(file = "data/clinical/gdc_sample_sheet.2021-12-03.tsv", sep = "\t", header = T)
# changing .gz to .txt
TCGA_OV_SAMPLE_INFO$File.Name <- gsub(".gz", ".txt", TCGA_OV_SAMPLE_INFO$File.Name)

# Removing dataframe 
if(exists("TCGA_OV_COUNTS_DF")){
  rm(TCGA_OV_COUNTS_DF) 
} 


# Building combined table
for (i in TCGA_OV_SAMPLE_INFO$File.Name) {
  # Reading file 
  print(paste0("Adding ", i, "  == ", round(match(i, TCGA_OV_SAMPLE_INFO$File.Name)*100/length(TCGA_OV_SAMPLE_INFO$File.Name), 2), "%"))
  file_path <- paste("./data/RNAseq/read_counts/", i, sep = "")
  input_df <- read.delim(file_path, header = F)
  
  if(exists("TCGA_OV_COUNTS_DF")){
    TCGA_OV_COUNTS_DF <- left_join(TCGA_OV_COUNTS_DF, input_df, by="V1")
  } else {
    TCGA_OV_COUNTS_DF <- input_df
  }
}



# Renaming columns
colnames(TCGA_OV_COUNTS_DF) <- c("ENSEMBL_ID", TCGA_OV_SAMPLE_INFO$File.Name)
# Removing version numbers from ENSEMBL IDs
TCGA_OV_COUNTS_DF$ENSEMBL_ID <- TCGA_OV_COUNTS_DF$ENSEMBL_ID %>%  
  str_replace(pattern = ".[0-9]+$",
              replacement = "") 


# Conducting TMM normalisation
source(file = "scripts/edgeR_TMM_count_normalisation.R")

# Creating annotation DF
# 1. Convert from ensembl.gene to gene.symbol
ensembl <- TCGA_OV_COUNTS_DF$ENSEMBL_ID %>% as.character()
# Creating DF for gene names matched with ensemble IDs
library(EnsDb.Hsapiens.v79)
TCGA_OV_GENE_ANNOT <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))


# Saving dfs
dir.create("./data/chris_project/", recursive = T, showWarnings = F)
write_delim(TCGA_OV_GENE_ANNOT, file = "./data/chris_project/TCGA_OV_GENE_ANNOT.txt")
write_delim(TCGA_OV_CLINICAL_DF, file = "./data/chris_project/TCGA_OV_SAMPLE_INFO.txt")
write_delim(TCGA_OV_SAMPLE_INFO, file = "./data/chris_project/TCGA_OV_SAMPLE_INFO.txt")
write_delim(TCGA_OV_COUNTS_DF, file = "./data/chris_project/TCGA_OV_COUNTS_DF.txt")
write_delim(TCGA_OV_CPM_DF, file = "./data/chris_project/TCGA_OV_TMM_CPM_DF.txt")
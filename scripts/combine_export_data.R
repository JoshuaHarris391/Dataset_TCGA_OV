# Setting wd
# setwd("~/Dropbox/Research/PhD/Bioinformatics/Datasets/Dataset_TCGA_BRCA/data/RNAseq")
library(tidyverse)
# getting list of files 
system("ls -1 ./data/RNAseq/read_counts > ./data/RNAseq/count_filenames.txt")

# loading list of files
count_filenames <- read.delim("./data/RNAseq/count_filenames.txt", header = F)
count_filenames <- count_filenames$V1 %>% as.character()

# Getting sample names
count_sample_names <- gsub(".htseq.counts.txt", "", count_filenames)

# Removing dataframe 
if(exists("TCGA_OV_COUNTS_DF")){
  rm(TCGA_OV_COUNTS_DF) 
} 


# Building combined table
for (i in count_filenames) {
  # Reading file 
  print(paste0("Adding ", i))
  file_path <- paste("./data/RNAseq/read_counts/", i, sep = "")
  input_df <- read.delim(file_path, header = F)
  
  if(exists("TCGA_OV_COUNTS_DF")){
    TCGA_OV_COUNTS_DF <- left_join(TCGA_OV_COUNTS_DF, input_df, by="V1")
  } else {
    TCGA_OV_COUNTS_DF <- input_df
  }
}

# Creating sample reference table
ID <- seq(1, length(count_sample_names)) %>% as.character()
ID <- paste0("ID_", ID)
CASE_ID <- count_sample_names
FILENAME <- count_filenames
ID_REF_TABLE <- data.frame(ID, CASE_ID, FILENAME)
dir.create("./data/RNAseq/combined_count_data/", recursive = T, showWarnings = F)
write_delim(ID_REF_TABLE, file = "./data/RNAseq/combined_count_data/ID_REF_TABLE.txt")

# Renaming columns
colnames(TCGA_OV_COUNTS_DF) <- c("ENSEMBL_ID", count_sample_names)
# Removing version numbers from ENSEMBL IDs
TCGA_OV_COUNTS_DF$ENSEMBL_ID <- TCGA_OV_COUNTS_DF$ENSEMBL_ID %>%  
  str_replace(pattern = ".[0-9]+$",
              replacement = "") 

# Saving count df
write_delim(TCGA_OV_COUNTS_DF, path = "./data/RNAseq/combined_count_data/TCGA_OV_COUNTS_DF.txt")

# Conducting TMM normalisation
source(file = "scripts/edgeR_TMM_count_normalisation.R")
# Saving count df
write_delim(TCGA_OV_CPM_DF, path = "./data/RNAseq/combined_count_data/TCGA_OV_TMM_CPM_DF.txt")

# Creating annotation DF
# 1. Convert from ensembl.gene to gene.symbol
ensembl <- TCGA_OV_COUNTS_DF$ENSEMBL_ID %>% as.character()
# Creating DF for gene names matched with ensemble IDs
library(EnsDb.Hsapiens.v79)
GENE_ANNOT_DF <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
# Saving annot df
write_delim(GENE_ANNOT_DF, path = "./data/RNAseq/combined_count_data/GENE_ANNOT_DF.txt")


# Saving Rdata
save(TCGA_OV_COUNTS_DF, TCGA_OV_CPM_DF, GENE_ANNOT_DF, ID_REF_TABLE,  file = "./data/RNAseq/combined_count_data/TCGA_OV_RNAseq.RData")


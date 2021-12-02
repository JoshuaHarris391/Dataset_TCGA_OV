# creating DGE object 
library(tidyverse)
library(edgeR)
TCGA_DGE <- edgeR::DGEList(counts=TCGA_OV_COUNTS_DF[, 2:ncol(TCGA_OV_COUNTS_DF)], genes=rownames(TCGA_OV_COUNTS_DF$ENSEMBL_ID))
# Calculating normalisation factors 
TCGA_DGE <- edgeR::calcNormFactors(TCGA_DGE, method="TMM")
# Calculating TMM normalised CPM
TCGA_OV_CPM_DF <- edgeR::cpm(TCGA_DGE) %>% as.data.frame()
rownames(TCGA_OV_CPM_DF) <- TCGA_OV_COUNTS_DF$ENSEMBL_ID

# References:
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# https://www.biostars.org/p/317701/
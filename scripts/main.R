# Processing Clinical data
source(file = "scripts/extracting_clinical.R")
# Downloading RNAseq data
system('bash scripts/download_tcga_ov_rnaseq.sh')
# Extracting count data to new folder
system('bash scripts/combine_counts.sh')
# Running combine data script
source("scripts/combine_export_data.R")

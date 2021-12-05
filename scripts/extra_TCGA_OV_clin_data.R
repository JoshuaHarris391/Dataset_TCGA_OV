library(curatedOvarianData)
data(package="curatedOvarianData")
data(TCGA_eset)
TCGA_eset$vital_status
featureNames(TCGA_eset)
sampleNames(TCGA_eset)
varLabels(TCGA_eset)
exp_mat <- exprs(TCGA_eset)

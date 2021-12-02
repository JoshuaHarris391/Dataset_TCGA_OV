#!/bin/bash
# defining path to gdc transfer tool
GDC_CLIENT=~/tools/gdc-client
# Setting wd
cd /Users/joshua_harris/Dropbox/Research/Bioinformatics/Datasets/Dataset_TCGA_OV
# Making dir
mkdir -p ./data/RNAseq/counts
# Changing dir
cd ./data/RNAseq/counts
# downloading files from manifest
$GDC_CLIENT download -m /Users/joshua_harris/Dropbox/Research/Bioinformatics/Datasets/Dataset_TCGA_OV/gdc_manifest/gdc_manifest.2021-11-26.txt

# This script is meant to download the emsembl dataset to the local machine,
# so that dataset updates do not make pathways dissappear
# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2022.07.2 Build 576
# R version: 4.2.2

library(biomaRt)

# Connect to the Ensembl biomart database
ensembl_mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Save the database to the local machine
save(ensembl_mart, file = "mmusculus_gene_ensembl.RData")
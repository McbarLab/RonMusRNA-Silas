# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2022.07.2 Build 576
# R version: 4.2.2

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")
if(!require("KEGGREST", quietly = TRUE))
  BiocManager::install("KEGGREST")

library(tidyverse)
library(org.Mm.eg.db)
library(KEGGREST)
library(AnnotationDbi)

pathway_dict <- function(species){
  pathways <- keggList("pathway", organism=species)
  pathway_dict <- setNames(names(pathways), pathways)
  return(pathway_dict)
}


# Extract a list of genes from a KEGG pathway
get_genes <- function(pathway_id){
  pathway <- keggGet(paste0("path:", pathway_id))
  genes <- mapIds(org.Mm.eg.db, keys=as.character(pathway[[1]][["GENE"]]),
                  keytype="ENTREZID", column="SYMBOL") %>% na.omit()
  return(genes)
}



mmu_dict <- pathway_dict("mmu")
ribosome_genes <- get_genes(
  mmu_dict["Ribosome - Mus musculus (house mouse)"])



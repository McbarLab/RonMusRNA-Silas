# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2022.12.0 Build 353
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

read_enrichment <- function(filename){
  dataset <- read_csv(filename, col_names = TRUE)
}

enrichment_compare <- function(df1, df2, df1_name, df2_name){
  merged_pathway <- merge(df1[,c("ID", "Description", "enrichmentScore", "gene_symbol")],
                          df2[,c("ID", "Description", "enrichmentScore", "gene_symbol")],
                          by = c("ID", "Description"),
                          suffixes = c(paste0("_",df1_name),
                                       paste0("_",df2_name)))
  write.csv(merged_pathway,
            file = paste("Enrichment_Comparison/",
                         title," enrichment_comparison.csv",sep=""))
}

heatmap_plot <- function(df1, df2, df1_name, df2_name, genes){
  groups <- c(df1_name, df2_name)
  heatmap_data <- expand.grid(X=groups, Y=genes)
  # Test-phrase, not real data
  data$count <- runif(440, 0, 6)
  heatmap_draw <- ggplot(heatmap_data, aes(X, Y, fill=count) +
    geom_tile() + theme_minimal() +
    scale_fill_gradient(low="white", high="red")+
    labs(title = "Heatmap", subtitle = "A simple heatmap using geom_tile()",
         x ="Group", y ="Genes"))
  return(heatmap_draw)
}

mmu_dict <- pathway_dict("mmu")
pathway_genes <- get_genes(
  mmu_dict["Ribosome - Mus musculus (house mouse)"])

enrichment_folder_name <- "GSEA_Pathway_csv"
comparison_list <- c(
  "01_mo06_AR_F_vs_mo06_C_F","03_mo22_AR_F_vs_mo22_C_F",
  "03_mo22_AR_F_vs_mo22_C_F","05_mo28_AR_F_vs_mo28_C_F"
  )
pair_count <- length(comparison_list) / 2

for(pair_index in 1:pair_count){
  df1_name <- comparison_list[2*pair_index - 1]
  df2_name <- comparison_list[2*pair_index]
  df1 <- read_enrichment(paste0(enrichment_folder_name, "/",
                                df1_name, " gseaKEGG_all.csv"))
  df2 <- read_enrichment(paste0(enrichment_folder_name, "/",
                                df2_name, " gseaKEGG_all.csv"))
  
  title <- paste0(df1_name,"_vs_", df2_name)
  enrichment_compare(df1, df2, df1_name, df2_name)
}

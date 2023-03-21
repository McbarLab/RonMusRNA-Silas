# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2022.12.0 Build 353
# R version: 4.2.2
# IMPORTANT: Restart R before running this script, if you have used any other script in this repo

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if(!require("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")
if(!require("EnhancedVolcano", quietly = TRUE))
  BiocManager::install('EnhancedVolcano')
if(!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(AnnotationDbi)
library(biomaRt)


# Local database is used to ensure pathway does not change due to updates
emsembl_mart <- load("mmusculus_gene_ensembl.RData")

import_dataset <- function(filename){
  rawDGE <- read_csv(filename)
  # The following line is for Mark Berres dataset
  curatedDGE <- rawDGE[,c(1,2,5,6)] %>% na.omit()
  colnames(curatedDGE) <- c("Symbol", "logFC", "pval", "FDR")
  
  # Add a column named direction, to show whether this gene is UP-regulated
  # or DOWN-regulated
  # Take into consideration of the FDR. Do not simply look at the sign.
  curatedDGE <- curatedDGE %>% 
    mutate(direction = case_when(FDR < 0.05 & logFC > 0 ~ "UP",
                                 FDR < 0.05 & logFC < 0 ~ "DOWN",
                                 FDR >= 0.05 ~ "NS")) 
  return(curatedDGE)
}

category_plot <- function(curatedDGE, title){
  totExpressed <- curatedDGE %>% 
    tally() %>% 
    mutate(direction = "all trans")
  sigDirection <- curatedDGE %>% 
    group_by(direction) %>% 
    tally() 
  allSig <- curatedDGE %>%
    filter(direction != "NS") %>%
    tally() %>% 
    mutate(direction = "normal Sig")
  restrictSig <- curatedDGE %>% 
    filter(FDR <= 0.01) %>% 
    tally() %>% 
    mutate(direction = "restrictive Sig")
  
  # Plot transcript categories
  transcriptBreakdown <- bind_rows(
    totExpressed, sigDirection, allSig, restrictSig)
  transcriptNumsPlot <- transcriptBreakdown %>% 
    arrange(match(direction, c("all trans", "DOWN", "UP", 
                               "normal Sig", "NS", 
                               "restrictive Sig")), 
            desc(direction)) %>% 
    mutate(direction = factor(direction, levels = direction)) %>% 
    ggplot(aes(x = direction, y = n, 
               fill = direction)) + geom_col(color = "black", 
                                             size = 0.25) +
    xlab("Transcript Category") + ylab("count") + theme_bw() + 
    theme(axis.line = element_line(color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") + 
    # 10000 is the number of rows in curatedDGE
    scale_y_continuous(expand = c(0,0), limits = c(0, nrow(curatedDGE)*1.1)) + 
    geom_text(aes(label = n), vjust = -0.5) + 
    scale_fill_manual(values = c("black", "#808080", 
                                 "#524fa1", "#fdb913", 
                                 "red", "cyan"))
  ggsave(path = "./Transcript_Categories",
         filename = paste(title,"transcriptCategories.pdf",sep=" "), 
         transcriptNumsPlot, 
         height = 4, width = 4)
}

volcano_plot <- function(curatedDGE, title){
  volcanoPlot <- 
    curatedDGE %>% EnhancedVolcano(
      lab = curatedDGE$Symbol,
      x = 'logFC',
      # IMPORTANT: Do NOT put logFDR here, it takes negative log by default
      y = 'FDR', 
      xlab = expression(paste("Log"[2],"(FC)")),
      ylab = expression(paste("-log"[10],"(FDR)")),
      FCcutoff = 1.5,
      cutoffLineType = 'twodash',
      # col = c("#999999", "#ed1c24", "#662d91"),
      # pointSize = 3.0,
      labSize = 4.0,
      boxedLabels = TRUE,
      colAlpha = 4/5,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      colConnectors = 'grey50',
      # Since dot colors are easy to understand, hide legends
      legendPosition = 'none',
      title = paste(title,"Volcano",sep=" "),
      subtitle = '',
      caption = '',
      max.overlaps = 100
    )
  ggsave(path = "./Volcano_Plots",
         filename = paste(title,"Volcano.pdf",sep=" "), 
         volcanoPlot)
}

sig_gene <- function(curatedDGE, title){
  # Extract restrictively significant genes
  sigGenes <- curatedDGE[,c("Symbol", "logFC","FDR")] %>% 
    filter(FDR < 0.01) %>% 
    # Calculate the non-log fold change and filter out rows with values smaller than 50%
    # log2(1.5) is about 0.5849, we take 0.58 here
    filter(abs(logFC)>= 0.58)
  
  write.csv(sigGenes, file = paste("Significant_Genes/",title,"sigGenes.csv",sep=" "))
  return(sigGenes)
}

path_generate <- function(curatedDGE, title){
  # Map the EMSEMBL IDs to their ENTREZID
  sameGenesEntrez <- mapIds(org.Mm.eg.db, 
                            keys = as.character(curatedDGE$Symbol),
                            keytype = "SYMBOL", column = "ENTREZID") %>% na.omit()
  # Function for enrichKEGG analysis and pathway dataset generation
  # By default, this works on mouse (mmu) only
  allPaths <- enrichKEGG(gene = sameGenesEntrez, organism = "mmu")
  #write.csv(allPaths, file = paste("All_Pathways/","allPaths_",title,".csv",sep = ""))
  
  allPaths <- as.data.frame(allPaths) %>%
    mutate(bkgdSize = 
             as.numeric(substring(BgRatio, 
                                  regexpr("/", BgRatio) + 1))) %>%
    mutate(pathBkgd = 
             as.numeric(substring(BgRatio, 1,
                                  regexpr("/", BgRatio)-1))) %>%
    mutate(bkgdPerc = pathBkgd/bkgdSize) %>%
    mutate(GeneRatTotal = 
             as.numeric(substring(GeneRatio, 
                                  regexpr("/", GeneRatio) + 1))) %>%
    mutate(percPath = Count/GeneRatTotal) %>%
    mutate(Enrichment = percPath/bkgdPerc)
  return(allPaths)
} 

ORA_plot <- function(pathway, title){
  plot_pathway <- pathway %>% 
    filter(p.adjust < 0.01) %>% 
    arrange(desc(Count)) %>%
    slice(1:10) %>%
    ggplot(aes(x = Enrichment, y = Description, 
               color = p.adjust, size = Count)) + 
    geom_point() + expand_limits(x = 0) + 
    labs(x = "Enrichment", y = "KEGG pathway", 
         color = "FDR", size = "Count") +
    theme_bw() + scale_color_gradient(low = "#B72668", 
                                      high = "#dba3b2") + 
    ggtitle(paste(title,"enrichKEGG", sep = " "))
  
  ggsave(paste("ORA_Pathways/",title,"enrichKEGG.pdf",sep=" "), plot_pathway)
}

GSEA_plot <- function(curatedDGE, title){
  GSEA_genes <- curatedDGE$logFC
  # GSEA can take in both Ensembl or Symbol
  names(GSEA_genes) <- curatedDGE$Symbol
  # ENTREZ ids for more accurate result
  GSEA_id <- bitr(names(GSEA_genes), fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # Remove duplicated entries
  dedup_id <- GSEA_id[!duplicated(GSEA_id[c("SYMBOL")]), ]
  dedup_id <- dedup_id[!duplicated(dedup_id[c("ENTREZID")]), ]
  
  # Remove genes that have more than 1 rows
  # Note: This removes ALL rows for those genes from the dataset
  # Extract these deduplicated entries from the original curatedDGE
  mappedIDs <- curatedDGE %>% 
    filter(!duplicated(Symbol)) %>% 
    filter(Symbol %in% dedup_id$SYMBOL)
  mappedIDs$entrez <- dedup_id$ENTREZID
  
  # Further distill the data by removing nulls
  kegg_genes <- mappedIDs$logFC
  names(kegg_genes) <- mappedIDs$entrez
  kegg_genes <- kegg_genes %>% na.omit()
  kegg_genes <- sort(kegg_genes, decreasing = T)
  
  # Set the timeout limit for gseKEGG() and useMart()
  options(timeout = Inf)
  
  # This is the very step of running GSEA, takes VERY LONG TIME
  gsea_result <- gseKEGG(geneList = kegg_genes, 
                         organism = "mmu",
                         minGSSize = 4, 
                         maxGSSize = 500, 
                         pvalueCutoff = 0.05, 
                         pAdjustMethod = "fdr",
                         keyType = "ncbi-geneid",
                         verbose = TRUE)

  # Output all pathways as text into a csv file
  # Translate the Entrez IDs into gene symbols for easier interpretation
  
  # Get the gene symbols and add them to gsea_result@result
  for (i in seq_along(gsea_result@result$core_enrichment)) {
    input_string = gsea_result@result$core_enrichment[i]
    input_ids = strsplit(input_string, "/")[[1]]
    output_ids = "mgi_symbol"
    
    # Get the gene symbols
    gsea_genes = getBM(attributes = c(output_ids), 
                       filters = c("entrezgene_id"), 
                       values = input_ids, 
                       mart = ensembl_mart,
                       )
    
    # Collapse the gene symbols into a single string separated by backslashes
    gene_string = paste(gsea_genes$mgi_symbol, collapse = "/")
    
    # Assign the gene symbol to the corresponding row of gsea_result@result
    gsea_result@result[i, "gene_symbol"] <- gene_string
  }
  
  # Write the dataframe into an external csv file
  write.csv(gsea_result@result,
            file = paste("GSEA_Pathway_csv/", title," gseaKEGG_top10.csv",sep=""),
            row.names = FALSE)
  
  gsea_dotplot <- dotplot(gsea_result, 
                          showCategory = 10,
                          title = paste(title,"gseaKEGG_top10",sep=" "), 
                          split = ".sign") + facet_grid(.~.sign)
  
  ggsave(path = "./GSEA_Pathways",
         filename = paste(title,"gseaKEGG_top10.pdf",sep=" "), 
         gsea_dotplot)
}

# Write all DGE spreadsheet names into a list
DGE_folder_name <- "Joe_DGE"
DGE_list <- list.files(path = paste("./",DGE_folder_name,sep = ""))
DGE_count <- length(DGE_list)

# Iterate through all files to generate their figures in 1 click
for(DGE_index in 1:DGE_count){
  title <- sub(".csv*.","",
                        sub(".*GENE_","",DGE_list[DGE_index]))
  curatedDGE <- import_dataset(paste(DGE_folder_name, DGE_list[DGE_index],sep="/"))
  category_plot(curatedDGE, title)
  volcano_plot(curatedDGE, title)
  sig_gene_list <- sig_gene(curatedDGE, title)
  allPaths <- path_generate(curatedDGE, title)
  ORA_plot(allPaths, title)
  GSEA_plot(curatedDGE, title)
}

# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu

if(!require(clusterProfiler))
  BiocManager::install("clusterProfiler")
if(!require(org.Mm.eg.db))
  BiocManager::install("org.Mm.eg.db")
if(!require(EnhancedVolcano))
  BiocManager::install('EnhancedVolcano')

library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)

import_dataset <- function(filename){
  rawDGE <- read_csv(filename)
  curatedDGE <- rawDGE[,c(2,3,4,7,8)]
  colnames(curatedDGE) <- c("Ensembl", "Symbol", "logFC", "pval", "FDR")
  
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
                            keys = as.character(curatedDGE$Ensembl),
                            keytype = "ENSEMBL", column = "ENTREZID") %>% na.omit()
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
  names(GSEA_genes) <- curatedDGE$Ensembl
  # ENTREZ ids for more accurate result
  GSEA_id <- bitr(names(GSEA_genes), fromType = "ENSEMBL",
                  toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # Remove duplicated entries
  dedup_id <- GSEA_id[!duplicated(GSEA_id[c("ENSEMBL")]), ]
  dedup_id <- dedup_id[!duplicated(dedup_id[c("ENTREZID")]), ]
  
  # Extract these deduplicated entries from the original curatedDGE
  mappedIDs <- curatedDGE[curatedDGE$Ensembl %in% 
                            dedup_id$ENSEMBL, ]
  mappedIDs$entrez <- dedup_id$ENTREZID
  
  # Further distill the data by removing nulls
  kegg_genes <- mappedIDs$logFC
  names(kegg_genes) <- mappedIDs$entrez
  kegg_genes <- kegg_genes %>% na.omit()
  kegg_genes <- sort(kegg_genes, decreasing = T)
  
  # This is the very step of running GSEA, takes VERY LONG TIME
  gsea_result <- gseKEGG(geneList = kegg_genes, 
                         organism = "mmu",
                         # IMPORTANT: Discuss the GSSize param with Tim!!
                         minGSSize = 4, 
                         maxGSSize = 500, 
                         pvalueCutoff = 0.05, 
                         pAdjustMethod = "fdr",
                         keyType = "ncbi-geneid")
  
  # Output all pathways as text into a csv file
  write.csv(gsea_result@result, file = paste("GSEA_Pathway_csv/", title," gseaKEGG_top10.csv",sep=""))
  
  gsea_dotplot <- dotplot(gsea_result, 
                          showCategory = 10, 
                          title = paste(title,"gseaKEGG_top10",sep=" "), 
                          split = ".sign") + facet_grid(.~.sign)
  
  ggsave(path = "./GSEA_Pathways",
         filename = paste(title,"gseaKEGG_top10.pdf",sep=" "), 
         gsea_dotplot)
}

# Write all DGE spreadsheet names into a list
DGE_folder_name <- "Corrected_DGE"
DGE_list <- list.files(path = paste("./",DGE_folder_name,sep = ""))
DGE_count <- length(DGE_list)

# Iterate through all files to generate their figures in 1 click
for(DGE_index in 1:DGE_count){
  title <- sub(".csv*.","",
               sub(".*edgeRglm_GENE_","",DGE_list[DGE_index]))
  curatedDGE <- import_dataset(paste(DGE_folder_name, DGE_list[DGE_index],sep="/"))
  category_plot(curatedDGE, title)
  volcano_plot(curatedDGE, title)
  sig_gene_list <- sig_gene(curatedDGE, title)
  allPaths <- path_generate(curatedDGE, title)
  ORA_plot(allPaths, title)
  GSEA_plot(curatedDGE, title)
}

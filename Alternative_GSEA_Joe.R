# Modified Joe's version of GSEA
# Deprecated by Tim's more elegant method
GSEA_plot <- function(curatedDGE, title){
  msig_mart <- msigdbr(species = "Mus musculus",
                       category = "C2", subcategory = "CP:KEGG")
  curatedDGE_Entrez <- curatedDGE %>% 
    mutate(Entrez = ensembl_to_entrez(curatedDGE))
  entrez_list <- msig_mart %>%
    dplyr::select(gs_name, entrez_gene) %>%
    group_by(gs_name) %>%
    summarise(all.genes = list(unique(entrez_gene))) %>%
    deframe()
  entrez_logFC <- curatedDGE_Entrez %>% 
    dplyr::select(Entrez, logFC) %>%
    na.omit(entrez_logFC$Entrez)
  
  # Form a vector(required format of input) for GSEA
  entrez_logFC.vec <- entrez_logFC$logFC
  names(entrez_logFC.vec) <- entrez_logFC$Entrez
  entrez_logFC.vec <- na.omit(entrez_logFC.vec)
  
  # Run GSEA
  entrez_logFC.output <- fgseaSimple(pathway = entrez_list,
                                     stats = entrez_logFC.vec,
                                     scoreType = "std",
                                     nperm = 10^10)
  
  # Plot the figure
  figure_plot <- as.data.frame(entrez_logFC.output) %>% 
    filter(p.adjust < 1*10^-20) %>% 
    ggplot(aes(x = Enrichment, y = Description, 
               color = p.adjust, size = Count)) + 
    geom_point() + expand_limits(x = 0) + 
    labs(x = "Enrichment", y = "KEGG pathway", 
         color = "FDR", size = "Count") +
    theme_bw() + scale_color_gradient(low = "#B72668", 
                                      high = "#dba3b2") + 
    ggtitle(paste(title,"KEGG pathways p < 1*10^-20", sep = " "))
  ggsave(path = "./GSEA_Pathways",
         filename = paste(title,"KEGG pathways p less than 0.01.pdf",sep=" "), 
         figure_plot)
}
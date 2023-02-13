# Silas's original ORA
# Deprecated by a more advanced and elegant version
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
  pdf(paste("ORA_Pathways/",title,"enrichKEGG_top10.pdf",sep=" "))  
  plot(pathway %>% 
         filter(p.adjust < 1*10^-20) %>% 
         ggplot(aes(x = Enrichment, y = Description, 
                    color = p.adjust, size = Count)) + 
         geom_point() + expand_limits(x = 0) + 
         labs(x = "Enrichment", y = "KEGG pathway", 
              color = "FDR", size = "Count") +
         theme_bw() + scale_color_gradient(low = "#B72668", 
                                           high = "#dba3b2") + 
         ggtitle(paste(title,"enrichKEGG_top10.pdf", sep = " ")))
  dev.off()
}
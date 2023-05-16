
  
  
  
  
  
  
  ##Network analysis is next.
  # Network analysis using WGCNA.
  #Use log2 transformed/z-scored data for this analysis: redImpGeneTPM.df
  
  
  
  #The following packages will be necessary for the downstream analysis and
  #visualization of this data.
  if (!requireNamespace("BiocManager", quietly = T))
    install.packages("BiocManager")
  if (!require(RColorBrewer))
    install.packages("RColorBrewer")
  if (!require(tidyverse))
    install.packages("tidyverse")
  if (!require(gplots))
    install.packages("gplots")
  
  library(BiocManager)
  library(RColorBrewer)
  library(tidyverse)
  library(gplots)
  
  if (!require(WGCNA))
    BiocManager::install("WGCNA")
  #This is for Mouse
  ##if(!require(org.Mm.eg.db))
  ##BiocManager::install("org.Mm.eg.db")
  
  if (!require(org.Hs.eg.db))
    BiocManager::install("org.Hs.eg.db")
  library(WGCNA)
  library(org.Hs.eg.db)
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("KEGG.db")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("AnnotationDbi")
  
  if (!require(clusterProiler))
    BiocManager::install("clusterProfiler")
  library(clusterProfiler)
  
  
  
  
  #Data integration. First step here will be collecting the gene membership of the modules
  #determined by WGCNA and then correlating with lipid composition.
  #Relies on data structures created in the NetworkAnalysis.R script (which in turn relies on the
  #RNAscript.R). Will also rely on the LipidScript.R
  
  
  if (!require(bnstruct))
    install.packages("bnstruct")
  library(bnstruct)
  
  if (!require(showtext))
    install.packages("showtext")
  library(showtext)
  
  font_add_google("Poppins", "Poppins")
  showtext_auto()
  
  
  
  ###__________________
  # 4.a Output gene lists for usei with online software
  # Read in the probe annotation
  annot = read.csv(file = "RawTPM_All_18288genes_AllSamples.csv")
  
  # Match probes in the data set to the probe IDs in the annotation file
  probes = names(datExpr)
  probes2annot = match(probes, annot$gene)
  # Get the corresponding Locuis Link IDs
  allLLIDs = annot$LocusLinkID[probes2annot]
  
  # $ Choose interesting modules
  intModules = c("turquoise", "red", "black", "green", "skyblue", "tan")
  for (module in intModules)
  {
    # Select module probes
    modGenes = (moduleColors == module)
    # Get their entrez ID codes
    modLLIDs = allLLIDs[modGenes]
    
    # Write them into a file
    fileName = paste("AllSamples_All18288genes_Remove626genes_CRvsC-",
                     module,
                     ".txt",
                     sep = "")
    
    write.table(
      as.data.frame(modLLIDs),
      file = fileName,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  # As background in the enrichment analysis, we will use all probes in the analysis.
  fileName = paste("AllSamples_All18288genes_Remove626genes_CRvsC.txt", sep =
                     "")
  
  write.table(
    as.data.frame(allLLIDs),
    file = fileName,
    row.names = FALSE,
    col.names = FALSE
  )
  
  ###___________
  #Gene Module membership
  MMredGenes <- names(datExpr)[moduleColors == "red"]
  MMcyanGenes <- names(datExpr)[moduleColors == "cyan"]
  MMblueGenes <- names(datExpr)[moduleColors == "blue"]
  MMpurpleGenes <- names(datExpr)[moduleColors == "purple"]
  MMplum1Genes <- names(datExpr)[moduleColors == "plum1"]
  MMyellowGenes <- names(datExpr)[moduleColors == "yellow"]
  MMdarkorangeGenes <- names(datExpr)[moduleColors == "darkorange"]
  MMturquoiseGenes <- names(datExpr)[moduleColors == "turquoise"]
  MMsienna3Genes <- names(datExpr)[moduleColors == "sienna3"]
  MMdarkslateblueGenes <-
    names(datExpr)[moduleColors == "darkslateblue"]
  MMbisque4Genes <- names(datExpr)[moduleColors == "bisque4"]
  
  
  MMmediumpurple3Genes <-
    names(datExpr)[moduleColors == "mediumpurple3"]
  MMbisque4Genes <- names(datExpr)[moduleColors == "bisque4"]
  MMdarkorangeGenes <- names(datExpr)[moduleColors == "darkorange"]
  MMorangered4Genes <- names(datExpr)[moduleColors == "orangered4"]
  MMdarkolivegreenGenes <-
    names(datExpr)[moduleColors == "darkolivegreen"]
  MMskyblue3Genes <- names(datExpr)[moduleColors == "skyblue3"]
  MMyellowGenes <- names(datExpr)[moduleColors == "yellow"]
  MMlightcyan1Genes <- names(datExpr)[moduleColors == "lightcyan1"]
  MMplum1Genes <- names(datExpr)[moduleColors == "plum1"]
  MMgreenyellowGenes <-
    names(datExpr)[moduleColors == "greenyellow"]
  MMdarkmagentaGenes <-
    names(datExpr)[moduleColors == "darkmagenta"]
  MMpurpleGenes <- names(datExpr)[moduleColors == "purple"]
  MMfloralwhiteGenes <-
    names(datExpr)[moduleColors == "floralwhite"]
  MMvioletGenes <- names(datExpr)[moduleColors == "violet"]
  MMdarkslateblueGenes <-
    names(datExpr)[moduleColors == "darkslateblue"]
  MMpaleturquoiseGenes <-
    names(datExpr)[moduleColors == "paleturquoise"]
  MMturquoiseGenes <- names(datExpr)[moduleColors == "turquoise"]
  MMtanGenes <- names(datExpr)[moduleColors == "tan"]
  MMblackGenes <- names(datExpr)[moduleColors == "black"]
  MMblueGenes <- names(datExpr)[moduleColors == "blue"]
  MMdarkgreyGenes <- names(datExpr)[moduleColors == "darkgrey"]
  MMsienna3Genes <- names(datExpr)[moduleColors == "sienna3"]
  MMdarkredGenes <- names(datExpr)[moduleColors == "darkred"]
  MMcyanGenes <- names(datExpr)[moduleColors == "cyan"]
  MMredGenes <- names(datExpr)[moduleColors == "red"]
  MMyellowgreenGenes <-
    names(datExpr)[moduleColors == "yellowgreen"]
  MMlightcyanGenes <- names(datExpr)[moduleColors == "lightcyan"]
  MMlightsteelblue1Genes <-
    names(datExpr)[moduleColors == "lightsteelblue1"]
  MMgreyGenes <- names(datExpr)[moduleColors == "grey"]
  
  
  #Pathway Analysis for gene modules
  
  ##mediumpurple3 Module. Contains 157 gene members.
  modmediumpurple3Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMmediumpurple3Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modmediumpurple3Entrez <- na.exclude(modmediumpurple3Entrez)
  modmediumpurple3Paths <-
    enrichKEGG(gene = modmediumpurple3Entrez, organism = 'hsa')
  modmediumpurple3Paths <-
    as.data.frame(modmediumpurple3Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                       1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                                 GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modmediumpurple3Paths, file = "modmediumpurple3Paths_blockwise.csv")
  modmediumpurple3Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modmediumpurple3KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##bisque4 Module. Contains 528 gene members.
  modbisque4Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMbisque4Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modbisque4Entrez <- na.exclude(modbisque4Entrez)
  modbisque4Paths <-
    enrichKEGG(gene = modbisque4Entrez, organism = 'hsa')
  modbisque4Paths <-
    as.data.frame(modbisque4Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                 1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                           GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modbisque4Paths, file = "modbisque4Paths_blockwise.csv")
  modbisque4Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modbisque4KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##darkorange Module. Contains 1662 gene members.
  moddarkorangeEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMdarkorangeGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  moddarkorangeEntrez <- na.exclude(moddarkorangeEntrez)
  moddarkorangePaths <-
    enrichKEGG(gene = moddarkorangeEntrez, organism = 'hsa')
  moddarkorangePaths <-
    as.data.frame(moddarkorangePaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                    1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                              GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(moddarkorangePaths, file = "moddarkorangePaths_blockwise.csv")
  moddarkorangePaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "moddarkorangeKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##orangered4 Module. Contains 726 gene members.
  modorangered4Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMorangered4Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modorangered4Entrez <- na.exclude(modorangered4Entrez)
  modorangered4Paths <-
    enrichKEGG(gene = modorangered4Entrez, organism = 'hsa')
  modorangered4Paths <-
    as.data.frame(modorangered4Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                    1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                              GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modorangered4Paths, file = "modorangered4Paths_blockwise.csv")
  modorangered4Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modorangered4KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##darkolivegreen Module. Contains 210 gene members.
  moddarkolivegreenEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMdarkolivegreenGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  moddarkolivegreenEntrez <- na.exclude(moddarkolivegreenEntrez)
  moddarkolivegreenPaths <-
    enrichKEGG(gene = moddarkolivegreenEntrez, organism = 'hsa')
  moddarkolivegreenPaths <-
    as.data.frame(moddarkolivegreenPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                        1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                                  GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(moddarkolivegreenPaths, file = "moddarkolivegreenPaths_blockwise.csv")
  moddarkolivegreenPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "moddarkolivegreenKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##skyblue3 Module. Contains 188 gene members.
  modskyblue3Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMskyblue3Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modskyblue3Entrez <- na.exclude(modskyblue3Entrez)
  modskyblue3Paths <-
    enrichKEGG(gene = modskyblue3Entrez, organism = 'hsa')
  modskyblue3Paths <-
    as.data.frame(modskyblue3Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                  1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                            GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modskyblue3Paths, file = "modskyblue3Paths_blockwise.csv")
  modskyblue3Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modskyblue3KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##yellow Module. Contains 687 gene members.
  modyellowEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMyellowGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modyellowEntrez <- na.exclude(modyellowEntrez)
  modyellowPaths <-
    enrichKEGG(gene = modyellowEntrez, organism = 'hsa')
  modyellowPaths <-
    as.data.frame(modyellowPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                          GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modyellowPaths, file = "modyellowPaths_blockwise.csv")
  modyellowPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modyellowKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##lightcyan1 Module. Contains 138 gene members.
  modlightcyan1Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMlightcyan1Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modlightcyan1Entrez <- na.exclude(modlightcyan1Entrez)
  modlightcyan1Paths <-
    enrichKEGG(gene = modlightcyan1Entrez, organism = 'hsa')
  modlightcyan1Paths <-
    as.data.frame(modlightcyan1Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                    1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                              GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modlightcyan1Paths, file = "modlightcyan1Paths_blockwise.csv")
  modlightcyan1Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modlightcyan1KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##plum1 Module. Contains 175 gene members.
  modplum1Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMplum1Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modplum1Entrez <- na.exclude(modplum1Entrez)
  modplum1Paths <-
    enrichKEGG(gene = modplum1Entrez, organism = 'hsa')
  modplum1Paths <-
    as.data.frame(modplum1Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                               1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                         GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modplum1Paths, file = "modplum1Paths_blockwise.csv")
  modplum1Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modplum1KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##greenyellow Module. Contains 1051 gene members.
  modgreenyellowEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMgreenyellowGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modgreenyellowEntrez <- na.exclude(modgreenyellowEntrez)
  modgreenyellowPaths <-
    enrichKEGG(gene = modgreenyellowEntrez, organism = 'hsa')
  modgreenyellowPaths <-
    as.data.frame(modgreenyellowPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                     1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                               GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modgreenyellowPaths, file = "modgreenyellowPaths_blockwise.csv")
  modgreenyellowPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modgreenyellowKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##darkmagenta Module. Contains 1019 gene members.
  moddarkmagentaEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMdarkmagentaGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  moddarkmagentaEntrez <- na.exclude(moddarkmagentaEntrez)
  moddarkmagentaPaths <-
    enrichKEGG(gene = moddarkmagentaEntrez, organism = 'hsa')
  moddarkmagentaPaths <-
    as.data.frame(moddarkmagentaPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                     1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                               GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(moddarkmagentaPaths, file = "moddarkmagentaPaths_blockwise.csv")
  moddarkmagentaPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "moddarkmagentaKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##purple Module. Contains 1892 gene members.
  modpurpleEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMpurpleGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modpurpleEntrez <- na.exclude(modpurpleEntrez)
  modpurplePaths <-
    enrichKEGG(gene = modpurpleEntrez, organism = 'hsa')
  modpurplePaths <-
    as.data.frame(modpurplePaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                          GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modpurplePaths, file = "modpurplePaths_blockwise.csv")
  modpurplePaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modpurpleKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##floralwhite Module. Contains 119 gene members.
  modfloralwhiteEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMfloralwhiteGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modfloralwhiteEntrez <- na.exclude(modfloralwhiteEntrez)
  modfloralwhitePaths <-
    enrichKEGG(gene = modfloralwhiteEntrez, organism = 'hsa')
  modfloralwhitePaths <-
    as.data.frame(modfloralwhitePaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                     1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                               GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modfloralwhitePaths, file = "modfloralwhitePaths_blockwise.csv")
  modfloralwhitePaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modfloralwhiteKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##violet Module. Contains 219 gene members.
  modvioletEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMvioletGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modvioletEntrez <- na.exclude(modvioletEntrez)
  modvioletPaths <-
    enrichKEGG(gene = modvioletEntrez, organism = 'hsa')
  modvioletPaths <-
    as.data.frame(modvioletPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                          GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modvioletPaths, file = "modvioletPaths_blockwise.csv")
  modvioletPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modvioletKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##darkslateblue Module. Contains 1020 gene members.
  moddarkslateblueEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMdarkslateblueGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  moddarkslateblueEntrez <- na.exclude(moddarkslateblueEntrez)
  moddarkslatebluePaths <-
    enrichKEGG(gene = moddarkslateblueEntrez, organism = 'hsa')
  moddarkslatebluePaths <-
    as.data.frame(moddarkslatebluePaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                       1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                                 GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(moddarkslatebluePaths, file = "moddarkslatebluePaths_blockwise.csv")
  moddarkslatebluePaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "moddarkslateblueKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##paleturquoise Module. Contains 236 gene members.
  modpaleturquoiseEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMpaleturquoiseGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modpaleturquoiseEntrez <- na.exclude(modpaleturquoiseEntrez)
  modpaleturquoisePaths <-
    enrichKEGG(gene = modpaleturquoiseEntrez, organism = 'hsa')
  modpaleturquoisePaths <-
    as.data.frame(modpaleturquoisePaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                       1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                                 GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modpaleturquoisePaths, file = "modpaleturquoisePaths_blockwise.csv")
  modpaleturquoisePaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modpaleturquoiseKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##turquoise Module. Contains 1331 gene members.
  modturquoiseEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMturquoiseGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modturquoiseEntrez <- na.exclude(modturquoiseEntrez)
  modturquoisePaths <-
    enrichKEGG(gene = modturquoiseEntrez, organism = 'hsa')
  modturquoisePaths <-
    as.data.frame(modturquoisePaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                   1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                             GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modturquoisePaths, file = "modturquoisePaths_blockwise.csv")
  modturquoisePaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modturquoiseKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##tan Module. Contains 463 gene members.
  modtanEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMtanGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modtanEntrez <- na.exclude(modtanEntrez)
  modtanPaths <- enrichKEGG(gene = modtanEntrez, organism = 'hsa')
  modtanPaths <-
    as.data.frame(modtanPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                             1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                       GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modtanPaths, file = "modtanPaths_blockwise.csv")
  modtanPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modtanKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##black Module. Contains 958 gene members.
  modblackEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMblackGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modblackEntrez <- na.exclude(modblackEntrez)
  modblackPaths <-
    enrichKEGG(gene = modblackEntrez, organism = 'hsa')
  modblackPaths <-
    as.data.frame(modblackPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                               1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                         GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modblackPaths, file = "modblackPaths_blockwise.csv")
  modblackPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modblackKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##blue Module. Contains 1371 gene members.
  modblueEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMblueGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modblueEntrez <- na.exclude(modblueEntrez)
  modbluePaths <- enrichKEGG(gene = modblueEntrez, organism = 'hsa')
  modbluePaths <-
    as.data.frame(modbluePaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                              1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                        GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modbluePaths, file = "modbluePaths_blockwise.csv")
  modbluePaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modblueKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##darkgrey Module. Contains 316 gene members.
  moddarkgreyEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMdarkgreyGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  moddarkgreyEntrez <- na.exclude(moddarkgreyEntrez)
  moddarkgreyPaths <-
    enrichKEGG(gene = moddarkgreyEntrez, organism = 'hsa')
  moddarkgreyPaths <-
    as.data.frame(moddarkgreyPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                  1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                            GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(moddarkgreyPaths, file = "moddarkgreyPaths_blockwise.csv")
  moddarkgreyPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "moddarkgreyKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##sienna3 Module. Contains 201 gene members.
  modsienna3Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMsienna3Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modsienna3Entrez <- na.exclude(modsienna3Entrez)
  modsienna3Paths <-
    enrichKEGG(gene = modsienna3Entrez, organism = 'hsa')
  modsienna3Paths <-
    as.data.frame(modsienna3Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                 1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                           GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modsienna3Paths, file = "modsienna3Paths_blockwise.csv")
  modsienna3Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modsienna3KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##darkred Module. Contains 353 gene members.
  moddarkredEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMdarkredGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  moddarkredEntrez <- na.exclude(moddarkredEntrez)
  moddarkredPaths <-
    enrichKEGG(gene = moddarkredEntrez, organism = 'hsa')
  moddarkredPaths <-
    as.data.frame(moddarkredPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                 1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                           GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(moddarkredPaths, file = "moddarkredPaths_blockwise.csv")
  moddarkredPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "moddarkredKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##cyan Module. Contains 1054 gene members.
  modcyanEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMcyanGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modcyanEntrez <- na.exclude(modcyanEntrez)
  modcyanPaths <- enrichKEGG(gene = modcyanEntrez, organism = 'hsa')
  modcyanPaths <-
    as.data.frame(modcyanPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                              1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                        GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modcyanPaths, file = "modcyanPaths_blockwise.csv")
  modcyanPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modcyanKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##red Module. Contains 586 gene members.
  modredEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMredGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modredEntrez <- na.exclude(modredEntrez)
  modredPaths <- enrichKEGG(gene = modredEntrez, organism = 'hsa')
  modredPaths <-
    as.data.frame(modredPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                             1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                       GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modredPaths, file = "modredPaths_blockwise.csv")
  modredPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##yellowgreen Module. Contains 193 gene members.
  modyellowgreenEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMyellowgreenGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modyellowgreenEntrez <- na.exclude(modyellowgreenEntrez)
  modyellowgreenPaths <-
    enrichKEGG(gene = modyellowgreenEntrez, organism = 'hsa')
  modyellowgreenPaths <-
    as.data.frame(modyellowgreenPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                     1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                               GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modyellowgreenPaths, file = "modyellowgreenPaths_blockwise.csv")
  modyellowgreenPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modyellowgreenKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##lightcyan Module. Contains 656 gene members.
  modlightcyanEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMlightcyanGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modlightcyanEntrez <- na.exclude(modlightcyanEntrez)
  modlightcyanPaths <-
    enrichKEGG(gene = modlightcyanEntrez, organism = 'hsa')
  modlightcyanPaths <-
    as.data.frame(modlightcyanPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                   1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                             GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modlightcyanPaths, file = "modlightcyanPaths_blockwise.csv")
  modlightcyanPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modlightcyanKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##lightsteelblue1 Module. Contains 148 gene members.
  modlightsteelblue1Entrez <-
    mapIds(org.Hs.eg.db,
           keys = MMlightsteelblue1Genes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modlightsteelblue1Entrez <- na.exclude(modlightsteelblue1Entrez)
  modlightsteelblue1Paths <-
    enrichKEGG(gene = modlightsteelblue1Entrez, organism = 'hsa')
  modlightsteelblue1Paths <-
    as.data.frame(modlightsteelblue1Paths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                                         1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                                   GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modlightsteelblue1Paths, file = "modlightsteelblue1Paths_blockwise.csv")
  modlightsteelblue1Paths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modlightsteelblue1KEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  ##grey Module. Contains 5 gene members.
  modgreyEntrez <-
    mapIds(org.Hs.eg.db,
           keys = MMgreyGenes,
           keytype = "SYMBOL",
           column = "ENTREZID")
  modgreyEntrez <- na.exclude(modgreyEntrez)
  modgreyPaths <- enrichKEGG(gene = modgreyEntrez, organism = 'hsa')
  modgreyPaths <-
    as.data.frame(modgreyPaths) %>% mutate(bkgdSize = as.numeric(substring(BgRatio, regexpr("/", BgRatio) + 1))) %>% mutate(pathBkgd = as.numeric(substring(BgRatio, 1, regexpr("/", BgRatio) -
                                                                                                                                                              1))) %>% mutate(bkgdPerc = pathBkgd / bkgdSize) %>%     mutate(GeneRatTotal = as.numeric(substring(GeneRatio, regexpr("/", GeneRatio) + 1))) %>% mutate(percPath = Count /
                                                                                                                                                                                                                                                                                                                        GeneRatTotal) %>% mutate(Enrichment = percPath / bkgdPerc)
  write.csv(modgreyPaths, file = "modgreyPaths_blockwise.csv")
  modgreyPaths %>% ggplot(aes(
    x = Enrichment,
    y = Description,
    color = p.adjust,
    size = Count
  )) + geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                                  y = "KEGG pathway",
                                                  color = "FDR",
                                                  size = "Count") + theme_bw() + scale_color_gradient(low = "#000000", high = "#ffffff")
  ggsave(filename = "modgreyKEGGpaths_blockwise.pdf", useDingbats = F)
  ggsave(filename = "modredKEGGpaths_blockwise.pdf", useDingbats = F)
  
  modKEGGpaths_blockwise <- rbind(
    modbisque4Paths,
    modblackPaths,
    modbluePaths,
    modcyanPaths,
    moddarkgreyPaths,
    moddarkmagentaPaths,
    moddarkolivegreenPaths,
    moddarkorangePaths,
    moddarkslatebluePaths,
    modfloralwhitePaths,
    modgreenyellowPaths,
    modgreyPaths,
    modlightcyanPaths,
    modplum1Paths,
    modpurplePaths,
    modredPaths,
    modsienna3Paths,
    modskyblue3Paths,
    modtanPaths,
    modturquoisePaths,
    modvioletPaths,
    modyellowPaths
  )
  write.csv(modKEGGpaths_blockwise, file = "modKEGGtotalPaths_blockwise.csv")
  
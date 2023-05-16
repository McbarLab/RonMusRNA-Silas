

# 2.a Automatic network construction and module detection
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(
  sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",
  type = "n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.90, col = "red")
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)

# 2.a.2 One-step network construction and module detection
## Setting Soft Threshold to 9
net = blockwiseModules(
  datExpr,
  power = 9,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AllSamples_CRvsC_TOM", # Change this title
  verbose = 3
)

## See how many modules and genes/module
table(net$colors)

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs

geneTree = net$dendrograms[[1]]

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "AllSamples_networkConstruction-auto.RData")




###________________________________________


# 3.a Quantifying module-trait associations
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10, 6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2),
                   "\n(",
                   signif(moduleTraitPvalue, 1),
                   ")",
                   sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3)) # Discuss on how parameters were determined

# Display the correlation values within a heatmap plot
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)



### _________________

# Define variable weight containing the weight column of datTrait
## Fasting Glucose End
glu_end = as.data.frame(datTraits$Fasting.Glucose.End)

names(glu_end) = "glu_end"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep = "")

names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(datExpr, glu_end, use = "p"))

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(glu_end), sep =
                                       "")

names(GSPvalue) = paste("p.GS.", names(glu_end), sep = "")

# 3.c Intramodular analysis: identifying genes with high GS and MM
module = "orange"
  column = match(module, modNames)
  
  moduleGenes = moduleColors == module
  
  sizeGrWindow(7, 7)
  
  par(mfrow = c(1, 1))
  
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = "Gene significance for triglycerides",
    main = paste("Module membership vs. gene significance\n"),
    cex.main = 1.2,
    cex.lab = 1.2,
    cex.axis = 1.2,
    col = module
  )
  
  # 3.d Summary output of network analysis results
  names(datExpr)
  names(datExpr)[moduleColors == "orange"]
  annot = read.csv(file = "RawTPM_All_18288genes_AllSamples.csv")
  
  dim(annot)
  names(annot)
  probes = names(datExpr)
  probes2annot = match(probes, annot$gene)
  # The following is the number or probes without annotation:
  sum(is.na(probes2annot))
  # Should return 0.
  
  # Create the starting data frame
  geneInfo0 = data.frame(gene = probes,
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  
  # Order modules by their significance for weight
  modOrder = order(-abs(cor(MEs, triglycerides, use = "p")))
  
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]])
    
    names(geneInfo0) = c(oldNames,
                         paste("MM.", modNames[modOrder[mod]], sep = ""),
                         paste("p.MM.", modNames[modOrder[mod]], sep =
                                 ""))
  }
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.triglycerides))
  
  geneInfo = geneInfo0[geneOrder, ]
  
  ## Remove specific columns to condense data
  geneInfo_Si = geneInfo[, -c(5:86)] #this removes these specific columns
  dim(geneInfo_Si)
  names(geneInfo_Si)
  
  write.csv(geneInfo, file = "AllSamples_All18288genes_Remove626genes_CRvsC_geneInfo_blockwise.csv")
  geneInfo_Gene_Gene = geneInfo
  geneInfo_Gene_Gene$gene = paste0(geneInfo_Gene_Gene$gene, "_")
  write.csv(geneInfo_Gene_Gene, file = "AllSamples_All18288genes_Remove626genes_CRvsC_geneInfo_blockwise.csv")
  
  # 4.a Output gene lists for usei with online software
  # Read in the probe annotation
  annot = read.csv(file = "RawTPM_All_18288genes_AllSamples.csv")
  
  # Match probes in the data set to the probe IDs in the annotation file
  probes = names(datExpr)
  probes2annot = match(probes, annot$gene)
  # Get the corresponding Locuis Link IDs
  allLLIDs = annot$LocusLinkID[probes2annot]
  
  # $ Choose interesting modules
  intModules = c("turquoise", "red", "black")
  for (module in intModules)
  {
    # Select module probes
    modGenes = (moduleColors == module)
    # Get their entrez ID codes
    modLLIDs = allLLIDs[modGenes]
    
    # Write them into a file
    fileName = paste("LocusLinkIDs-", module, ".txt", sep = "")
    
    write.table(
      as.data.frame(modLLIDs),
      file = fileName,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  # As background in the enrichment analysis, we will use all probes in the analysis.
  fileName = paste("LocusLinkIDs-all.txt", sep = "")
  
  write.table(
    as.data.frame(allLLIDs),
    file = fileName,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # 5.a Visualizing the gene network
  # Calculate topological overlap anew: this could be done more efficiently by saving the TOM
  # calculated during module detection, but let us do it again here.
  dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = 6)
  
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM ^ 7
  
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA
  
  # Call the plot function
  sizeGrWindow(9, 9)
  TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  
  ## Do this one to minimize time of computation
  dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = 9)
  
  nSelect = 400
  # For reproducibility, we set the random seed
  set.seed(10)
  
  select = sample(nGenes, size = nSelect)
  
  selectTOM = dissTOM[select, select]
  
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select]
  
  # Open a graphical window
  sizeGrWindow(9, 9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM ^ 7
  
  diag(plotDiss) = NA
  
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
  
  # 5.b Visualizing the network of eigengenes
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  # Isolate weight from the clinical traits
  triglycerides = as.data.frame(datTraits$Triglycerides_mgPERdl)
  
  names(triglycerides) = "triglycerides"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, triglycerides))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5, 7.5)
  
  par(cex = 0.9)
  plotEigengeneNetworks(
    MET,
    "",
    marDendro = c(0, 4, 1, 2),
    marHeatmap = c(3, 4, 1, 2),
    cex.lab = 0.8,
    xLabelsAngle
    = 90
  )
  # Plot the dendrogram
  sizeGrWindow(6, 6)
  
  par(cex = 1.0)
  plotEigengeneNetworks(
    MET,
    "Eigengene dendrogram",
    marDendro = c(0, 4, 2, 0),
    plotHeatmaps = FALSE
  )
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  plotEigengeneNetworks(
    MET,
    "Eigengene adjacency heatmap",
    marHeatmap = c(3, 4, 2, 2),
    plotDendrograms = FALSE,
    xLabelsAngle = 90
  )
  
  
  
  
  
  
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
  
load("01-dataInput.RData")
load("02_networkConstr.RData")

# 3.a Quantifying module-trait associations
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = WGCNA::cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2),
                   "\n(",
                   signif(moduleTraitPvalue, 1),
                   ")",
                   sep = "")

dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
pdf("Module-trait relationships.pdf", width = 20, height = 12)
par(mar = c(6, 8.5, 3, 3))
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
dev.off()

feature_plot <- function(feature){
  feature_df <- as.data.frame(datTraits$feature)
  #names(feature_df) <- sub(".*\\$", "", as.character(feature))
  names(feature_df) <- deparse(substitute(feature))
  return(feature_df)
}

test_feature <- feature_plot(Fasting.Glucose.End)

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
  ylab = "Gene significance for fasting glucose end",
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = module
)

# 3.d Summary output of network analysis results
names(datExpr)
names(datExpr)[moduleColors == "orange"]
annot = read.csv(file = "rsubread_GENE_tpm.csv")

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
modOrder = order(-abs(cor(MEs, glu_end, use = "p")))

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
geneOrder = order(geneInfo0$moduleColor,-abs(geneInfo0$GS.glu_end))

geneInfo = geneInfo0[geneOrder,]

## Remove specific columns to condense data
geneInfo_Si = geneInfo[,-c(5:86)] #this removes these specific columns
dim(geneInfo_Si)
names(geneInfo_Si)

write.csv(geneInfo, file = "geneInfo.csv")
geneInfo_Gene = geneInfo
geneInfo_Gene$gene = paste0(geneInfo_Gene$gene, "_")
write.csv(geneInfo_Gene, file = "geneInfo_genes.csv")
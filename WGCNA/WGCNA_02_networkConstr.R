source("WGCNA_01_dataInput.R")

load("01-dataInput.RData")

cor <- WGCNA::cor

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

sft_power <- sft$powerEstimate
# 2.a.2 One-step network construction and module detection
## Setting Soft Threshold to a suitable number
net = blockwiseModules(
  datExpr,
  power = sft_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  # saveTOMFileBase = "AllSamples_TOM",
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

save(MEs, moduleLabels, moduleColors, geneTree, sft_power,
     file = "02_networkConstr.RData")

cor <- stats::cor
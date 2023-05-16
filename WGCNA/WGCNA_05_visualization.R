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
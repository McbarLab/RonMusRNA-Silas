source("WGCNA_03_relateModsToExt.R")

load("01-dataInput.RData")
load("02_networkConstr.RData")

# 5.b Visualizing the network of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

eigenAdjHeatmap_plot <- function(trait) {
  # Isolate a trait from the clinical traits
  trait_df <- as.data.frame(datTraits[, trait])
  names(trait_df) <- trait
  # Add the trait to existing module eigengenes
  MET <- orderMEs(cbind(MEs, trait_df))
  
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5, 7.5)
  
  #par(cex = 0.9)
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
  pdf(
    paste0("./EigenNetwork_Plot/", trait, ".pdf"),
    width = 6,
    height = 6
  )
  
  #par(cex = 1.0)
  plotEigengeneNetworks(
    MET,
    "Eigengene dendrogram",
    marDendro = c(0, 4, 2, 0),
    plotHeatmaps = FALSE
  )
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  #par(cex = 1.0)
  plotEigengeneNetworks(
    MET,
    "Eigengene adjacency heatmap",
    marHeatmap = c(3, 4, 2, 2),
    plotDendrograms = FALSE,
    xLabelsAngle = 90
  )
  dev.off()
}

trait_list <- names(datTraits)

# Parallelize code with 4 physical cores, which is all i7-6700 has
eigen_cluster <- makeCluster(4)
registerDoParallel(eigen_cluster)

# Iterate through all traits
foreach(i = seq_along(trait_list), .packages = c("WGCNA")) %dopar% {
  eigenAdjHeatmap_plot(trait_list[i])
}

stopCluster(eigen_cluster)
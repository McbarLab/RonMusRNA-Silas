source("WGCNA_00_setEnvir.R")
load("01-dataInput.RData")
load("02_networkConstr.RData")

# 3.a Quantifying module-trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)



# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Retrive a list of all traits and one of all MEs
trait_list <- names(datTraits)
MEs_list <- names(MEs)


moduleTraitCor = WGCNA::cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2),
                   "\n(",
                   signif(moduleTraitPvalue, 1),
                   ")",
                   sep = "")

dim(textMatrix) = dim(moduleTraitCor)

heatmap_plot <- function(matrix, title){
  pdf(paste0("./", title,".pdf"),
      width = 20,
      height = 12)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(
    Matrix = matrix,
    xLabels = trait_list,
    yLabels = MEs_list,
    ySymbols = MEs_list,
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = title
  )
  dev.off()
}

# Display the correlation values within a heatmap plot
heatmap_plot(moduleTraitCor, "Module-trait relationships")

# Set 2 empty tables to store correlation estimates and p values
cor_table <- data.frame(replace(moduleTraitCor, TRUE, NA))
p_table <- cor_table

# Extract names of all modules for iterations
modNames = substring(MEs_list, 3)
num_module <- length(modNames)

trait_cor <- function(trait) {
  # Define variable weight containing the weight column of datTrait
  trait_df <- as.data.frame(datTraits[,trait])
  names(trait_df) <- trait
  
  # names (colors) of the modules
  # modNames, MEs, datExpr, nSamples are global variables
  
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep = "")
  names(MMPvalue) = paste("p.MM", modNames, sep = "")
  geneTraitSignificance = as.data.frame(cor(datExpr, trait_df, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(trait_df), sep =
                                         "")
  names(GSPvalue) = paste("p.GS.", names(trait_df), sep = "")
  
  # 3.c Intramodular analysis: identifying genes with high GS and MM
  for (i in 1:num_module) {
    module <- modNames[i]
    moduleGenes <- moduleColors == module
    cor_test <- cor.test(abs(geneModuleMembership[moduleGenes, i]),
                         abs(geneTraitSignificance[moduleGenes, 1]))
    ME_name <- paste0("ME", module)
    cor_table[ME_name, trait] <<- cor_test$estimate
    p_table[ME_name, trait] <<- cor_test$p.value
  }
}

# Iterate through all traits
for (i in 1:length(trait_list)){
  trait_cor(trait_list[i])
  print(paste0("Correlation analysis for ", trait_list[i], " is done..."))
  if(i==length(trait_list)){
    print("Correlation anallyses all done!")
  }
}

heatmap_plot(cor_table, "Module_Trait_Correlation")

write.csv(cor_table, "Correlation_Estimate.csv")
write.csv(p_table, "Correlation_PValue.csv")

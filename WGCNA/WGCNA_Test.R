trait <- "Sex.M.vs.F"
trait_df <- as.data.frame(datTraits[,trait])
names(trait_df) <- trait

# names (colors) of the modules
# modNames, MEs, datExpr, nSamples are global variables


geneTraitSignificance = as.data.frame(cor(datExpr, trait_df, use = "p"))
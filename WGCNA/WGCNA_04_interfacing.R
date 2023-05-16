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
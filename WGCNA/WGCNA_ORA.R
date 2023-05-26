# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2023.03.0 Build 386
# R version: 4.3.0
source("WGCNA_05_visualization.R")


load("01-dataInput.RData")
load("02_networkConstr.RData")
load("mmusculus_gene_ensembl.RData")

fdrCutOff = 0.05

# Extract names of genes
gene_list <- colnames(datExpr)

# Create an empty list with all modules as keys
gene_module <- setNames(vector("list", length(modNames)), modNames)

# Iterate through the gene list and module list
for (i in seq_along(gene_list)) {
  gene <- gene_list[i]
  module <- moduleColors[i]
  
  # Append the gene to the corresponding module entry
  gene_module[[module]] <- c(gene_module[[module]], gene)
}

# Iterate through the vector to run ORA on all modules
for (j in seq_along(gene_module)) {
  module_name <- names(gene_module)[j]
  module_gene <- gene_module[[module_name]]
  module_entrez <- mapIds(
    org.Mm.eg.db,
    keys = as.character(module_gene),
    keytype = "SYMBOL",
    column = "ENTREZID"
  ) %>% na.omit()
  module_ora <- ora_result <- enrichKEGG(
    gene = module_entrez,
    organism = "mmu",
    minGSSize = 4,
    maxGSSize = 500,
    pvalueCutoff = fdrCutOff,
    pAdjustMethod = "fdr"
  )
  
  # Save the result into csv files
  if (is.null(module_ora@result) == FALSE) {
    # Write the dataframe into an external csv file
    write.csv(
      module_ora@result,
      file = paste("./ORA_Pathway_csv/", module_name, " ORA.csv", sep = ""),
      row.names = TRUE
    )
  }
}

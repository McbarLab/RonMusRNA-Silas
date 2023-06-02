# This analysis is coded by Di "Silas" Kuang
# Email: dkuang5@wisc.edu
# RStudio version: 2023.03.0 Build 386
# R version: 4.3.0
source("WGCNA_05_visualization.R")


load("01-dataInput.RData")
load("02_networkConstr.RData")
load("03_relatedModsToExt.RData")
load("mmusculus_gene_ensembl.RData")

# Set the start time
start_time <- Sys.time()

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


# Parallelize code with 4 physical cores, which is all i7-6700 has
ora_cluster <- makeCluster(4)
registerDoParallel(ora_cluster)

# Iterate through the vector to run ORA on all modules
foreach(
  j = seq_along(gene_module),
  .packages = c("org.Mm.eg.db", "clusterProfiler", "biomaRt", "doParallel")
) %dopar% {
  module_name <- names(gene_module)[j]
  module_gene <- gene_module[[module_name]]
  module_entrez <- mapIds(
    org.Mm.eg.db,
    keys = as.character(module_gene),
    keytype = "SYMBOL",
    column = "ENTREZID"
    
  ) %>% na.omit()
  module_ora <- enrichKEGG(
    gene = module_entrez,
    organism = "mmu",
    minGSSize = 4,
    maxGSSize = 500,
    pvalueCutoff = fdrCutOff,
    pAdjustMethod = "fdr"
  )
  
  module_result <- module_ora@result
  
  # Save the result into csv files
  if (is.null(module_result) == FALSE) {
    # Get the gene symbols and add them to gsea_result@result
    foreach(k = seq_along(module_result$geneID),
            .packages = c("biomaRt")) %do% {
              input_string = module_result$geneID[k]
              input_ids = strsplit(input_string, "/")[[1]]
              output_ids = "mgi_symbol"
              
              # Get the gene symbols
              ora_genes = getBM(
                attributes = c(output_ids),
                filters = c("entrezgene_id"),
                values = input_ids,
                mart = ensembl_mart,
              )
              
              # Collapse the gene symbols into a single string separated by backslashes
              gene_string = paste(ora_genes$mgi_symbol, collapse = "/")
              
              # Assign the gene symbol to the corresponding row of gsea_result@result
              module_result[k, "gene_symbol"] <- gene_string
            }
    # Write the dataframe into an external csv file
    write.csv(
      module_result,
      file = paste("./ORA_Pathway_csv/", module_name, " ORA.csv", sep = ""),
      row.names = TRUE
    )
  }
}

stopCluster(ora_cluster)

save(gene_module,
     file = "ORA_geneModule.RData")


# Calculate the elapsed time
end_time <- Sys.time()
elapsed_time <- end_time - start_time

# Print the elapsed time
cat("Elapsed time:", format(elapsed_time, units = "mins"), "\n")

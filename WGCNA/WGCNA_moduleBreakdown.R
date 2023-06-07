load("01-dataInput.RData")
load("02_networkConstr.RData")
load("03_relatedModsToExt.RData")
load("ORA_geneModule.RData")

module_interest <- c("brown")

tpm_table <- t(datExpr)

for (i in seq_along(module_interest)) {
  module <- module_interest[i]
  module_genes <- gene_module[[module]]
  module_breakdown <-
    tpm_table[rownames(tpm_table) %in% module_genes, , drop = FALSE]
  normalized_breakdown <-
    t(apply(module_breakdown, 1, function(x)
      base::scale(x, center = TRUE, scale = TRUE)))
  
  colnames(normalized_breakdown) <- colnames(module_breakdown)
  distance_matrix <-
    dist(normalized_breakdown, method = "euclidean")
  hclust_result <- hclust(distance_matrix)
  sorted_breakdown <- normalized_breakdown[hclust_result$order, ]
  sorted_breakdown <- sorted_breakdown[,-AR_list]
  # heatmap(sorted_breakdown,
  #         col = colorRampPalette(c("blue", "white", "red"))(100),
  #         Colv = NA)
  par(mar = c(1, 1, 1, 1))
  heatmap.2(x = sorted_breakdown,
            Colv=FALSE, 
            dendrogram="row",
            scale="row",
            col="bluered",
            trace="none",
            density.info = "none",
            labRow=FALSE,
            key = TRUE
            )
}

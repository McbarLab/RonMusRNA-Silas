source("WGCNA_ORA.R")

load("01-dataInput.RData")
load("02_networkConstr.RData")
load("03_relatedModsToExt.RData")
load("ORA_geneModule.RData")

AR_list <-
  c(4:6, 10:12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 31, 33, 36)

module_interest <- c("brown")

tpm_table <- t(datExpr)[,-AR_list]


for (i in seq_along(module_interest)) {
  module <- module_interest[i]
  module_genes <- gene_module[[module]]
  module_breakdown <-
    tpm_table[rownames(tpm_table) %in% module_genes, , drop = FALSE]
  normalized_breakdown <-
    t(apply(module_breakdown, 1, function(x)
      scale(x, center = TRUE, scale = TRUE)))
  colnames(normalized_breakdown) <- colnames(module_breakdown)
  distance_matrix <-
    dist(normalized_breakdown, method = "euclidean")
  hclust_result <- hclust(distance_matrix)
  sorted_breakdown <- normalized_breakdown[hclust_result$order, ]
  heatmap(sorted_breakdown,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          Colv = NA)
}

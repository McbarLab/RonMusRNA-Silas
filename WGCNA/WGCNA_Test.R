module <- module_interest[1]
module_genes <- gene_module[[module]]
module_breakdown <-
  tpm_table[rownames(tpm_table) %in% module_genes, , drop = FALSE]
normalized_breakdown <-
  scale(module_breakdown, center = TRUE, scale = TRUE)
distance_matrix <- dist(normalized_breakdown, method = "euclidean")
hclust_result <- hclust(distance_matrix)
sorted_breakdown <- normalized_breakdown[hclust_result$order, ]

heatmap(sorted_breakdown, col = colorRampPalette(c("blue", "white", "red"))(100), Colv = NA)

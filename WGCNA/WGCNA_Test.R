counter <- 0
for(i in seq_along(gene_module)){
  counter <- counter + length(gene_module[[i]])
}
print(counter)
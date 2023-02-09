# This is Silas's personal version of the analysis
library(tidyverse)

if(!require(org.Mm.eg.db))
  BiocManager::install("org.Mm.eg.db")
if(!require(WGCNA))
  BiocManager::install("WGCNA")


library(readxl)
library(clusterProfiler)
library(pathfindR)
library(org.Mm.eg.db)
library(WGCNA)


#Previous analysis DGE (via Mark Berres and edgeR)
oMdrug.dge <- read_excel("DGE/edgeRglm_GENE_s28mo_C_M-s28mo_AR_M.xlsx")
mMdrug.dge <- read_excel("DGE/edgeRglm_GENE_s22mo_C_M-s22mo_AR_M.xlsx")
yMdrug.dge <- read_excel("DGE/edgeRglm_GENE_s06mo_C_M-s06mo_AR_M.xlsx")

oFdrug.dge <- read_excel("DGE/edgeRglm_GENE_s28mo_C_F-s28mo_AR_F.xlsx")
mFdrug.dge <- read_excel("DGE/edgeRglm_GENE_s22mo_C_F-s22mo_AR_F.xlsx")
yFdrug.dge <- read_excel("DGE/edgeRglm_GENE_s06mo_C_F-s06mo_AR_F.xlsx")

#TPM (via Mark Berres and RSEM)

allTpm_cols <- c("Ensembl", "G25-06mo-C-M",
                 "G26-06mo-C-M", "G27-06mo-C-M",
                 "J33-06mo-AR-M", "J34-06mo-AR-M",
                 "J35-06mo-AR-M", "C1-06mo-C-F",
                 "C2-06mo-C-F", "C11-06mo-C-F",
                 "YF3-06mo-AR-F", "YF10-06mo-AR-F",
                 "YF12-06mo-AR-F", "M4-22mo-C-M",
                 "M7-22mo-AR-M", "M12-22mo-C-M",
                 "M15-22mo-AR-M", "M19-22mo-C-M", 
                 "M23-22mo-AR-M", "F3-22mo-C-F",
                 "F7-22mo-AR-F", "F12-22mo-C-F",
                 "F17-22mo-AR-F", "F27-22mo-C-F",
                 "F28-22mo-AR-F", "M33-28mo-C-M",
                 "M35-28mo-AR-M", "M37-28mo-C-M",
                 "M50-28mo-AR-M", "M65-28mo-C-M",
                 "M72-28mo-AR-M", "F32-28mo-AR-F",
                 "F41-28mo-C-F", "F47-28mo-AR-F",
                 "F49-28mo-C-F", "F62-28mo-C-F",
                 "F66-28mo-AR-F")

allTpm <- read_table("rsem_GENE_TPM.txt", col_names = allTpm_cols, skip = 1)


#Note - these fold changes were calculated with the control in the numerator,
#so they go in the opposite direction of what is typically expected.

colnames(oMdrug.dge) <- c("Ensembl", "symbol", "logFC", "logCPM", "LR", "pval",
                          "FDR", "sig", "descript")
colnames(mMdrug.dge) <- c("Ensembl", "symbol", "logFC", "logCPM", "LR", "pval",
                          "FDR", "sig", "descript")
colnames(yMdrug.dge) <- c("Ensembl", "symbol", "logFC", "logCPM", "LR", "pval",
                          "FDR", "sig", "descript")

colnames(oFdrug.dge) <- c("Ensembl", "symbol", "logFC", "logCPM", "LR", "pval",
                          "FDR", "sig", "descript")
colnames(mFdrug.dge) <- c("Ensembl", "symbol", "logFC", "logCPM", "LR", "pval",
                          "FDR", "sig", "descript")
colnames(yFdrug.dge) <- c("Ensembl", "symbol", "logFC", "logCPM", "LR", "pval",
                          "FDR", "sig", "descript")

tidy_yMdrug.dge <- yMdrug.dge[,-c(4:5,9)] %>% 
  mutate(logFDR = -log10(FDR))

yMdrug_volc <- tidy_yMdrug.dge %>% 
  ggplot(aes(x = logFC, y = logFDR, color = sig)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("blue", "gray", "red")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,25)) + 
  ggtitle("Young Male C vs AR") + theme(legend.position = c(0.82, 0.8),
                                        legend.background = 
                                          element_rect(fill = "white",
                                                       color = "black"))
plot(yMdrug_volc)

ggsave("yMdrg_volc.pdf", yMdrug_volc, height = 4, width = 4)

tidy_yFdrug.dge <- yFdrug.dge[,-c(4:5,9)] %>% 
  mutate(logFDR = -log10(FDR))

yFdrug_volc <- tidy_yFdrug.dge %>% 
  ggplot(aes(x = logFC, y = logFDR, color = sig)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("blue", "gray", "red")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,17)) +
  scale_x_continuous(expand = c(0,0), limits = c(-3,4)) +
  ggtitle("Young Female C vs AR") + theme(legend.position = c(0.83, 0.82),
                                          legend.background = 
                                            element_rect(fill = "white", 
                                                         color = "black"))
plot(yFdrug_volc)

ggsave("yFdrg_volc.pdf", yFdrug_volc, height = 4, width = 4)

tidy_mMdrug.dge <- mMdrug.dge[,-c(4:5,9)] %>% 
  mutate(logFDR = -log10(FDR), remove = F)

mMdrug_volc <- tidy_mMdrug.dge %>% 
  ggplot(aes(x = logFC, y = logFDR, color = sig)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("blue", "gray", "red")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6)) + 
  ggtitle("Middle Male C vs AR")

plot(mMdrug_volc)

ggsave("mMdrg_volc.pdf", mMdrug_volc, height = 4, width = 4)

tidy_mFdrug.dge <- mFdrug.dge[,-c(4:5,9)] %>% 
  mutate(logFDR = -log10(FDR), remove = F)

mFdrug_volc <- tidy_mFdrug.dge %>% 
  ggplot(aes(x = logFC, y = logFDR, color = sig)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("blue", "gray", "red")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,10)) +
  scale_x_continuous(limits = c(-3,3)) +
  ggtitle("Middle Female C vs AR") + theme(legend.position = c(0.87, 0.87),
                                          legend.background = 
                                            element_rect(fill = "white", 
                                                         color = "black"))

plot(mFdrug_volc)

ggsave("mFdrug_volc.pdf", mFdrug_volc, height = 4, width = 4)

tidy_oMdrug.dge <- oMdrug.dge[,-c(4:5,9)] %>% 
  mutate(logFDR = -log10(FDR), remove = F)

oMdrug_volc <- tidy_oMdrug.dge %>% 
  ggplot(aes(x = logFC, y = logFDR, color = sig)) + geom_point() + theme_bw() + 
  scale_color_manual(values = c("blue", "gray", "red")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5)) +
  scale_x_continuous(limits = c(-3,3)) +
  ggtitle("Old Male C vs AR")

plot(oMdrug_volc)

ggsave("oMdrug_volc.pdf", oMdrug_volc, height = 4, width = 4)

tidy_oFdrug.dge <- oFdrug.dge[,-c(4:5,9)] %>% 
  mutate(logFDR = -log10(FDR), remove = F)

oFdrug_volc <- tidy_oFdrug.dge %>% 
  ggplot(aes(x = logFC, y = logFDR, color = sig)) + geom_point() + theme_bw() + 
  scale_color_manual(values = c("blue", "gray", "red")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,2.5)) +
  #scale_x_continuous(limits = c(-3,3)) +
  ggtitle("Old Female C vs AR")

plot(oFdrug_volc)

ggsave("oFdrug_volc.pdf", oFdrug_volc, height = 4, width = 4)


#conventional pathway analysis
#males
DE_yMdrug <- tidy_yMdrug.dge %>% 
  filter(FDR < 0.05)

DE_yMdrug_entrez <- mapIds(org.Mm.eg.db, keys = DE_yMdrug$Ensembl,
                           keytype = "ENSEMBL", column = "ENTREZID")
DE_yMdrug_entrez <- na.exclude(DE_yMdrug_entrez)

DE_yMdrug_paths <- enrichKEGG(gene = DE_yMdrug_entrez, organism = "mmu",
                              pAdjustMethod = "BH", minGSSize = 4,
                              maxGSSize = 500, qvalueCutoff = 0.05)
write.csv(DE_yMdrug_paths, "YmalePaths.csv")
#qvalue < 0.007

DE_yMdrug_paths <- as.data.frame(DE_yMdrug_paths) %>% 
  mutate(bkgdSize = as.numeric(substring(BgRatio, 
                                         regexpr("/", BgRatio) + 1))) %>%
  mutate(pathBkgd = as.numeric(substring(BgRatio, 1, 
                                         regexpr("/", BgRatio)-1))) %>%
  mutate(bkgdPerc = pathBkgd/bkgdSize) %>%
  mutate(GeneRatTotal = as.numeric(substring(GeneRatio, 
                                             regexpr("/", GeneRatio) + 1))) %>%
  mutate(percPath = Count/GeneRatTotal) %>%
  mutate(Enrichment = percPath/bkgdPerc)

DE_yMdrug_paths %>%
  filter(qvalue < 0.007) %>% 
  ggplot(aes(x = Enrichment, y = Description, color = p.adjust, size = Count)) +
  geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                             y = "KEGG pathway",
                                             color = "FDR", size = "Count") +
  theme_bw() + scale_color_gradient(low = "#f7941d", high = "#fed09e") + 
  ggtitle("Males C vs AR")

ggsave(filename = "yMpaths.pdf", height = 4, width = 6)

#females
DE_yFdrug <- tidy_yFdrug.dge %>% 
  filter(FDR < 0.05)

DE_yFdrug_entrez <- mapIds(org.Mm.eg.db, keys = DE_yFdrug$Ensembl,
                           keytype = "ENSEMBL", column = "ENTREZID")
DE_yFdrug_entrez <- na.exclude(DE_yFdrug_entrez)

DE_yFdrug_paths <- enrichKEGG(gene = DE_yFdrug_entrez, organism = "mmu",
                              pAdjustMethod = "BH", minGSSize = 4,
                              maxGSSize = 500, qvalueCutoff = 0.05)
write.csv(DE_yFdrug_paths, "YfemalePaths.csv")
#p < 0.0006
DE_yFdrug_paths <- as.data.frame(DE_yFdrug_paths) %>% 
  mutate(bkgdSize = as.numeric(substring(BgRatio, 
                                         regexpr("/", BgRatio) + 1))) %>%
  mutate(pathBkgd = as.numeric(substring(BgRatio, 1, 
                                         regexpr("/", BgRatio)-1))) %>%
  mutate(bkgdPerc = pathBkgd/bkgdSize) %>%
  mutate(GeneRatTotal = as.numeric(substring(GeneRatio, 
                                             regexpr("/", GeneRatio) + 1))) %>%
  mutate(percPath = Count/GeneRatTotal) %>%
  mutate(Enrichment = percPath/bkgdPerc)

DE_yFdrug_paths %>% 
  filter(pvalue < 0.0006) %>% 
  ggplot(aes(x = Enrichment, y = Description, color = p.adjust, size = Count)) +
  geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                             y = "KEGG pathway",
                                             color = "FDR", size = "Count") +
  theme_bw() + scale_color_gradient(low = "#d0417e", high = "#e7b0c1") + 
  ggtitle("Females C vs AR")

ggsave(filename = "yFpaths.pdf", height = 4, width = 5.55)







#pathfindR pathway analysis

genes_kegg <- mmu_kegg_genes
descriptions_kegg <- mmu_kegg_descriptions

##acquire STRING PIN file
url <- "https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz"
path2file <- file.path(tempdir(check = T), "STRING.txt.gz")
download.file(url, path2file)

mmu_string_df <- read.table(path2file, header = T)

mmu_string_df <- mmu_string_df[mmu_string_df$combined_score >= 800, ]

mmu_string_pin <- data.frame(Interactor_A = sub("^10090\\.","", mmu_string_df$protein1),
                             Interactor_B = sub("^10090\\.","", mmu_string_df$protein2))

library(biomaRt)

mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

converted <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol"),
                   filters = "ensembl_peptide_id",
                   values = unique(unlist(mmu_string_pin)),
                   mart = mmu_ensembl)

mmu_string_pin$Interactor_A <- converted$mgi_symbol[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
mmu_string_pin$Interactor_B <- converted$mgi_symbol[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]

mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

self_intr_cond <- mmu_string_pin$Interactor_A == mmu_string_pin$Interactor_B
mmu_string_pin <- mmu_string_pin[!self_intr_cond, ]

#remove duplicated interactions (including symmetric ones)
mmu_string_pin <- unique(t(apply(mmu_string_pin, 1, sort))) #returns matrix object

mmu_string_pin <- data.frame(A = mmu_string_pin[, 1],
                             pp = "pp",
                             B = mmu_string_pin[, 2])

path2SIF <- file.path(tempdir(), "mmusculusPIN.sif")
write.table(mmu_string_pin,
            file = path2SIF,
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)
path2SIF <- normalizePath(path2SIF)

#young Male drug treatment

tidy_yMdrug_subset <- tidy_yMdrug.dge[,c(2:3,5)] %>% 
  filter(FDR < 0.05)

write_csv(tidy_yMdrug_subset, "yM_sigTrans.csv")

tidy_yMdrug_subset <- data.frame(tidy_yMdrug_subset)

output_yMdrug <- run_pathfindR(input = tidy_yMdrug_subset,
                               convert2alias = F,
                               gene_sets = "Custom",
                               p_val_threshold = 0.05,
                               custom_genes = genes_kegg,
                               min_gset_size = 5,
                               max_gset_size = 500,
                               output_dir = "yMdrug",
                               adj_method = "fdr",
                               custom_descriptions = descriptions_kegg,
                               pin_name_path = path2SIF)

clustered_yMdrug <- cluster_enriched_terms(output_yMdrug, plot_dend = F,
                                           plot_clusters_graph = F)
knitr::kable(clustered_yMdrug[clustered_yMdrug$Status == "Representative", ])

write_csv(clustered_yMdrug, "yMdrug_PathClustered.csv")
top10_clustered_yMdrug <- subset(clustered_yMdrug, Cluster %in% 1:10)

top10enrich_yMdrug <- enrichment_chart(top10_clustered_yMdrug, plot_by_cluster = T)

ggsave("top10_yMdrug.pdf", top10enrich_yMdrug, height = 10, width = 8, useDingbats = F)

#young female drug treatment

tidy_yFdrug_subset <- tidy_yFdrug.dge[,c(2:3,5)] %>% 
  filter(FDR < 0.05)

write_csv(tidy_yFdrug_subset, "yF_sigTrans.csv")

tidy_yFdrug_subset <- data.frame(tidy_yFdrug_subset)

output_yFdrug <- run_pathfindR(input = tidy_yFdrug_subset,
                               convert2alias = F,
                               gene_sets = "Custom",
                               p_val_threshold = 0.05,
                               custom_genes = genes_kegg,
                               min_gset_size = 5,
                               max_gset_size = 500,
                               output_dir = "yMdrug",
                               adj_method = "fdr",
                               custom_descriptions = descriptions_kegg,
                               pin_name_path = path2SIF)

clustered_yFdrug <- cluster_enriched_terms(output_yFdrug, plot_dend = F,
                                           plot_clusters_graph = F)
knitr::kable(clustered_yFdrug[clustered_yFdrug$Status == "Representative", ])

write_csv(clustered_yFdrug, "yFdrug_PathClustered.csv")
top10_clustered_yFdrug <- subset(clustered_yFdrug, Cluster %in% 1:10)

top10enrich_yFdrug <- enrichment_chart(top10_clustered_yFdrug, plot_by_cluster = T)

ggsave("top10_yFdrug.pdf", top10enrich_yFdrug, height = 10, width = 8, useDingbats = F)

#middle male drug treatment
tidy_mMdrug_subset <- tidy_mMdrug.dge[,c(2:3,5)] %>% 
  filter(FDR < 0.05)

write_csv(tidy_mMdrug_subset, "mM_sigTrans.csv")

tidy_mMdrug_subset <- data.frame(tidy_mMdrug_subset)

output_mMdrug <- run_pathfindR(input = tidy_mMdrug_subset,
                               convert2alias = F,
                               gene_sets = "Custom",
                               p_val_threshold = 0.05,
                               custom_genes = genes_kegg,
                               min_gset_size = 5,
                               max_gset_size = 500,
                               output_dir = "mMdrug",
                               adj_method = "fdr",
                               custom_descriptions = descriptions_kegg,
                               pin_name_path = path2SIF)

clustered_mMdrug <- cluster_enriched_terms(output_mMdrug, plot_dend = F,
                                           plot_clusters_graph = F)
knitr::kable(clustered_mMdrug[clustered_mMdrug$Status == "Representative", ])

write_csv(clustered_mMdrug, "mMdrug_PathClustered.csv")
top10_clustered_mMdrug <- subset(clustered_mMdrug, Cluster %in% 1:10)

top10enrich_mMdrug <- enrichment_chart(top10_clustered_mMdrug, plot_by_cluster = T)

ggsave("top10_mMdrug.pdf", top10enrich_mMdrug, height = 10, width = 8, useDingbats = F)

#Middle Female drug treatment
tidy_mFdrug_subset <- tidy_mFdrug.dge[,c(2:3,5)] %>% 
  filter(FDR < 0.05)

write_csv(tidy_mFdrug_subset, "mF_sigTrans.csv")

tidy_mFdrug_subset <- data.frame(tidy_mFdrug_subset)

output_mFdrug <- run_pathfindR(input = tidy_mFdrug_subset,
                               convert2alias = F,
                               gene_sets = "Custom",
                               p_val_threshold = 0.05,
                               custom_genes = genes_kegg,
                               min_gset_size = 5,
                               max_gset_size = 500,
                               output_dir = "mFdrug",
                               adj_method = "fdr",
                               custom_descriptions = descriptions_kegg,
                               pin_name_path = path2SIF)

clustered_mFdrug <- cluster_enriched_terms(output_mFdrug, plot_dend = F,
                                           plot_clusters_graph = F)
knitr::kable(clustered_mFdrug[clustered_mFdrug$Status == "Representative", ])

write_csv(clustered_mFdrug, "mFdrug_PathClustered.csv")
top10_clustered_mFdrug <- subset(clustered_mFdrug, Cluster %in% 1:10)

top10enrich_mFdrug <- enrichment_chart(top10_clustered_mFdrug, plot_by_cluster = T)

ggsave("top10_mFdrug.pdf", top10enrich_mFdrug, height = 10, width = 8, useDingbats = F)

#Old Male drug treatment
tidy_oMdrug_subset <- tidy_oMdrug.dge[,c(2:3,5)] %>% 
  filter(FDR < 0.05)

write_csv(tidy_oMdrug_subset, "oM_sigTrans.csv")

tidy_oMdrug_subset <- data.frame(tidy_oMdrug_subset)

output_oMdrug <- run_pathfindR(input = tidy_oMdrug_subset,
                               convert2alias = F,
                               gene_sets = "Custom",
                               p_val_threshold = 0.05,
                               custom_genes = genes_kegg,
                               min_gset_size = 5,
                               max_gset_size = 500,
                               output_dir = "oMdrug",
                               adj_method = "fdr",
                               custom_descriptions = descriptions_kegg,
                               pin_name_path = path2SIF)

clustered_oMdrug <- cluster_enriched_terms(output_oMdrug, plot_dend = F,
                                           plot_clusters_graph = F)
knitr::kable(clustered_oMdrug[clustered_oMdrug$Status == "Representative", ])

write_csv(clustered_oMdrug, "oMdrug_PathClustered.csv")
top10_clustered_oMdrug <- subset(clustered_oMdrug, Cluster %in% 1:10)

top10enrich_oMdrug <- enrichment_chart(top10_clustered_oMdrug, plot_by_cluster = T)

ggsave("top10_oMdrug.pdf", top10enrich_oMdrug, height = 10, width = 8, useDingbats = F)

#Old female drug treatment
tidy_oFdrug_subset <- tidy_oFdrug.dge[,c(2:3,5)] %>% 
  filter(FDR < 0.05)

write_csv(tidy_oFdrug_subset, "oF_sigTrans.csv")

tidy_oFdrug_subset <- data.frame(tidy_oFdrug_subset)

output_oFdrug <- run_pathfindR(input = tidy_oFdrug_subset,
                               convert2alias = F,
                               gene_sets = "Custom",
                               p_val_threshold = 0.05,
                               custom_genes = genes_kegg,
                               min_gset_size = 5,
                               max_gset_size = 500,
                               output_dir = "oFdrug",
                               adj_method = "fdr",
                               custom_descriptions = descriptions_kegg,
                               pin_name_path = path2SIF)

clustered_oFdrug <- cluster_enriched_terms(output_oFdrug, plot_dend = F,
                                           plot_clusters_graph = F)
knitr::kable(clustered_oFdrug[clustered_oFdrug$Status == "Representative", ])

write_csv(clustered_oFdrug, "oFdrug_PathClustered.csv")
top10_clustered_oFdrug <- subset(clustered_oFdrug, Cluster %in% 1:10)

top10enrich_oFdrug <- enrichment_chart(top10_clustered_oFdrug, plot_by_cluster = T)

ggsave("top10_oFdrug.pdf", top10enrich_oFdrug, height = 10, width = 8, useDingbats = F)



### WGCNA analysis
library(gplots)

## data format

filtTPM <- allTpm %>% 
  pivot_longer(2:37, names_to = "sample", values_to = "tpm") %>% 
  extract(sample, into = c("animalID", "age", "diet", "sex"), regex = "([[:alnum:]]{2,})-([[:alnum:]]{4})-([[:upper:]]{1,})-([[:upper:]]{1})") %>% 
  filter(tpm > 0.1) %>% 
  group_by(Ensembl, age, diet, sex) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(n > 2) %>% 
  group_by(Ensembl, age, diet, sex) %>% 
  tally() %>% 
  tally(n > 2) %>% 
  tally(n > 1) %>% 
  tally(n > 1) %>% 
  filter(n > 2)

filtData_tpm <- allTpm %>% 
  filter(Ensembl %in% filtTPM$Ensembl)

wideTPM <- as.data.frame(t(filtData_tpm[,-1]))
names(wideTPM) <- filtData_tpm$Ensembl

gsg <- goodSamplesGenes(wideTPM, verbose = 3)
gsg$allOK

#cluster samples to look for outliers
sampleTree <- hclust(dist(wideTPM), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
pdf("clusterTree.pdf")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "",
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

disableWGCNAThreads()

sft <- pickSoftThreshold(wideTPM, powerVector = powers, verbose = 5)

pdf("NetworkTopology.pdf")
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")

abline(h = 0.9, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1,
     col = "red")
dev.off()

#thresholding power is 14.
cor <- WGCNA::cor

net <- blockwiseModules(wideTPM, power = 14, TOMType = "signed",
                        minModuleSize = 50, reassignThreshold = 0,
                        maxBlockSize = 14500, mergeCutHeight = 0.25,
                        numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                        saveTOMFileBase = "weasMus", verbose = 3)

moduleSizes <- data.frame(table(net$colors))
colnames(moduleSizes) <- c("ME#", "Freq")
MEColors <- data.frame(table(labels2colors(net$colors)))

moduleStats <- merge (moduleSizes, MEColors)
colnames(moduleStats) <- c("# genes", "ME#", "colors")

write_csv(moduleStats, "moduleSizes.csv")

#dendrogram
mergedColors <- labels2colors(net$colors)

pdf("moduleDendrogram.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module Colors", dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)
dev.off()

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)

MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

MEs0 <- moduleEigengenes(wideTPM, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

write_csv(MEs0, "moduleNames.csv")

pdf("modRelate.pdf", width = 6, height = 10)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,2,0),
                      marHeatmap = c(3,4,2,2), cex.lab = 0.8)
dev.off()

#pathway analysis for modules
modGreenYellowEntrez <- mapIds(org.Mm.eg.db,
                               keys = names(wideTPM)[moduleColors == "greenyellow"],
                               keytype = "ENSEMBL", column = "ENTREZID")
modGreenYellowEntrez <- na.exclude(modGreenYellowEntrez)

modGreenYellowPaths <- enrichKEGG(gene = modGreenYellowEntrez, organism = "mmu",
                                  pAdjustMethod = "BH", minGSSize = 4,
                                  maxGSSize = 500, qvalueCutoff = 0.05)

write.csv(modGreenYellowPaths, "modGreenYellowPaths.csv")

#visualization
modGreenYellowPaths <- as.data.frame(modGreenYellowPaths) %>% 
  mutate(bkgdSize = as.numeric(substring(BgRatio,
                                         regexpr("/", BgRatio) + 1))) %>% 
  mutate(pathBkgd = as.numeric(substring(BgRatio, 1,
                                         regexpr("/", BgRatio) - 1))) %>% 
  mutate(bkgdPerc = pathBkgd/bkgdSize) %>% 
  mutate(GeneRatTotal = as.numeric(substring(GeneRatio,
                                             regexpr("/", GeneRatio) + 1))) %>% 
  mutate(percPath = Count/GeneRatTotal) %>% 
  mutate(Enrichment = percPath/bkgdPerc)

modGreenYellowPaths %>% 
  filter(Enrichment > 10) %>% 
  ggplot(aes(x = Enrichment, y = Description, color = p.adjust, size = Count)) +
  geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                             y = "KEGG pathway",
                                             color = "FDR", size = "Count") + 
  theme_bw() + scale_color_gradient(low = "#E0E721", high = "#F0F1B1") +
  ggtitle("Green-Yellow Mod Paths")
ggsave(filename = "modGreenYellowPaths.pdf")


#Purple
modPurpleEntrez <- mapIds(org.Mm.eg.db,
                          keys = names(wideTPM)[moduleColors == "purple"],
                          keytype = "ENSEMBL", column = "ENTREZID")
modPurpleEntrez <- na.exclude(modPurpleEntrez)

modPurplePaths <- enrichKEGG(gene = modPurpleEntrez, organism = "mmu",
                             pAdjustMethod = "BH", minGSSize = 4,
                             maxGSSize = 500, qvalueCutoff = 0.05)

write.csv(modPurplePaths, "modPurplePaths.csv")
#No pathways enriched

#Magenta
modMagentaEntrez <- mapIds(org.Mm.eg.db,
                           keys = names(wideTPM)[moduleColors == "magenta"],
                           keytype = "ENSEMBL", column = "ENTREZID")
modMagentaEntrez <- na.exclude(modMagentaEntrez)

modMagentaPaths <- enrichKEGG(gene = modMagentaEntrez, organism = "mmu",
                             pAdjustMethod = "BH", minGSSize = 4,
                             maxGSSize = 500, qvalueCutoff = 0.05)

write.csv(modMagentaPaths, "modMagentaPaths.csv")

modMagentaPaths <- as.data.frame(modMagentaPaths) %>% 
  mutate(bkgdSize = as.numeric(substring(BgRatio,
                                         regexpr("/", BgRatio) + 1))) %>% 
  mutate(pathBkgd = as.numeric(substring(BgRatio, 1,
                                         regexpr("/", BgRatio) - 1))) %>% 
  mutate(bkgdPerc = pathBkgd/bkgdSize) %>% 
  mutate(GeneRatTotal = as.numeric(substring(GeneRatio,
                                             regexpr("/", GeneRatio) + 1))) %>% 
  mutate(percPath = Count/GeneRatTotal) %>% 
  mutate(Enrichment = percPath/bkgdPerc)

modMagentaPaths %>% 
  filter(pvalue < 0.0018) %>% 
  ggplot(aes(x = Enrichment, y = Description, color = p.adjust, size = Count)) +
  geom_point() + expand_limits(x = 0) + labs(x = "Enrichment",
                                             y = "KEGG pathway",
                                             color = "FDR", size = "Count") + 
  theme_bw() + scale_color_gradient(low = "#d0417e", high = "#e7b0c1") +
  ggtitle("Magenta Mod Paths")
ggsave(filename = "modMagentaPaths.pdf", height = 4, width = 6.5)

normTPM <- filtData_tpm %>% 
  pivot_longer(2:37, names_to = "sample", values_to = "tpm") %>% 
  extract(sample, into = c("animalID", "age", "diet", "sex"), 
          regex = "([[:alnum:]]{2,})-([[:alnum:]]{4})-([[:upper:]]{1,})-([[:upper:]]{1})") %>%
  mutate(logTpm  = log2(tpm)) %>% 
  group_by(animalID) %>% 
  mutate(normTpm = scale(logTpm, center = T, scale = T)[,1]) 



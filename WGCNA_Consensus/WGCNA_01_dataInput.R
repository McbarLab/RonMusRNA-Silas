trimList <- 1:37
AR_list <-
  c(4:6, 10:12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 31, 33, 36)
maleList <- c(1:7, 14:19, 26:31)

# 1.a Loading expression data
RawTPM = read.csv("rsubread_GENE_tpm.csv")

# Trim the data by sex
maleData <- RawTPM[,maleList]
femList <- trimList[!trimList %in% maleList]
femData <- RawTPM[,c(1,femList)]

# We work with two sets:
nSets = 2

# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female", "Male")
shortLabels = c("Female", "Male")

# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(femData[-c(1)])))
names(multiExpr[[1]]$data) = femData$gene
rownames(multiExpr[[1]]$data) = names(femData)[-c(1)]
multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1)])))
names(multiExpr[[2]]$data) = maleData$gene
rownames(multiExpr[[2]]$data) = names(maleData)[-c(1)]

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# Print whether the vector is valid
if (exprSize$structureOK == TRUE) {
  print("Multiple-set Vector Successfully Created!")
} else{
  stop("Multiple-set Vector Creation Failed!")
}

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
    printFlush(paste("In set", setLabels[set], "removing samples",
                     paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

# Cluster the 2 sets to see if there are any outliers
pdf(file = "SampleClustering.pdf", width = 12, height = 12)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)
dev.off()

# 1.c Loading clinical trait data
traitData <- read.csv("Weasley_Biometrics_ALL.csv")
allTraits = traitData[,-c(2,5:22)] #this removes these specific columns

Traits = vector(mode="list", length = nSets)

for(set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data)
  traitRows = match(setSamples, allTraits$Animal.ID)
  Traits[[set]] = list(data = allTraits[traitRows, -1])
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1]
}
collectGarbage()
# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples

save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "Consensus-dataInput.RData")
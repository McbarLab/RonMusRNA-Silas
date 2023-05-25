# 1.a Loading expresssion data
RawTPM = read.csv("rsubread_GENE_tpm.csv", check.names = FALSE)
datExpr = as.data.frame(t(RawTPM[,-c(1)]))
names(datExpr) = RawTPM$gene
rownames(datExpr) = names(RawTPM)[-c(1)]

# 1.b Checking data for excessive missing values and identification of outlier microarray samples
## Make sure that the numbers to not contain commas in the thousands position (e.g.- 1,234=NO; 1234=YES)
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

## Remove offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
  
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg$allOK

## Cluster the samples
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# 1.c Loading clinical trait data
traitData = read.csv("Weasley_Biometrics_07March2023.csv")

allTraits = traitData[,-c(2:32)] #this removes these specific columns

ALLsamples = rownames(datExpr)
traitRows = match(ALLsamples, allTraits$Animal.ID)
datTraits = allTraits[traitRows,-1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

## Visualize how the clinical traits relate to the sample dendrogram
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2,
                    traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "01-dataInput.RData")
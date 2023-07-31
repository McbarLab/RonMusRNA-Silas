# mode <- "Male"
mode <- "ALL"
# mode <- "Female"

trimList <- 1:36
AR_list <-
  c(4:6, 10:12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 31, 33, 36)
maleList <- c(1:6, 13:18, 25:30)
if (mode == "ALL") {
  trimList <- trimList
} else if (mode == "Male") {
  trimList <- maleList
} else if (mode == "Female") {
  trimList <- trimList[!trimList %in% maleList]
} else{
  stop("The mode can only be one of the following: ALL, Male, Female")
}

# 1.c Loading clinical trait data
traitData = read.csv("Weasley_Biometrics_ALL.csv")
allTraits = traitData[, -c(5:22)] #this removes these specific columns

# Trim the data by sex
allTraits <- allTraits[trimList, ]

# 1.a Loading expression data
RawTPM = read.csv("rsubread_GENE_tpm.csv", check.names = FALSE)
datExpr = as.data.frame(t(RawTPM[, -c(1)]))
names(datExpr) = RawTPM$gene
rownames(datExpr) = names(RawTPM)[-c(1)]

# Trim dataExpr to match dimension of traits
datExpr <- datExpr[trimList, ]

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


# Split age/group into pairwise comparisons
allTraits <-
  binarizeCategoricalColumns(
    data = allTraits,
    convertColumns = c("Sex", "Age.Group", "Diet"),
    includePairwise = TRUE,
    includeLevelVsAll = FALSE
  )

ALLsamples = rownames(datExpr)
traitRows = match(ALLsamples, allTraits$Animal.ID)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

## Visualize how the clinical traits relate to the sample dendrogram
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2,
                    traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(mode, AR_list, datExpr, datTraits, file = "01-dataInput.RData")
# Define variable Sex.M.vs.F containing the Sex.M.vs.F column of datTrait
Sex.M.vs.F = as.data.frame(datTraits$Sex.M.vs.F);
names(Sex.M.vs.F) = "Sex.M.vs.F"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Sex.M.vs.F, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Sex.M.vs.F), sep="");
names(GSPvalue) = paste("p.GS.", names(Sex.M.vs.F), sep="");

module = "brown"
column = match(module, modNames)

moduleGenes = moduleColors == module

sizeGrWindow(7, 7)

par(mfrow = c(1, 1))

verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for sex",
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = module
)
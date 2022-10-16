#install.packages("BiocManager")
#BiocManager::install("WGCNA")
#BiocManager::install("tximport")
library(WGCNA)
library(readr)
library(tximport)
library(tximportData)
library(DESeq2)
library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(apeglm)
library(pcaExplorer)
library(dplyr)
#install.packages('tximportData')

options(stringsAsFactors = FALSE)
enableWGCNAThreads()


#### READ IN GENE COUNTS TABLE FROM SALMON ########
sample = read.table('metadata_salmon.txt', header = TRUE)
sample = unique(sample)
sample_files = 
  paste0(as.vector(sample$sample), '/quant.sf')
names(sample_files) = 
  as.vector(sample$sample)
gene_map = read.table('gene_map.txt') 
colnames(gene_map) <- c('transid', 'geneid')
map_gene = data_frame(gene_map) 

count_data = tximport(files = sample_files, 
                      type = 'salmon',
                      tx2gene = map_gene)

cts = (count_data$counts)
dim(cts)
metadata = sample


#### UPLOAD TRAIT DATA ####
traitData = metadata
dim(traitData)
names(traitData)
IctiSamples = rownames(datExpr0)
IctiSamples
traitRows = match(IctiSamples, traitData$sample)
datTraits = as.data.frame(cbind(traitData$metabolic_state, 
                                traitData$sex), 
                          stringsAsFactors = TRUE)# keep only sex and IBA or SA
colnames(datTraits) = c('metabolic_state', 'sex')
datTraits
rownames(datTraits) = traitData[traitRows, 3] # Specify sample names
rownames(datTraits)


### CHECK FOR MISSING VALUES FROM COUNTS DATA ####
datExpr0 = as.data.frame(t(cts))
names(datExpr0)
rownames(datExpr0)
#check genes for for missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
#since false...
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


### CLUSTER SAMPLES TO LOOK FOR OUTLIERS ####
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(20,10)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 80000, col = "red");



#### CLUSTER SAMPLES WITH TRAIT DATA ####
sampletree1 = hclust(dist(datExpr0), method = 'average')
#Convert traits to color
traitColors = labels2colors(datTraits)
plotDendroAndColors(sampletree1, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


#### FIND AND CHOOSING THRESHOLD POWER ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2))
cex1=1.0
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



#### NETWORK AND MODULE CONSTRUCTION ####
net = blockwiseModules(datExpr0, power = 14,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "IctiTOM",
                       verbose = 3)
table(net$blocks)
# open a graphics window
par(mfrow = c(3,3))
cex1=1.0
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath for each block
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for each block
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for each block
plotDendroAndColors(net$dendrograms[[3]], mergedColors[net$blockGenes[[3]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for each block
plotDendroAndColors(net$dendrograms[[4]], mergedColors[net$blockGenes[[4]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for each block
plotDendroAndColors(net$dendrograms[[5]], mergedColors[net$blockGenes[[5]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for each block
plotDendroAndColors(net$dendrograms[[6]], mergedColors[net$blockGenes[[6]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#Save for subsequent analysis
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction-auto.RData")



##### FIND AND QUANTIFY MODULE-TRAIT ASSOCIATIONS ####

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEsnew = sapply(MEs, unclass)

datTraitsnew = sapply(datTraits, unclass) #Turn categorical variables to numeric
moduleTraitCor = cor(MEs, datTraitsnew, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(15,15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 1, 4))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = greenWhiteRed(12),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-.25,.25),
               main = paste("Module-trait relationships"))



#### GENE SIGNIFICANCE AND MODULE MEMBERSHIP ####
# Define variable weight containing the weight column of datTrait
metState = as.data.frame(datTraits$metabolic_state)
metStatenew = sapply(metState, unclass) # characters to numeric
names(metStatenew) = "metState"

sex = as.data.frame(datTraits$sex)
sexNew = sapply(sex, unclass)
names(sexNew) = "sex"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#Gene sig for metState
geneTraitSignificance = as.data.frame(cor(datExpr0, metStatenew, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(metState), sep="");
names(GSPvalue) = paste("p.GS.", names(metState), sep="");

#Genesig for sex
geneTraitSignificance1 = as.data.frame(cor(datExpr0, sexNew, use = "p"));
GSPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance1), nSamples));
names(geneTraitSignificance1) = paste("GS.", names(sex), sep="");
names(GSPvalue1) = paste("p.GS.", names(sex), sep="");


#### INTRAMODULAR ANALYSIS: ID GENES W/ HIGH GENE SIG AND MODULE MEMBERSHIP ####
module1 = "red" # Shows high correlation with metState
column = match(module1, modNames);
moduleGenes = moduleColors==module1;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module1, "module"),
                   ylab = "Gene significance for metabolic state",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module1)

module2 = "darkslateblue" # Shows high correlation with sex
column = match(module2, modNames);
moduleGenes = moduleColors==module2;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   
                   abs(geneTraitSignificance1[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module2, "module"),
                   ylab = "Gene significance for metabolic state",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module2)


module3 = "salmon4" # Shows high correlation with sex
column = match(module3, modNames);
moduleGenes = moduleColors==module3;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance1[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module3, "module"),
                   ylab = "Gene significance for metabolic state",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module3)


module4 = "paleturquoise" # Shows high correlation with sex
column = match(module4, modNames);
moduleGenes = moduleColors==module4;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance1[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module4, "module"),
                   ylab = "Gene significance for metabolic state",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module4)



#### SUMMARY OUTPUT OF NETWORK ANALYSIS RESULTS #### 

names(datExpr0)
write.csv((names(datExpr0) [moduleColors=='red']), file = 'redGene.csv')
write.csv((names(datExpr0) [moduleColors=='darkslateblue']), file = 'darkslateblueGene.csv')
write.csv(as.data.frame(names(datExpr0) [moduleColors=='salmon4']), file = 'salmon4Gene.csv')
write.csv(as.data.frame(names(datExpr0) [moduleColors=='paleturquoise']), file = 'paleturquoise.csv')





#### VISULIZING THE NETWORK OF EIGENGENES for metstate and sex ####
# Recalculate module eigengenes
MEs1 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
# Isolate metstate from the clinical traits
metState1 = as.data.frame(datTraits$metabolic_state)
metState1new = sapply(metState1, unclass)
names(metState1) = "metStatet"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs1, metState1new))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Recalculate module eigengenes
MEs2 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
# Isolate sex from the clinical traits
metState2 = as.data.frame(datTraits$sex)
metState2new = sapply(metState2, unclass)
names(metState2) = "Sex"
# Add the weight to existing module eigengenes
MET1 = orderMEs(cbind(MEs2, metState2new))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET1, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


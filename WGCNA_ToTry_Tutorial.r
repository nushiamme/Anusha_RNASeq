#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
library(here)
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
#setwd(workingDir); 
# Load the package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set
femData = read.csv(here("WGCNATutorial", "FemaleLiver-Data", "LiverFemale3600.csv"))
# Read in the male liver data set
# maleData = read.csv(here("WGCNATutorial", "MaleLiver-Data", "LiverMale3600.csv"))
# Take a quick look at what is in the data sets (caution, longish output):
dim(femData)
names(femData)
# dim(maleData)
# names(maleData)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# # We work with two sets:
# nSets = 2
# # For easier labeling of plots, create a vector holding descriptive names of the two sets.
# setLabels = c("Female liver", "Male liver")
# shortLabels = c("Female", "Male")
# # Form multi-set expression data: columns starting from 9 contain actual expression data.
# multiExpr = vector(mode = "list", length = nSets)
# 
# multiExpr[[1]] = list(data = as.data.frame(t(femData[-c(1:8)])));
# names(multiExpr[[1]]$data) = femData$substanceBXH;
# rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];
# multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1:8)])));
# names(multiExpr[[2]]$data) = maleData$substanceBXH;
# rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)];
# # Check that the data has the correct format for many functions operating on multiple sets:
# exprSize = checkSets(multiExpr)
# 
# 
## Alternative
datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH;
rownames(datExpr0) = names(femData)[-c(1:8)];



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# # Check that all genes and samples have sufficiently low numbers of missing values.
# gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
# gsg$allOK


## Alternative
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# if (!gsg$allOK)
# {
#   # Print information about the removed genes:
#   if (sum(!gsg$goodGenes) > 0)
#     printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
#                                               collapse = ", ")))
#   for (set in 1:exprSize$nSets)
#   {
#     if (sum(!gsg$goodSamples[[set]]))
#       printFlush(paste("In set", setLabels[set], "removing samples",
#                        paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
#     # Remove the offending genes and samples
#     multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
#   }
#   # Update exprSize
#   exprSize = checkSets(multiExpr)
# }


## Alternative
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# sampleTrees = list()
# for (set in 1:nSets)
# {
#   sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
# }


## Alternative
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# pdf(file = here::here("WGCNATutorial", "SampleClustering.pdf"), width = 12, height = 12);
# par(mfrow=c(2,1))
# par(mar = c(0, 4, 2, 0))
# for (set in 1:nSets)
#   plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
#        xlab="", sub="", cex = 0.7);
# dev.off();

## Alternative
# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Choose the "base" cut height for the female data set
# baseHeight = 16
# # Adjust the cut height for the male data set for the number of samples
# cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# # Re-plot the dendrograms including the cut lines
# pdf(file = here::here("WGCNATutorial", "SampleClustering.pdf"), width = 12, height = 12);
# par(mfrow=c(2,1))
# par(mar = c(0, 4, 2, 0))
# for (set in 1:nSets)
# {
#   plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
#        xlab="", sub="", cex = 0.7);
#   abline(h=cutHeights[set], col = "red");
# }
# dev.off();


## Ignore
#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# for (set in 1:nSets)
# {
#   # Find clusters cut by the line
#   labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
#   # Keep the largest one (labeled by the number 1)
#   keep = (labels==1)
#   multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
# }
# collectGarbage();
# # Check the size of the leftover data
# exprSize = checkSets(multiExpr)
# exprSize

##Ignore

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# traitData = read.csv(here::here("WGCNATutorial", "FemaleLiver-Data", "ClinicalTraits.csv"))
# dim(traitData)
# names(traitData)
# # remove columns that hold information we do not need.
# allTraits = traitData[, -c(31, 16)];
# allTraits = allTraits[, c(2, 11:36) ];
# # See how big the traits are and what are the trait and sample names
# dim(allTraits)
# names(allTraits)
# allTraits$Mice
# # Form a multi-set structure that will hold the clinical traits.
# Traits = vector(mode="list", length = nSets);
# for (set in 1:nSets)
# {
#   setSamples = rownames(multiExpr[[set]]$data);
#   traitRows = match(setSamples, allTraits$Mice);
#   Traits[[set]] = list(data = allTraits[traitRows, -1]);
#   rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
# }
# collectGarbage();
# # Define data set dimensions
# nGenes = exprSize$nGenes;
# nSamples = exprSize$nSamples;


##Alternative
traitData = read.csv(here("WGCNATutorial", "FemaleLiver-Data", "ClinicalTraits.csv"));
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Mice);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();

#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

# Re-cluster samples
# sampleTree2 = hclust(dist(datExpr), method = "average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datTraits, signed = FALSE);
# # Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits),
#                     main = "Sample dendrogram and trait heatmap")
# save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, 
#      file = here::here("WGCNATutorial", "Consensus-dataInput.RData"))

## Alternative
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = here("WGCNATutorial", "FemaleLiver-01-dataInput.RData"))


#=====================================================================================
#
#  PART 2
#
#=====================================================================================



#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
#setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = here("WGCNATutorial", "FemaleLiver-01-dataInput.RData"));
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
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


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = here("WGCNATutorial", "femaleMouseTOM"), 
                       verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = here("WGCNATutorial", "FemaleLiver-02-networkConstruction-auto.RData"))


#=====================================================================================
#=====================================================================================
#
#  Part 3
#
#=====================================================================================
#=====================================================================================


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = ".";
# setwd(workingDir); 
# Load the WGCNA package
# library(WGCNA) #run
# The following setting is important, do not omit.
# options(stringsAsFactors = FALSE); #run
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
# enableWGCNAThreads() #run
# Load the data saved in the first part
# lnames = load(file = "FemaleLiver-01-dataInput.RData"); #run
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
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


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


softPower = 6;
adjacency = adjacency(datExpr, power = softPower);


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = here("WGCNATutorial", "FemaleLiver-02-networkConstruction-stepByStep.RData"))


rm(list = ls())

#=====================================================================================
#=====================================================================================
#
#  Part 4
#
#=====================================================================================
#=====================================================================================

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = ".";
# setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
library(here)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
# enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = here("WGCNATutorial", "FemaleLiver-01-dataInput.RData"));
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
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


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 6, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = here("WGCNATutorial", "femaleMouseTOM-blockwise"),
                         verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Load the results of single-block analysis
load(file = here("WGCNATutorial", "FemaleLiver-02-networkConstruction-auto.RData"));
# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


sizeGrWindow(12,9)
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes;
blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes;


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)

#=====================================================================================
#=====================================================================================
#
#  Part 5 
#
#=====================================================================================
#=====================================================================================


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
# getwd();
# # If necessary, change the path below to the directory where the data files are stored. 
# # "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = ".";
# setwd(workingDir); 
# # Load the WGCNA package
# library(WGCNA)
# # The following setting is important, do not omit.
# options(stringsAsFactors = FALSE);
# # Load the expression and trait data saved in the first part
lnames = load(file = here("WGCNATutorial", "FemaleLiver-01-dataInput.RData"));
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = here("WGCNATutorial", "FemaleLiver-02-networkConstruction-auto.RData"));
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


names(datExpr)[moduleColors=="brown"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


annot = read.csv(file = here("WGCNATutorial", "FemaleLiver-Data", "GeneAnnotation.csv"));
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


write.csv(geneInfo, file = here("WGCNATutorial", "geneInfo.csv"))

#=====================================================================================
#=====================================================================================
#
#  Part 6
#
#=====================================================================================
#=====================================================================================



#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
# getwd();
# # If necessary, change the path below to the directory where the data files are stored. 
# # "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = ".";
# setwd(workingDir); 
# # Load the WGCNA package
library(WGCNA)
library(here)
# # The following setting is important, do not omit.
# options(stringsAsFactors = FALSE);
#BiocManager::install("org.Mm.eg.db") # For versions of R before 4.2
library(org.Mm.eg.db)
# # Load the expression and trait data saved in the first part
lnames = load(file = here("WGCNATutorial", "FemaleLiver-01-dataInput.RData"));
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = here("WGCNATutorial", "FemaleLiver-02-networkConstruction-auto.RData"));
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Read in the probe annotation
annot = read.csv(file = here("WGCNATutorial", "FemaleLiver-Data", "GeneAnnotation.csv"));
# Match probes in the data set to the probe IDs in the annotation file 
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
intModules = c("brown", "red", "salmon")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = here("WGCNATutorial", "LocusLinkIDs-all.txt");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


tab = GOenr$bestPTerms[[4]]$enrichment


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


names(tab)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


write.table(tab, file = here("WGCNATutorial", "GOEnrichmentTable.csv"), sep = ",", quote = TRUE, row.names = FALSE)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

#=====================================================================================
#=====================================================================================
#
#  Part 7 - visualization
#
#=====================================================================================
#=====================================================================================


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
# getwd();
# # If necessary, change the path below to the directory where the data files are stored. 
# # "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = ".";
# setwd(workingDir); 
# # Load the WGCNA package
# library(WGCNA)
# # The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = here("WGCNATutorial", "FemaleLiver-02-networkConstruction-auto.RData"));
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


#=====================================================================================
#
#  Code chunk 3 - optional - restricted number of genes to 400
#
#=====================================================================================


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


#=====================================================================================
#
#  Code chunk 4 - Visualizing the network of eigengenes
#
#=====================================================================================


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


#=====================================================================================
#=====================================================================================
#
#  Part 8 - exporting for Cytoscape etc.
#
#=====================================================================================
#=====================================================================================

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
# getwd();
# # If necessary, change the path below to the directory where the data files are stored. 
# # "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = ".";
# setwd(workingDir); 
# # Load the WGCNA package
# library(WGCNA)
# # The following setting is important, do not omit.
# options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = here("WGCNATutorial", "FemaleLiver-01-dataInput.RData"));
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = here("WGCNATutorial", "FemaleLiver-02-networkConstruction-auto.RData"));
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = here("WGCNATutorial", "GeneAnnotation.csv"));
# Select module
module = "brown";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = here("WGCNATutorial", paste("VisANTInput-", module, ".txt", sep="")),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = here("WGCNATutorial", paste("VisANTInput-", module, "-top30.txt", sep="")),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = here("WGCNATutorial", "FemaleLiver-Data", "GeneAnnotation.csv"));
# Select modules
modules = c("brown", "red");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = here("WGCNATutorial", paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep="")),
                               nodeFile = here("WGCNATutorial", paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep="")),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);




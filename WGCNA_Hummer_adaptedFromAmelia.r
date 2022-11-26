## Oct 31 - code from Amelia Demery, trying out on my liver data

# WCGNA - starlings - data preprocessing - BEAK ---------------------------------------
#setwd("C:/Users/Lukyr/OneDrive/Desktop/github/starlings/wgcna-tutorial")
#library(BiocManager)
library(WGCNA) # BiocManager::install("WGCNA")
library(here)
library(dplyr)
options(stringsAsFactors = F)

# clean the data set ----
## loading expression data
star <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                        "final", "RNASeq_NormCounts.csv"))
#star$gene_symbol <- star$gene
dim(star)
#rownames(star) <- star$gene
head(star)

#### Read in metadata
metadat <- read.table(here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                           "final", "AnushaShankar_Liver_RNASeq_metadata.csv"), sep = ",", header = T)
liversamples <- metadat$X #[metadat$Tissue=="Liver"] 

# remove auxiliary data and only keep the columns for expression
datExpr_pre <- star %>% 
  dplyr::select(any_of(liversamples), gene)
datExpr0 <- t(datExpr_pre[,-19])
colnames(datExpr0) <- datExpr_pre$gene
rownames(datExpr0) <- names(datExpr_pre[1:18])

# check data for excessive missing values and identification of outlier microarray samples
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # if true then we are all set

# if false, then we remove offending genes and samples from data
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

# cluster the samples (not the same as clustering genes) to look for any obvious outliers
sampleTree <- hclust(dist(datExpr0), method = "average")

# plot the sample tree
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# one outlier, T2.S24.BE. You have options to remove by hand OR use an automatic approach
# Plot a line to show the cut
# AS: using normalized data, not going to remove outliers
# abline(h = 1500000, col = "red")

# Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 1500000, minSize = 10)
# table(clust)
# 
# # clust 1 contains the samples we want to keep.
# keepSamples <-  (clust==1)
datExpr <-  datExpr0 #[keepSamples, ]
nGenes <-  ncol(datExpr)
nSamples <-  nrow(datExpr)

## load the clinical trait data
# read in data, match the samples for which they were measured to the expression samples
# traitData <- read.table("../data/samples_ALL_with-color.txt", header = T, sep = ",") %>% 
#   pivot_wider(names_from = trait,
#               values_from = mean.prop)

traitData <- metadat #%>% 
#   pivot_wider(names_from = trait,
#               values_from = mean.prop)
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits <-  traitData %>% 
  #filter(sample.b != "T2.S24.BE") %>% 
  select(X, BirdID_Tissue, Metabolic_State, CaptureWeight_g, Tmax_euthanasia, Temp_diff)

# Form a data frame analogous to expression data that will hold the clinical traits.
livSamples <-  rownames(datExpr)
traitRows <- match(livSamples, allTraits$X)
datTraits <-  allTraits #[traitRows, -1]
rownames(datTraits) <- datTraits$X
collectGarbage()

# Re-cluster samples
sampleTree2 <-  hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
# datTraits <- as.matrix(datTraits)
# names(datTraits) <- rownames(datTraits)
datTraits[7] <- seq(1:18)
traitColors <-  numbers2colors(as.numeric(as.matrix(datTraits[7])), signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                                     "final", "livWGCNA-data.RData"))


# WCGNA - starlings - network analysis - BEAK ------------------------------------
# setwd("C:/Users/Lukyr/OneDrive/Desktop/github/starlings/data/")
# library(BiocManager)
# library(WGCNA) # BiocManager::install("WGCNA")
# options(stringsAsFactors = F)

### CLEANING WITH THE FEMALE DATASET
# ## loading expression data
# star <- read.table("4043D_rawCounts.txt", sep = "\t", header = T)
# star$gene_symbol <- star$X
# dim(star)
# names(star)

lnames <- load(file = here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                           "final", "livWGCNA-data.RData"))

# choose a set of soft-thresholding poweers
powers <- c(1:20)

# call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

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
# choose the lowest power for which the scale-free topology fit index reaches 0.9. That would be 18.

## co-expression similarity and adjacency
softPower <- 18

# one-step network construction and module detection
datExpr <- datExpr*1.0 # got this from the web to turn integers into numerics...it worked!!!

net <-blockwiseModules(datExpr, power = 18,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                                              "final","LiverTOM"),
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                 "final","hummer-modules_liver.RData"))

# WGCNA - exporting a gene network to external visualization softw --------
# library(WGCNA)
# options(stringsAsFactors = F)
## load in the data
lnames <- load(file = here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                            "final", "livWGCNA-data.RData"))
datExpr <- datExpr*1.0 # got this from the web to turn integers into numerics...it worked!!!
### exporting to VisANT
# recalculate topological overlap
TOM <- TOMsimilarityFromExpr(datExpr, power = 18)
# read in annotation
annot <- read.csv(file = "FemaleLiver-Data/GeneAnnotation.csv")
head(annot)
module <- "grey"
# select module probe
probes <- names(datExpr)
inModule <- (moduleColors == module)
modProbes <- probes[inModule]
# select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)#,
                            #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )
# restrict the genes in output to top 30 hub genes in brown module
nTop = 30
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("brown", "red");
# Select module probes
probes = colnames(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
## modGenes = datExpr$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                                               "final",paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep="")),
                               nodeFile = here("DESeq_Data_Mar2022_all tissues", "all tissues", 
                                               "final",paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep="")),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

library(biomaRt)
library(biomartr)
biomartr::getMarts()
head(biomaRt::listMarts(host = "https://www.ensembl.org"), 10)

head(biomaRt::listDatasets(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")), 10)     

head(biomaRt::listAttributes(biomaRt::useDataset(
  dataset = "hsapiens_gene_ensembl",         
  mart    = useMart("ENSEMBL_MART_ENSEMBL",      
                    host    = "https://www.ensembl.org"))), 10)


#### DeSeq for GSEA analysis - brain RNA data
## Anusha Shankar, github: nushiamme
## Code started June 17, 2022

## Initially followed this tutorial:
## https://www.youtube.com/watch?v=OzNzO8qwwp0
## And then for visualizations
## https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md
## For GSEA: https://www.youtube.com/watch?v=KY6SS4vRchY

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('DESeq2')
# BiocManager::install('EnhancedVolcano') 
# BiocManager::install('apeglm') 
# Enhanced Volcano plot help: 
# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(here)
library(viridis)
library(gridExtra)
library(apeglm)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(dendextend) # For making gene trees for heatmaps
library(ComplexHeatmap)


## If you opened the .Rproj file from the OneDrive folder, reading in these files will work- 
## the relative paths should be the same
data <- read.csv(here("..//DESeq_Data_Brains_Jun2022//RNASeqBrain_rawcounts_data.csv"))
meta <- read.csv(here("..//DESeq_Data_Brains_Jun2022//AnushaShankar_brainRNASeq_metadata.csv"))
#star_alignment <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//star_alignment_plot.csv"))
#star <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//star_gene_counts.csv"))
#foldchange <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//TopFoldChanges.csv"))
## Made this one in this script
norm_counts_df <- read.csv(here("..//DESeq_Data_Brains_Jun2022//RNASeqBrain_NormCounts.csv"))

my_theme <- theme_classic(base_size = 15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.key.height = unit(2, "line"))

## Viridis colors
my_gradient <- c("#823de9", "#7855ce", "#6e6eb2", "#648697", "#599e7c", "#4fb760", "#45cf45")
my_col_rainbows <- c("#f94144", "#f3722c", "#f8961e", "#f9844a",
                     "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1")

#### Don't need to re-do - this is a check, not used in analyses really. #####                    
## Checking on the star alignment
head(star)
star <- star %>%
  rename(Sample = ï..Sample) %>%
  mutate(Sample = gsub("c", "", Sample))

#meta2 <- meta %>%
 # rename(Sample = X)

meta2 <- meta

starlong_meta <- star %>%
  gather(key="Mapped", value="count", -Sample) %>%
  left_join(., meta2, by="Sample")

starlong_meta %>%
  mutate(Mapped_relevel = 
           fct_relevel(Mapped, 
                       "Overlapping_Genes", "Multimapping", "Unmapped", "Ambiguous_Features", "No_Feature")) %>%
  ggplot(., aes(Tissue, count, fill=Mapped_relevel)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis(discrete = T) + my_theme + ylab("Percent counts")


starlong_meta %>%
  ggplot(., aes(Tissue, count, fill=Mapped)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis(discrete = T)

starlong_meta %>%
  ggplot(., aes(Tissue, count, fill=Mapped)) + 
  geom_bar(stat="identity", position="fill") +
  scale_fill_viridis(discrete = T)

### Start on the deSeq part of analysis ####
rownames(meta) <- unique(meta$BirdID_Tissue)
rownames(data) <- unique(data$X)
data <- subset(data, select = -c(X))
#meta <- subset(meta, select = -c(X))

meta <- meta[match(colnames(data), rownames(meta)),]

## Check that all column names from the data set are present as rownames in the metadata file
all(colnames(data) %in% rownames(meta))

## Check that the cols are in the same order in the data as the rownames in the metadata
all(colnames(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = meta, design = Tissue ~ Metabolic_State)

## Pre-filtering
## Taking out rows that have counts less than 10 reads total - recommended, not required, step
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

## Set the factor level to compare against (in our case, normothermy)
dds$Metabolic_State <- relevel(dds$Metabolic_State, ref="N")

## Run DESeq
dds <- DESeq(dds)

## Get normalized data, save it
normalized_counts <- counts(dds, normalized=TRUE)
norm_counts_df <- as.data.frame(normalized_counts)
gene <- rownames(norm_counts_df)
norm_counts_df <- cbind(gene, norm_counts_df)
head(norm_counts_df)

#### Needed for later analyses
## Make data long-form and merge with metadata file
meta2$Metabolic_State <- as.factor(meta2$Metabolic_State)
datlong <- norm_counts_df %>%
  #rename(gene = X) %>%
  gather(key = 'Sample', value= 'counts', -gene) %>%
  left_join(., meta2, by="Sample") %>%
  mutate(Metabolic_State = 
           fct_relevel(Metabolic_State, 
                       "N", "T", "D"))

## Make data long-form and merge with metadata file
foldlong <- foldchange %>%
  select(Tissue, gene, D_vs_N_FoldChange, D_vs_N_log2FoldChange, D_vs_N_padj,
         T_vs_D_FoldChange, T_vs_D_log2FoldChange, T_vs_D_padj, T_vs_N_FoldChange, T_vs_N_log2FoldChange, T_vs_N_padj) %>%
  gather(key = 'Measure', value= 'value', -c(gene,Tissue))


#### Don't re-run unless you need to make the norm counts df again
#write.csv(norm_counts_df, file = here("..//DESeq_Data_Brains_Jun2022/RNASeqBrain_NormCounts.csv"))

## For GSEA, make a file that has all the groups in a row in the order they are in in the columns of the the normalized 
## counts dataset. Then in Excel, add a top row that has e.g. <119 2 1> where 
## 119 is the number of samples; 2 is the number of groups and 1 is a constant for all files.
## Add a second row that has <# Trt1 Trt2 Cntrl> with spaces in between where each substring is the name of the treatment group 
## Encompassing all the cell values below
meta2$Tissue_State <- paste0(meta2$Tissue, "_", meta2$Metabolic_State)

phenotype_labs <- data.frame()

phenotype_labs_state <- rbind(phenotype_labs, meta2$Metabolic_State)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Brains_Jun2022//PhenoBrain_labs_state.csv"))

phenotype_labs_tissue <- rbind(phenotype_labs, meta2$Tissue)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Brains_Jun2022//PhenoBrain_labs_tissue.csv"))

phenotype_labs_state <- rbind(phenotype_labs, meta2$Tissue_State)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Brains_Jun2022//PhenoBrain_labs_tissue_state.csv"))


## If you have technical replicates, collapse them now. We don't.. we only have biological replicates. 
# Do not collapse those.


### can run here onwards every time
## Take a look at the data a bit
## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

## Save results
res_unshrunken <- results(dds)
res <- results(dds)

## Take a quick look at the results
res_unshrunken

## Summary of results
summary(res)

## Compare different pairs
res_ND <- results(dds, contrast=c("Metabolic_State", "D", "N"))
summary(res_ND)
res_NT <- results(dds, contrast=c("Metabolic_State", "T", "N"))
summary(res_NT)
res_TD <- results(dds, contrast=c("Metabolic_State", "D", "T"))
summary(res_TD)


### Subsetting results per tissue type
res_CB_ND <- results(dds[dds$Tissue=="CB",], contrast=c("Metabolic_State", "D", "N"))

res_DI_ND <- results(dds[dds$Tissue=="DI",], contrast=c("Metabolic_State", "D", "N"))

res_RT_ND <- results(dds[dds$Tissue=="RT",], contrast=c("Metabolic_State", "D", "N"))

res_DT_ND <- results(dds[dds$Tissue=="DT",], contrast=c("Metabolic_State", "D", "N"))


## Trying to make a res file per tissue type
results(dds, contrast=c("Metabolic_State", "D", "N"),)

# ### Adjusted p values of 0.01
# res_ND_0.01 <- results(dds, contrast=c("Metabolic_State", "D", "N"), alpha=0.01)
# summary(res_ND_0.01)
# res_NT_0.01 <- results(dds, contrast=c("Metabolic_State", "T", "N"), alpha=0.01)
# summary(res_NT_0.01)
# res_TD_0.01 <- results(dds, contrast=c("Metabolic_State", "D", "T"), alpha=0.01)
# summary(res_TD_0.01)

## Shrinking
res <- lfcShrink(dds, coef = 3, res = res)
#res

## Fold change results
mcols(res, use.names=T)
res %>% data.frame() %>% View()
## Summarize results
summary(res)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58 ## equal to fold change of 1.5

## We can easily subset the results table to only include 
## those that are significant using the filter() function, but first we will convert the results table into a tibble:
res_tb <- res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

## Just subset significantly different genes
sig <- res_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)


# Create tibbles including row names of normalized counts
normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# ## Order results by padj values
# top30_sig_genes <- res_tb %>% 
#   arrange(padj) %>% 	#Arrange rows by padj values
#   pull(gene) %>% 		#Extract character vector of ordered genes
#   head(n=30) 		#Extract the first 30 genes

# ## Order results by padj values, just ND and just upreg genes
# top30_sig_genes_upreg <- res_ND_tb %>% 
#   arrange(padj) %>% 	#Arrange rows by padj values
#   filter(log2FoldChange>1) %>% ### NOT WORKING YET
#   pull(gene) %>% 		#Extract character vector of ordered genes
#   head(n=30) 		#Extract the first 30 genes

## For enrichr analysis, just filtering out genes that have basemean >10, logFC >1 or < -1,
# And adjusted p-val of < 0.05

upreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

## For Transition vs. Normo
upreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

## For deep torpor vs. transition
upreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

## normalized counts for top significantly upregulated genes
top_upreg_ND_norm <- normalized_counts %>%
  filter(gene %in% upreg_ND)

# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_ND <- top_upreg_ND_norm %>%
  gather(colnames(top_upreg_ND_norm)[2:73], key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top_upreg_ND)

gathered_top_upreg_ND <- inner_join(meta2, gathered_top_upreg_ND, by="Sample")

## Plot this subset of these upreg genes
ggplot(gathered_top_upreg_ND) + #facet_grid(.~Metabolic_State) +
  geom_point(aes(x = gene, y = normalized_counts, color = Tissue)) +
  scale_y_log10() + 
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top sig upregulated genes Deep Torpor vs. Normo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis_d()


### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% upreg_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Upreg and downreg per tissue, just for ND comparison
upreg_CB_ND <-  res_CB_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_CB_ND <-  res_CB_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

upreg_DI_ND <-  res_DI_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_DI_ND <-  res_DI_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

upreg_RT_ND <-  res_RT_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_RT_ND <-  res_RT_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

upreg_DT_ND <-  res_DT_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_DT_ND <-  res_DT_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)


CB_samples <- meta2$Sample[meta2$Tissue=="CB"]
DI_samples <- meta2$Sample[meta2$Tissue=="DI"]
RT_samples <- meta2$Sample[meta2$Tissue=="RT"]
DT_samples <- meta2$Sample[meta2$Tissue=="DT"]


## CB ND UP
## normalized counts for top significantly upregulated genes
top_upreg_CB_ND_norm <- normalized_counts[,c("gene",CB_samples)] %>%
  filter(gene %in% upreg_CB_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_CB_ND <- top_upreg_CB_ND_norm %>%
  gather(colnames(top_upreg_CB_ND_norm)[2:ncol(top_upreg_CB_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top_upreg_CB_ND)

gathered_top_upreg_CB_ND <- inner_join(meta2, gathered_top_upreg_CB_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_CB <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% upreg_CB_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## CB DOWN
## normalized counts for top significantly downregulated genes
top_downreg_CB_ND_norm <- normalized_counts[,c("gene",CB_samples)] %>%
  filter(gene %in% downreg_CB_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_CB_ND <- top_downreg_CB_ND_norm %>%
  gather(colnames(top_downreg_CB_ND_norm)[2:ncol(top_downreg_CB_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top_downreg_CB_ND)

gathered_top_downreg_CB_ND <- inner_join(meta2, gathered_top_downreg_CB_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_CB <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% downreg_CB_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

## DI UP
## normalized counts for top significantly upregulated genes
top_upreg_DI_ND_norm <- normalized_counts[,c("gene",DI_samples)] %>%
  filter(gene %in% upreg_DI_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_DI_ND <- top_upreg_DI_ND_norm %>%
  gather(colnames(top_upreg_DI_ND_norm)[2:ncol(top_upreg_DI_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_upreg_DI_ND)

gathered_top_upreg_DI_ND <- inner_join(meta2, gathered_top_upreg_DI_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_DI <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% upreg_DI_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## DI DOWN
## normalized counts for top significantly downregulated genes
top_downreg_DI_ND_norm <- normalized_counts[,c("gene",DI_samples)] %>%
  filter(gene %in% downreg_DI_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_DI_ND <- top_downreg_DI_ND_norm %>%
  gather(colnames(top_downreg_DI_ND_norm)[2:ncol(top_downreg_DI_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_downreg_DI_ND)

gathered_top_downreg_DI_ND <- inner_join(meta2, gathered_top_downreg_DI_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_DI <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% downreg_DI_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## RT UP
## normalized counts for top significantly upregulated genes
top_upreg_RT_ND_norm <- normalized_counts[,c("gene",RT_samples)] %>%
  filter(gene %in% upreg_RT_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_RT_ND <- top_upreg_RT_ND_norm %>%
  gather(colnames(top_upreg_RT_ND_norm)[2:ncol(top_upreg_RT_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_upreg_RT_ND)

gathered_top_upreg_RT_ND <- inner_join(meta2, gathered_top_upreg_RT_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_RT <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% upreg_RT_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## RT DOWN
## normalized counts for top significantly downregulated genes
top_downreg_RT_ND_norm <- normalized_counts[,c("gene",RT_samples)] %>%
  filter(gene %in% downreg_RT_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_RT_ND <- top_downreg_RT_ND_norm %>%
  gather(colnames(top_downreg_RT_ND_norm)[2:ncol(top_downreg_RT_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_downreg_RT_ND)

gathered_top_downreg_RT_ND <- inner_join(meta2, gathered_top_downreg_RT_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_RT <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% downreg_RT_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

## DT UP
## normalized counts for top significantly upregulated genes
top_upreg_DT_ND_norm <- normalized_counts[,c("gene",DT_samples)] %>%
  filter(gene %in% upreg_DT_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_DT_ND <- top_upreg_DT_ND_norm %>%
  gather(colnames(top_upreg_DT_ND_norm)[2:ncol(top_upreg_DT_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_upreg_DT_ND)

gathered_top_upreg_DT_ND <- inner_join(meta2, gathered_top_upreg_DT_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_DT <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% upreg_DT_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## DT DOWN
## normalized counts for top significantly downregulated genes
top_downreg_DT_ND_norm <- normalized_counts[,c("gene",DT_samples)] %>%
  filter(gene %in% downreg_DT_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_DT_ND <- top_downreg_DT_ND_norm %>%
  gather(colnames(top_downreg_DT_ND_norm)[2:ncol(top_downreg_DT_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_downreg_DT_ND)

gathered_top_downreg_DT_ND <- inner_join(meta2, gathered_top_downreg_DT_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_DT <- normalized_counts[,c(1, 2:73)] %>% 
  filter(gene %in% downreg_DT_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Trying out heatmaps
## Tried a bit from but that didn't work 
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/


# my_sample_col <- data.frame(sample = rep(c("tumour", "normal"), c(4,2)))
# row.names(my_sample_col) <- colnames(data_subset)

annotation_col <- data.frame(
  Tissue = factor(meta2$Tissue), 
  Metab = factor(meta$Metabolic_State))
row.names(annotation_col) <- colnames(norm_sig_ND_upreg)

pheatmap(norm_sig_ND_upreg, annotation_col = annotation_col)


### Annotate our heatmap (optional)
annotation <- meta2 %>% 
  select(Sample, Metabolic_State, Tissue) %>% 
  data.frame(row.names = "Sample") %>%
  arrange(Metabolic_State)

### Set a color palette
heat_colors <- brewer.pal(9, "YlOrRd")

### Run pheatmap without col clusters
pheatmap(norm_sig_ND_upreg, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         cluster_cols = FALSE)


# create a `df` with the samples grouped in the same way you want to show
anno <- data.frame(SampleID = unique(gathered_top_upreg_ND$Sample), 
                   TissueState= unique(gathered_top_upreg_ND$Tissue_State))

# set rownames so that anno and your data can be matched
rownames(anno) <- gathered_top_upreg_ND$Sample


# set colours
anno_colors <- list(
  TissueState = c(CB_N = "#1B9E77", Gut1_N = "#D95F02", Gut2_N = "#823de9", Gut3_N = "#7855ce",
                  DI_N = "#6e6eb2", RT_N = "#648697", DT_N = "#599e7c", Gut3_D = "#f94144", 
                  Gut1_D = "#f3722c", Gut2_D = "#f8961e", Gut3_T = "#f9844a", RT_T = "#f9c74f",
                  RT_D = "#90be6d", DT_D = "#43aa8b", DI_D = "#4d908e", Gut2_T = "#577590", 
                  Gut1_T = "#277da1", CB_D = "#004e64", CB_T = "#ffba08", DT_T = "#f7b2bd", DI_T = "#c60f7b"))


## THIS is the one I used in the lab meeting presentation on 4/21/22
# set colours
anno_colors <- list(
  Tissue = c(CB = "#f3722c", Gut1 = "#D95F02", Gut2 = "#f8961e", Gut3 = "#577590",
             DI = "#ffba08", RT = "#004e64", DT = "#599e7c"),
  Metabolic_State = c(N = "#f9c74f", D = "#43aa8b", T = "#277da1"))

pheatmap(norm_sig_ND_upreg, 
         annotation_col = annotation, 
         annotation_colors = anno_colors,
         #color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         border_color = NA, 
         scale = "row", 
         fontsize_row = 10,
         fontsize_col = 10,
         height = 20,
         cluster_cols = F,
         legend = T,
         fontsize = 20)


annotation_col <- data.frame(
  Tissue = factor(meta2$Tissue), 
  Metab = factor(meta$Metabolic_State))

#rownames(annotation_col) = paste("Test", 1:10, sep = "")

# annotation_row <- data.frame(
#   GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
# )
# rownames(annotation_row) = paste("Gene", 1:20, sep = "")

ann_colors = list(
  Tissue = c(CB = "#1B9E77", Gut1 = "#D95F02", Gut2 = "#823de9", Gut3 = "#7855ce",
             DI = "#6e6eb2", RT = "#648697", DT = "#599e7c"),
  Metab = c(N = "red", T = "black", D = "purple")
)

#c("#823de9", "#7855ce", "#6e6eb2", "#648697", "#599e7c", "#4fb760", "#45cf45")

## Tissue-specific plots


# set colours
anno_colors <- list(
  # Tissue = c(CB = "#f3722c", Gut1 = "#D95F02", Gut2 = "#f8961e", Gut3 = "#577590",
  #            DI = "#ffba08", RT = "#004e64", DT = "#599e7c"),
  Metabolic_State = c(N = "#f9c74f", T = "violet", D = "#599e7c"))

## CB heatmap
annotation_CB <- meta2 %>%
  filter(Tissue == "CB") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## Ordering cols by metabolic state
## USE THIS
# 1) reorder the matrix based in the annotation
CB_upND_ordered <- top_upreg_CB_ND_norm[, rownames(annotation_CB)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(CB_upND_ordered,
                   annotation_col = annotation_CB, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   #cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## Downreg CB ND
# 1) reorder the matrix based in the annotation
CB_dnND_ordered <- norm_sig_ND_downreg_CB[, rownames(annotation_CB)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(CB_dnND_ordered,
                   annotation_col = annotation_CB, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## DI heatmaps
## Ordering cols by metabolic state
annotation_DI <- meta2 %>%
  filter(Tissue == "DI") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## USE THIS
# 1) reorder the matrix based in the annotation
DI_upND_ordered <- norm_sig_ND_upreg_DI[, rownames(annotation_DI)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(DI_upND_ordered,
                   annotation_col = annotation_DI, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)

## Downreg DI ND
# 1) reorder the matrix based in the annotation
DI_dnND_ordered <- norm_sig_ND_downreg_DI[, rownames(annotation_DI)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(DI_dnND_ordered,
                   annotation_col = annotation_DI, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## RT heatmaps
## Ordering cols by metabolic state
annotation_RT <- meta2 %>%
  filter(Tissue == "RT") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## USE THIS
# 1) reorder the matrix based in the annotation
RT_upND_ordered <- norm_sig_ND_upreg_RT[, rownames(annotation_RT)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(RT_upND_ordered,
                   annotation_col = annotation_RT, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)

## Downreg RT ND
# 1) reorder the matrix based in the annotation
RT_dnND_ordered <- norm_sig_ND_downreg_RT[, rownames(annotation_RT)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(RT_dnND_ordered,
                   annotation_col = annotation_RT, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## DT heatmaps
## Ordering cols by metabolic state
annotation_DT <- meta2 %>%
  filter(Tissue == "DT") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## USE THIS
# 1) reorder the matrix based in the annotation
DT_upND_ordered <- norm_sig_ND_upreg_DT[, rownames(annotation_DT)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(DT_upND_ordered,
                   annotation_col = annotation_DT, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)

## Downreg DT ND
# 1) reorder the matrix based in the annotation
DT_dnND_ordered <- norm_sig_ND_downreg_DT[, rownames(annotation_DT)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(DT_dnND_ordered,
                   annotation_col = annotation_DT, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## Volcano plot normothermy vs deep torpor
#### These are probab
ND_volcano <- EnhancedVolcano(res_ND,
                              lab = rownames(res_ND),
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'Deep torpor vs. Normothermy',
                              subtitle = "Enhanced volcano, p cutoff = 0.05",
                              pCutoff = 0.05)


## Volcano plot normothermy vs deep torpor
NT_volcano <- EnhancedVolcano(res_NT,
                              lab = rownames(res_NT),
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'Transition vs. Normothermy',
                              subtitle = "Enhanced volcano, p cutoff = 0.05",
                              pCutoff = 0.05)

## Volcano plot normothermy vs deep torpor
TD_volcano <- EnhancedVolcano(res_TD,
                              lab = rownames(res_TD),
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'Deep Torpor vs. Transition',
                              subtitle = "Enhanced volcano, p cutoff = 0.05",
                              pCutoff = 0.05)

grid.arrange(ND_volcano, NT_volcano, TD_volcano, ncol=3)

## Use adjusted p-vals to make plots instead of nominal p-values
ND_adj_volcano <- EnhancedVolcano(res_ND, lab = rownames(res_ND),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  #xlab = bquote(~Log[2]~ "fold change"),
                                  ylab = bquote(-Log[10]~adjusted~italic(P)),
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  #transcriptLabSize = 3.0,
                                  colAlpha = 1,
                                  title = 'Deep Torpor vs. Normothermy',
                                  subtitle = "Enhanced volcano, adj p cutoff = 0.05")
#legend=c("NS","Log2 FC","Adjusted p-value",
#        "Adjusted p-value & Log2 FC"),
#legendPosition = "bottom",
#legendLabSize = 10,
#legendIconSize = 3.0)

NT_adj_volcano <- EnhancedVolcano(res_NT, lab = rownames(res_NT),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  #xlab = bquote(~Log[2]~ "fold change"),
                                  ylab = bquote(-Log[10]~adjusted~italic(P)),
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  #transcriptLabSize = 3.0,
                                  colAlpha = 1,
                                  title = 'Transition vs. Normothermy',
                                  subtitle = "Enhanced volcano, adj p cutoff = 0.05")

TD_adj_volcano <- EnhancedVolcano(res_TD, lab = rownames(res_TD),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  #xlab = bquote(~Log[2]~ "fold change"),
                                  ylab = bquote(-Log[10]~adjusted~italic(P)),
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  #transcriptLabSize = 3.0,
                                  colAlpha = 1,
                                  title = 'Deep Torpor vs. Transition',
                                  subtitle = "Enhanced volcano, adj p cutoff = 0.05")

grid.arrange(ND_adj_volcano, NT_adj_volcano, TD_adj_volcano, ncol=3)

ND_adj_volcano
NT_adj_volcano
TD_adj_volcano


## Tissue specific
## Use adjusted p-vals to make plots instead of nominal p-values
ND_adj_volcano_CB <- EnhancedVolcano(res_CB_ND, lab = rownames(res_CB_ND),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        #xlab = bquote(~Log[2]~ "fold change"),
                                        ylab = bquote(-Log[10]~adjusted~italic(P)),
                                        pCutoff = 0.05,
                                        FCcutoff = 1.0,
                                        #transcriptLabSize = 3.0,
                                        colAlpha = 1,
                                        title = 'CB: Deep Torpor vs. Normothermy',
                                        subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_CB


ND_adj_volcano_DI <- EnhancedVolcano(res_DI_ND, lab = rownames(res_DI_ND),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        #xlab = bquote(~Log[2]~ "fold change"),
                                        ylab = bquote(-Log[10]~adjusted~italic(P)),
                                        pCutoff = 0.05,
                                        FCcutoff = 1.0,
                                        #transcriptLabSize = 3.0,
                                        colAlpha = 1,
                                        title = 'DI: Deep Torpor vs. Normothermy',
                                        subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_DI

ND_adj_volcano_RT <- EnhancedVolcano(res_RT_ND, lab = rownames(res_RT_ND),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        #xlab = bquote(~Log[2]~ "fold change"),
                                        ylab = bquote(-Log[10]~adjusted~italic(P)),
                                        pCutoff = 0.05,
                                        FCcutoff = 1.0,
                                        #transcriptLabSize = 3.0,
                                        colAlpha = 1,
                                        title = 'RT: Deep Torpor vs. Normothermy',
                                        subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_RT

ND_adj_volcano_DT <- EnhancedVolcano(res_DT_ND, lab = rownames(res_DT_ND),
                                       x = 'log2FoldChange',
                                       y = 'padj',
                                       #xlab = bquote(~Log[2]~ "fold change"),
                                       ylab = bquote(-Log[10]~adjusted~italic(P)),
                                       pCutoff = 0.05,
                                       FCcutoff = 1.0,
                                       #transcriptLabSize = 3.0,
                                       colAlpha = 1,
                                       title = 'DT: Deep Torpor vs. Normothermy',
                                       subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_DT


## From Helen Chmura's 2022 Communications Biology paper
#mRNA profiling by in-situ hybridization in the ventral hypothalamus and pituitary PT 
# revealed significant suppression of TSHβ, Dio2, and elevation of Dio3 during early hibernation in females. 
# Strikingly, in late hibernation, we observed significant elevation in TSHβ and Dio2 and suppression of Dio3.
## geom_density of DIO3
datlong %>%
  filter(gene %in% c("DIO2", "DIO3")) %>%
  ggplot(., aes(x=counts, fill=Metabolic_State)) +
  geom_density(alpha=0.5) + my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states")

# Boxplot/Violin plot of DIO3
datlong %>%
  filter(gene == "DIO2") %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + geom_point() +
  #geom_violin() +
  my_theme + facet_grid(.~Tissue, scales = "free") +
  ggtitle("DIO3 gene expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

datlong %>%
  filter(gene %in% c('DIO2', "DIO3")) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_grid(gene~Tissue, scales = "free") +
  ggtitle("DIO3 gene expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")


## geom_density of EYA3
datlong %>%
  filter(gene == "EYA3") %>%
  ggplot(., aes(x=counts, fill=Metabolic_State)) +
  geom_density(alpha=0.5) + my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states")

# Boxplot/Violin plot of EYA3
datlong %>%
  filter(gene == "FGF21") %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("DIO3 gene expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

## geom_density of clock genes (there is no BMAL in this dataset)
# and feeding behavior and obesity
clockgenes <- c('CLOCK', 'PER2', 'PER3', 'CRY1', 'CRY2')
## Feeding and obesity genes from hibernation paper https://onlinelibrary.wiley.com/doi/10.1111/gbb.12199
foodgenes <- c('VGF', 'TRH', 'LEPR', 'ADIPOR2', 'IRS2')
## Metabolism genes from Figure 3 in https://www.nature.com/articles/s41598-018-31506-2/figures/3
# This paper also has clock genes
metabgenes <- c('PPARA', 'SIRT1', 'LEPR')
## Insulin related genes from brown bear hibernation
## https://academic.oup.com/icb/advance-article/doi/10.1093/icb/icac093/6609442?login=true
insulingenes <- c('APPL1', 'ATP6AP1', 'ATP6V0A1', 'GRB14', 'IGF2') #AKT2 isn't present
insulingenes2 <- c('INSR', 'PIK3R1', 'PPFIBP1', 'SHC1', 'SLC18B1', 'SORBS1') # 'TCIRG1' isn't present

datlong %>%
  filter(gene == clockgenes) %>%
  ggplot(., aes(x=counts, fill=Metabolic_State)) +
  geom_density(alpha=0.5) + my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states")

# Boxplot/Violin plot of clock genes
datlong %>%
  filter(gene == clockgenes) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(Tissue~gene, scales = "free") +
  ggtitle("Clock genes expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

# Boxplot/Violin plot of clock genes just in DI
datlong %>%
  filter(gene == clockgenes, Tissue == 'DI') %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + geom_point(aes(col=BirdID)) +
  #geom_violin() +
  my_theme + facet_wrap(.~gene, scales = "free") +
  ggtitle("Clock genes expression in the diencephalon across metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

# Boxplot/Violin plot of insulin genes
## Not much going on, and low-ish levels of expression, low sample sizes
datlong %>%
  filter(gene == insulingenes) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_grid(gene~Tissue, scales = "free") +
  ggtitle("Insulin genes expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

## 6 other insulin genes
## Not much going on, and low-ish levels of expression, low sample sizes
# Except.. SLC18B1 for DI, which has some pattern, and PPFIBP1 for DT and RT, maybe SORBS1 for DT and RT
datlong %>%
  filter(gene == insulingenes2) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_grid(gene~Tissue, scales = "free") +
  ggtitle("Insulin genes expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

datlong %>%
  filter(gene == metabgenes) %>%
  ggplot(., aes(x=counts, fill=Metabolic_State)) +
  geom_density(alpha=0.5) + my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states")

# Boxplot/Violin plot of metabolism genes
datlong %>%
  filter(gene == metabgenes) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(gene~Tissue, scales = "free") +
  ggtitle("Metabolism genes expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")


#old, manually selected from each tissue type
#genes.of.interest <- c('NR1D1', 'SGK1', 'RSRP1', 'DUSP1', 'NME5', 'TMEM39B', 'BHLHE40', 'LOC115598688')

## pull together all the top upreg and downreg genes from each pairwise comparison
## and remove the ones which are not annotated
top_genes <- c(upreg_ND, upreg_NT, upreg_TD, downreg_ND, downreg_NT, downreg_TD)
top_genes_annotated <- top_genes[!grepl("LOC", top_genes)]

upreg_ND[!grepl("LOC", upreg_ND)]
upreg_NT[!grepl("LOC", upreg_NT)]
upreg_TD[!grepl("LOC", upreg_TD)]

downreg_ND[!grepl("LOC", downreg_ND)]
downreg_NT[!grepl("LOC", downreg_NT)]
downreg_TD[!grepl("LOC", downreg_TD)]


upreg_CB_ND[!grepl("LOC", upreg_CB_ND)]
upreg_DI_ND[!grepl("LOC", upreg_DI_ND)]
upreg_RT_ND[!grepl("LOC", upreg_RT_ND)]
upreg_DT_ND[!grepl("LOC", upreg_DT_ND)]

downreg_CB_ND[!grepl("LOC", downreg_CB_ND)]
downreg_DI_ND[!grepl("LOC", downreg_DI_ND)]
downreg_RT_ND[!grepl("LOC", downreg_RT_ND)]
downreg_DT_ND[!grepl("LOC", downreg_DT_ND)]

datlong %>%
  filter(gene %in% top_genes_annotated) %>%
  ggplot(., aes(x = Tissue, y = gene, fill = counts)) +
  geom_tile() + my_theme +
  scale_fill_gradient(low = 'white', high = 'red') + facet_grid(.~Metabolic_State, scales = "free")


## Examples
genes.of.interest <- c('NR1D1', 'SGK1', 'RSRP1', 'DUSP1', 'NME5', 'TMEM39B', 'BHLHE40', 'LOC115598688')
foldlong %>%
  filter(Measure %in% c('D_vs_N_log2FoldChange', 'T_vs_N_log2FoldChange', 'T_vs_D_log2FoldChange')) %>%
  ggplot(., aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')

foldlong %>%
  filter(Measure %in% c('D_vs_N_log2FoldChange', 'T_vs_N_log2FoldChange', 'T_vs_D_log2FoldChange')) %>%
  ggplot(., aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')


sig_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 | log2FoldChange < -1) %>%
  filter(padj < 0.05 & baseMean > 10) %>%
  left_join(., meta2, by="Sample") %>%
  mutate(Metabolic_State = 
           fct_relevel(Metabolic_State, 
                       "N", "T", "D")) 
arrange(padj) 

ggplot(sig_ND, aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')


res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene")
filter(Measure %in% c('D_vs_N_log2FoldChange', 'T_vs_N_log2FoldChange', 'T_vs_D_log2FoldChange')) %>%
  ggplot(., aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')


sig_top_all_comparisons  <- datlong %>%
  filter(gene %in% top_genes_annotated)

sig_top_all_comparisons <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% top_genes_annotated) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

pheatmap(sig_top_all_comparisons, 
         annotation_col = annotation, 
         annotation_colors = anno_colors,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         border_color = NA, 
         scale = "row", 
         fontsize_row = 15,
         fontsize_col = 10,
         height = 20,
         cluster_cols = F,
         legend = T,
         fontsize = 20)




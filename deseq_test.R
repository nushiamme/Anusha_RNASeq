#### DeSeq for GSEA analysis
## Anusha Shankar, github: nushiamme
## Code started April 5, 2022

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


## If you opened the .Rproj file from the OneDrive folder, reading in these files will work- 
## the relative paths should be the same
data <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//RNASeq_rawcounts_data.csv"))
meta <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//AnushaShankar_RNASeq_metadata.csv"))
star_alignment <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//star_alignment_plot.csv"))
star <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//star_gene_counts.csv"))
foldchange <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//TopFoldChanges.csv"))

my_theme <- theme_classic(base_size = 15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.key.height = unit(2, "line"))

## Viridis colors
my_gradient <- c("#823de9", "#7855ce", "#6e6eb2", "#648697", "#599e7c", "#4fb760", "#45cf45")

## Checking on the star alignment
head(star)
star <- star %>%
  rename(Sample = ï..Sample) %>%
  mutate(Sample = gsub("c", "", Sample))

meta2 <- meta %>%
  rename(Sample = X)

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

## Make data long-form and merge with metadata file
datlong <- data %>%
  rename(gene = X) %>%
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


rownames(meta) <- unique(meta$X)
rownames(data) <- unique(data$X)
data <- subset(data, select = -c(X))
meta <- subset(meta, select = -c(X))

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
write.csv(norm_counts_df, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//RNASeq_NormCounts.csv"))

## For GSEA, make a file that has all the groups in a row in the order they are in in the columns of the the normalized 
## counts dataset. Then in Excel, add a top row that has e.g. <119 2 1> where 
## 119 is the number of samples; 2 is the number of groups and 1 is a constant for all files.
## Add a second row that has <# Trt1 Trt2 Cntrl> with spaces in between where each substring is the name of the treatment group 
## Encompassing all the cell values below
meta2$Tissue_State <- paste0(meta2$Tissue, "_", meta2$Metabolic_State)

phenotype_labs <- data.frame()

phenotype_labs_state <- rbind(phenotype_labs, meta2$Metabolic_State)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//Pheno_labs_state.csv"))

phenotype_labs_tissue <- rbind(phenotype_labs, meta2$Tissue)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//Pheno_labs_tissue.csv"))

phenotype_labs_state <- rbind(phenotype_labs, meta2$Tissue_State)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//Pheno_labs_tissue_state.csv"))


## If you have technical replicates, collapse them now. We don't.. we only have biological replicates. 
# Do not collapse those.



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
res_tb <- res %>%
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

## Order results by padj values
top30_sig_genes <- res_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=30) 		#Extract the first 30 genes

# ## Order results by padj values, just ND and just upreg genes
# top30_sig_genes_upreg <- res_ND_tb %>% 
#   arrange(padj) %>% 	#Arrange rows by padj values
#   filter(log2FoldChange>1) %>% ### NOT WORKING YET
#   pull(gene) %>% 		#Extract character vector of ordered genes
#   head(n=30) 		#Extract the first 30 genes

## For enrichr analysis, just filtering out genes that have basemean >10, logFC >1 or < -1,
# And adjusted p-val of < 0.05

downreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10)

upreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10)

## For Transition vs. Normo
downreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10)

upreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10)

## For deep torpor vs. transition
downreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10)

upreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10)


## normalized counts for top 20 significant genes
top30_sig_norm <- normalized_counts %>%
  filter(gene %in% top30_sig_genes)

# Gathering the columns to have normalized counts to a single column
gathered_top30_sig <- top30_sig_norm %>%
  gather(colnames(top30_sig_norm)[2:120], key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top30_sig)

gathered_top30_sig <- inner_join(meta2, gathered_top30_sig, by="Sample")

## Meh plot using ggplot2
ggplot(gathered_top30_sig) + #facet_grid(.~Metabolic_State) +
  geom_point(aes(x = gene, y = normalized_counts, color = Tissue)) +
  scale_y_log10() + 
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 30 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis_d()



### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% sig$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

### Annotate our heatmap (optional)
annotation <- meta2 %>% 
  select(Sample, Tissue) %>% 
  data.frame(row.names = "Sample")

### Set a color palette
heat_colors <- brewer.pal(9, "YlOrRd")

### Run pheatmap without col clusters
pheatmap(norm_sig, 
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

### Run pheatmap with col clusters
pheatmap(norm_sig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         cluster_cols = T)

annotation_col <- data.frame(
  Tissue = factor(meta2$Tissue), 
  Metab = factor(meta$Metabolic_State))

#rownames(annotation_col) = paste("Test", 1:10, sep = "")

# annotation_row <- data.frame(
#   GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
# )
# rownames(annotation_row) = paste("Gene", 1:20, sep = "")

ann_colors = list(
  Tissue = c(Liver = "#1B9E77", Gut1 = "#D95F02", Gut2 = "#823de9", Gut3 = "#7855ce",
             Heart = "#6e6eb2", Lungs = "#648697", Pect = "#599e7c"),
  Metab = c(N = "#7570B3", T = "#E7298A", D = "#66A61E")
)

#c("#823de9", "#7855ce", "#6e6eb2", "#648697", "#599e7c", "#4fb760", "#45cf45")

pheatmap(norm_sig, annotation_col = annotation_col,#, #annotation_row = annotation_row, 
         annotation_colors = ann_colors,
         #cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         cluster_cols = FALSE)


pheatmap(norm_sig, annotation_col = annotation_col,
         annotation_colors = ann_colors,
         #cluster_rows = T, 
         color = heat_colors, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)



## Volcano plot normothermy vs deep torpor
#### MAKE SURE THIS IS ADJUSTED P-VAL
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

## geom_density
datlong %>%
  filter(gene == "RSRP1") %>%
  ggplot(., aes(x=counts, fill=Metabolic_State)) +
  geom_density(alpha=0.5) + my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states")

# Boxplot/Violin plot
datlong %>%
  filter(gene == "RSRP1") %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

genes.of.interest <- c('NR1D1', 'SGK1', 'RSRP1', 'DUSP1', 'NME5', 'TMEM39B', 'BHLHE40', 'LOC115598688')
datlong %>%
  filter(gene %in% genes.of.interest) %>%
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


## For enrichr analysis, just filtering out genes that have basemean >10, logFC >1 or < -1,
# And adjusted p-val of < 0.05

downreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10)

upreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10)

## For Transition vs. Normo
downreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10)

upreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10)

## For deep torpor vs. transition
downreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10)

upreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -1 & padj < 0.05 & baseMean > 10)



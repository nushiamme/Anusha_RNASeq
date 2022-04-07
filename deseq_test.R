#### DeSeq for GSEA analysis
## Anusha Shankar, github: nushiamme
## Code started April 5, 2022

## Initially followed this tutorial:
## https://www.youtube.com/watch?v=OzNzO8qwwp0

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('DESeq2')
# BiocManager::install('EnhancedVolcano') 
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

## If you opened the .Rproj file from the OneDrive folder, reading in these files will work- 
## the relative paths should be the same
data <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//RNASeq_rawcounts_data.csv"))
meta <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//AnushaShankar_RNASeq_metadata.csv"))
star_alignment <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//star_alignment_plot.csv"))
star <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//star_gene_counts.csv"))
foldchange <- read.csv(here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//TopFoldChanges.csv"))

my_theme <- theme_classic(base_size = 15) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.key.height = unit(2, "line"))


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

## If you have technical replicates, collapse them now. We don't.. we only have biological replicates. 
# Do not collapse those.

## Save results
res <- results(dds)

## Take a quick look at the results
res

## Summary of results
summary(res)

## Compare different pairs
res_ND <- results(dds, contrast=c("Metabolic_State", "D", "N"))
summary(res_ND)
res_NT <- results(dds, contrast=c("Metabolic_State", "T", "N"))
summary(res_NT)
res_TD <- results(dds, contrast=c("Metabolic_State", "D", "T"))
summary(res_TD)

### Adjusted p values of 0.01
res_ND_0.01 <- results(dds, contrast=c("Metabolic_State", "D", "N"), alpha=0.01)
summary(res_ND_0.01)
res_NT_0.01 <- results(dds, contrast=c("Metabolic_State", "T", "N"), alpha=0.01)
summary(res_NT_0.01)
res_TD_0.01 <- results(dds, contrast=c("Metabolic_State", "D", "T"), alpha=0.01)
summary(res_TD_0.01)

## Shrinking??
#res <- lfcShrink(dds, coef = 2, res = res)
#res

## MA plot
plotMA(res_ND)

## Volcano plot normothermy vs deep torpor
ND_volcano <- EnhancedVolcano(res_ND,
                lab = rownames(res_ND),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Deep torpor vs. Normothermy',
                subtitle = "Enhanced volcano, p cutoff = 10e-10",
                pCutoff = 10e-10)


## Volcano plot normothermy vs deep torpor
NT_volcano <- EnhancedVolcano(res_NT,
                lab = rownames(res_NT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Transition vs. Normothermy',
                subtitle = "Enhanced volcano, p cutoff = 10e-10",
                pCutoff = 10e-10)

## Volcano plot normothermy vs deep torpor
TD_volcano <- EnhancedVolcano(res_TD,
                lab = rownames(res_TD),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Deep Torpor vs. Transition',
                subtitle = "Enhanced volcano, p cutoff = 10e-10",
                pCutoff = 10e-10)

grid.arrange(ND_volcano, NT_volcano, TD_volcano, ncol=3)




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


genes.of.interest <- c('NR1D1', 'SGK1', 'RSRP1', 'DUSP1', 'NME5', 'TMEM39B', 'BHLHE40', 'LOC115598688')
foldlong %>%
  filter(Measure %in% c('D_vs_N_log2FoldChange', 'T_vs_N_log2FoldChange', 'T_vs_D_log2FoldChange')) %>%
  ggplot(., aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')




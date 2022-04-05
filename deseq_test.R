#### DeSeq for GSEA analysis
## Anusha Shankar, github: nushiamme
## Code started April 5, 2022

library(DESeq2)


data <- read.csv("C://Users//nushi//OneDrive - Cornell University//Anusha_personal//Cornell//TorporFieldwork_IR//RNASeq_andQC_data//DESeq_Data_Mar2022_all tissues//all tissues//final//RNASeq_rawcounts_data.csv")
meta <- read.csv("C://Users//nushi//OneDrive - Cornell University//Anusha_personal//Cornell//TorporFieldwork_IR//RNASeq_andQC_data//DESeq_Data_Mar2022_all tissues//all tissues//final//AnushaShankar_RNASeq_metadata.csv")

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
dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Metabolic_State)

# see vignette for suggestions on generating
# count tables from RNA-Seq data
cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
cond <- factor(rep(1:2, each=5))

# object construction
dds <- DESeq2::DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)

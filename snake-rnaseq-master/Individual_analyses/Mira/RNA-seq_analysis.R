# Rhesus Macaque brain data RNA-seq quantification analysis 

setwd("~/Desktop/extra_dennis_lab/rhesus_tag-seq")

install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
library(GenomicFeatures)
library(tximport)
library(DESeq2)

txdb = makeTxDbFromGFF("Macaca_mulatta.Mmul_8.0.1.97.gff3")
columns(txdb)

k = keys(txdb, keytype="TXNAME") # extract transcript names
tx2gene <- select(txdb, k, "GENEID", "TXNAME") # for each transcript extract the gene


samples = sapply(read.table("samples.txt"), as.character) # list of samples; defaults to factor
files = file.path("results", samples, "quant.sf") # list of filepaths
names(files) = samples

# transcripts and gene names used by salmon are ensemble names, while those in txdb are actual gene and tx names 
# convert gene and tx names to ensemble names 

# read in list of ensemble names and their corresponding tx and gene names
mart = read.table("mart_export.txt", sep = '\t', header = TRUE)

# subset mart list by transcript names found in tx2gene
mart_sub = mart[mart$Transcript.name %in% tx2gene$TXNAME,]

# merge mart list with tx2gene data frame by transcript name which they have in common
new_tx = merge(mart_sub, tx2gene, by.x = "Transcript.name", by.y = "TXNAME" )
head(new_tx)

newtx2gene = subset(new_tx, select=c("Transcript.stable.ID.version", "GENEID"))
names(newtx2gene) = c("TXNAME", "GENEID" )
head(newtx2gene)

# newtx2gene is tx2gene with ensemble name substituted for TXNAME, so that it can match with quant files from salmon 

txi = tximport(files, type = "salmon", tx2gene = newtx2gene)
names(txi)
head(txi)

# abundance = TPM 
# use count data as input for these methods
head(samples)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~V1)
vsd <- vst(dds)
plotPCA(vsd)


 
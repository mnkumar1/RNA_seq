#!/usr/bin/env Rscript
# convert Salmon output from multiple samples into gene-level TPM table
# Usage:
# Rscript tximport.R [list of Salmon output directory names] [path to dir containing Salmon output] ["human", "chimpanzee", or path to tx2gene table] [optional: "afterBar"]

library(tximport)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)

#afterBar = F
#if (args[3]=="afterBar"){
#	afterBar = T
#}

# read in conversion table
print("READING tx2gene...")
if (args[3]=="human") {
	tx2gene = read.table("/share/dennislab/users/cshew/dreg/transcriptome_update/tximport/tx2gene.tsv", header=T)
	} else if (args[3]=="chimpanzee") {
	tx2gene = read.table("/share/dennislab/users/cshew/dreg/interspp/4_orthologs/chimp_tx2gene.tsv", header=T)
	} else {
	tx2gene = read.table(args[3], header=T)
}

# compile Salmon output in tabular format per sample
print("GENERATING TABLE FOR ALL SAMPLES...")
samples = sapply(read.table(args[1]), as.character) # list of samples; defaults to factor
files = file.path(args[2], samples, "quant.sf") # list of filepaths
names(files) = samples

print(files)
# convert to gene level
print("TXIMPORT...")
txi = tximport(files, type = "salmon", tx2gene = tx2gene)

# generate TPM table
print("MAKING OUTPUT TABLE...")
tpm = txi$abundance
tpm = as.data.frame(cbind(gene=rownames(tpm), tpm, median=rowMedians(tpm), mean=rowMeans(tpm))) # add column for genes, mean, and median
tpm = transform(tpm, sd=apply(tpm[,2:(ncol(tpm)-1)],1, sd, na.rm = TRUE)) # add SD

# save output
print("SAVING...")
write.table(tpm, "allSamplesTPM.tsv", row.names=F, quote=F, sep='\t')

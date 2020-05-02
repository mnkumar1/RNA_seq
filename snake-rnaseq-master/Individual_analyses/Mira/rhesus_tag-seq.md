# Rhesus Tag-Seq analysis

The goal of this project is to analyze Rhesus macaque tag-seq expression data from different brain regions

The analysis can be found at
```bash
/share/dennislab/users/miramastoras/rhesus_tag-seq/
```

## Data Download

Rhesus macaque transcriptome was downloaded from the Ensembl website using the following command
```bash
wget ftp://ftp.ensembl.org/pub/release-97/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_8.0.1.cdna.all.fa.gz
```
Tag-seq data is located at:

"Failed" run:
```bash
/share/dennislab-backedup/illumina/noncoding/190712_rhesus_3prime-tag (edited)
```
"Good" run:
```bash
/share/dennislab-backedup/illumina/noncoding/190722_rhesus_3prime-tag
```

Naming scheme:
There are four individuals for this experiment:
60-1,
60-2,
60-3,
150-2

60-1-1 and 60-2-1 are biological replicates: the third number  corresponds to the region of the brain, first number corresponds to the individual.


## Description of Snakefile rules and parameter choices

Trimmomatic:

HEADCROP:22 the UMI barcodes are present in the first 22 bases, so those are trimmed off first
https://dnatech.genomecenter.ucdavis.edu/faqs/where-can-i-find-the-umis-in-the-tag-seq-data-when-and-how-should-i-trim-my-tag-seq-data/


## Running the snakemake

```snakemake -n -p --config filenames="fastq_files.txt" reference="reference/Macaca_mulatta.Mmul_8.0.1.cdna.all.fa.gz" -j 20```

filenames = text file containing list of all fastq files to run pipeline on

## Results 

Salmon output quant files are located here 
```/share/dennislab/users/miramastoras/rhesus_tag-seq/results/quant```

## Analysis in R 

See ``RNA-seq_analysis.R`` for R script generating initial PCA plot of the Salmon results


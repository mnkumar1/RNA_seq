## Generic RNA-seq pipeline for Dennis lab

Steps:
- qc
- trimming
- qc trimming
- salmon index
- salmon quant
- collate salmon quant

Steps to Run:
1. Working directory must contain
    - Snakefile
    - scripts/tximport.R (tximport.R contained in directory "scripts")
    - config file

2. Config File Format:
    ```yaml
    # EXAMPLE config file for snake-combine-rnaseq

    reference: path/to/reference/

    input_dir: path/to/fastq/files/  

    end: paired or single # must write either paired or single. case sensitive

    tx2gene: path to tx2gene table for summing TPM to gene level

    transcriptome: path to transcriptome

    samples: file containing sample names. if paired end, do not include the _1 or _2
    ```
3. Activate conda environment
```bash
source /share/dennislab/programs/dennis-miniconda/etc/profile.d/conda.sh
conda activate rna-seq
```

4. Run snakemake with the following command:
    ```bash
    /share/dennislab/programs/dennis-miniconda/bin/snakemake --configfile config.yaml -p -j 20
    ```

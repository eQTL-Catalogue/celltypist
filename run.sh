#!/bin/bash

#SBATCH --job-name=celltypist
#SBATCH --partition=amd
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

nextflow -log logs/.nextflow.log run main.nf -profile tartu_hpc \
    --samples assets/input_examples/samplesheet.tsv \
    --model Immune_All_Low \
    --outdir results \
    -resume

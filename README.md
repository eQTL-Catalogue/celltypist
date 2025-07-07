# celltypist
Nextflow workflow for scRNAseq cell type annotation with CellTypist

Prerequisites:
* Nextflow
* Singularity

An example of the nextflow run script is in [`run.sh`](run.sh)

Parameters:
* `--samples` TSV containing the sample IDs, paths to the normalized (scaled to counts per 10,000 and natural log with pseudocount) `.h5ad` count matrices and cell metadata (example: [`assets/input_examples/samplesheet.tsv`](assets/input_examples/samplesheet.tsv))
* `--model` name of the CellTypist model to use
* `--outdir` name of the output directory

Outputs:
* `annotate` contains the count matrices with majority voted cell types in `.h5ad` files, cell type distributions, UMAPs, and CellTypist outputs
* `celltype_concordance` contains comparisons of annotated cell types to the metadata - cell type distributions and a confusion matrix

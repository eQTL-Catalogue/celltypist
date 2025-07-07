#!/usr/bin/env python

import os
import sys
import celltypist
import pandas as pd
import scanpy as sc
from celltypist import models
import matplotlib.pyplot as plt

countmtx = sc.read(sys.argv[1])
model_name = sys.argv[2]
outpath = sys.argv[3]
cell_meta = pd.read_csv(sys.argv[4], sep='\t')
sample = sys.argv[5]
cell_meta_out = sys.argv[6]

#print('Input count matrix:')
#print(countmtx)

model = models.Model.load(model=f'{model_name}.pkl')
print(f'\nModel used for predicting celltypes: {model_name}:')
print(model)

predictions = celltypist.annotate(countmtx, model=model, majority_voting=True)  # majority_voting refines cell identities within local subclusters but increases runtime
print('\nPredictions:')
print(predictions)
#print(predictions.predicted_labels)

# Save the predictions DFs
predictions.predicted_labels.to_csv('predicted_labels.tsv.gz', sep='\t')
predictions.decision_matrix.to_csv(f'decision_matrix.tsv.gz', sep='\t')
predictions.probability_matrix.to_csv(f'probability_matrix.tsv.gz', sep='\t')

# Add the majority-voted cell type labels to the count matrix
countmtx.obs = (countmtx.obs
    .merge(predictions.predicted_labels[['majority_voting']], left_index=True, right_index=True, how='inner', validate='1:1')
    .rename(columns={'majority_voting': 'celltype'})
)

# Sanitize the cell type labels
countmtx.obs['celltype'] = countmtx.obs['celltype'].str.replace(' ', '_').str.replace('/', '_').str.replace('(', '_').str.replace(')', '_')

# Make the gene_ids column the index
countmtx.var = countmtx.var.rename(columns={'gene_id': 'gene_ids'})
countmtx.var.set_index('gene_ids', inplace=True, verify_integrity=True, drop=False)

countmtx.write_h5ad(outpath, compression='gzip')

# Make some plots
adata = predictions.to_adata()

sc.tl.umap(adata)
sc.pl.umap(adata, color=['predicted_labels', 'majority_voting'], legend_loc='on data', save=f'.png')

# Make two plots side by side: predicted_labels and majority_voting
fig, axs = plt.subplots(1, 2, figsize=(18, 9))

for i, colname in enumerate(['predicted_labels', 'majority_voting']):
    value_counts = adata.obs[colname].value_counts()
    bars = axs[i].barh(value_counts.index, value_counts)
    axs[i].set_title(colname)
    axs[i].set_xlabel('# cells')
    axs[i].invert_yaxis()

    # Add text with the exact count value to the end of each bar
    for bar in bars:
        axs[i].text(bar.get_width(), bar.get_y() + bar.get_height() / 2, 
                    f'{int(bar.get_width())}', 
                    va='center', ha='left')

plt.tight_layout()
plt.savefig(f'celltype_distribution.png')

# Add the pool label to metadata
cell_meta['pool'] = sample
cell_meta.to_csv(cell_meta_out, sep='\t', index=False)

#!/usr/bin/env python

import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import sys

countmtxs = sys.argv[1].split(' ')
cell_metas = sys.argv[2].split(' ')

# Concat the annotated cells into one df
annotated_cells = []
for countmtx in countmtxs:
    pool = countmtx.removeprefix('countmtx_w_celltypes.').removesuffix('.h5ad')
    adata = ad.read_h5ad(countmtx)
    adata.obs['pool'] = pool
    adata.obs.index.name = 'cell_id'
    annotated_cells.append(adata.obs.reset_index(drop=False).copy())
annotated_cells = pd.concat(annotated_cells, ignore_index=True)

# Concat the metadata cells into one df
meta_cells = []
for cell_meta in cell_metas:
    df = pd.read_csv(cell_meta, sep='\t')
    meta_cells.append(df)
meta_cells = pd.concat(meta_cells, ignore_index=True)

# Merge the shared cells' annotations
merged = (annotated_cells
    .merge(meta_cells, on=['cell_id', 'pool'], how='inner', validate='1:1', suffixes=('_annot', '_meta'))
    .sort_values(by=['celltype_meta', 'celltype_annot', 'pool', 'cell_id'])
)

# Save the merged data as well
merged[['cell_id', 'pool', 'celltype_meta', 'celltype_annot']].to_csv('celltype_concordance.tsv', sep='\t', index=False)

# Calculate the conf matrix values
conf_matrix = pd.crosstab(merged['celltype_meta'], merged['celltype_annot'])
conf_matrix_normalized = conf_matrix.div(conf_matrix.sum(axis=1), axis=0) * 100

# Plot the conf matrix
plt.figure(figsize=(len(merged['celltype_annot'].unique())*0.8, len(merged['celltype_meta'].unique())*0.7))
sns.heatmap(conf_matrix_normalized, annot=True, fmt='.1f', cmap='Greys', mask=(conf_matrix_normalized < 0.1), linecolor='lightgray', linewidths=0.5)
plt.title(f'Confusion matrix of metadata and annotated celltypes in {merged.shape[0]:,} shared cells.\nRows sum to 100%')
plt.xlabel('CellTypist')
plt.ylabel('Metadata')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('confusion_matrix.png', dpi=300)

# Also plot the counts of cells across celltypes in annotated and metadata cells
for suffix, label in zip(['annot', 'meta'], ['CellTypist', 'Metadata']):
    plt.figure(figsize=(8, max(2, len(merged[f'celltype_{suffix}'].unique()) * 0.2)))
    counts = merged[f'celltype_{suffix}'].value_counts(sort=False).sort_index()
    plt.barh(counts.index, counts.values, color='gray')
    plt.title(f'Distribution of cell types in {label} cells')
    plt.xlabel('Cells')
    plt.tight_layout()
    plt.savefig(f'celltype_distribution.{suffix}.png', dpi=300)

#!/usr/bin/env python

import scanpy as sc
import matplotlib.pyplot as plt

# Read 10x dataset
adata = sc.read_10x_h5('neuron_10k_v3_filtered_feature_bc_matrix.h5')

# Make variable names (in this case the genes) unique
adata.var_names_make_unique()

############
## STEP 1 ##
############

### PCA BEFORE FILTERING ###

fig, ax = plt.subplots(nrows=2, figsize=(5,10))
sc.tl.pca(adata)
sc.pl.pca(adata, na_color='blue', ax=ax[0], show=False)
ax[0].set_box_aspect(1)
ax[0].set_title('Before filtering')

### FILTERING ###

sc.pp.recipe_zheng17(adata)                            

### PCA AFTER FILTERING ###

sc.tl.pca(adata)
sc.pl.pca(adata, na_color='blue', ax=ax[1], show=False)
ax[1].set_box_aspect(1)
ax[1].set_title('After filtering')
plt.tight_layout()
plt.savefig('pca_plots.png', dpi=150)


############
## STEP 2 ##
############

### LEIDEN CLUSTERING ###

sc.pp.neighbors(adata)
sc.tl.leiden(adata)

## PLOTTING TSNE AND UMAP ##

fig, ax = plt.subplots(ncols = 2, figsize=(10,5))

# t-SNE
sc.tl.tsne(adata)
sc.pl.tsne(adata, ax=ax[0], show=False, color='leiden', legend_loc=None)
ax[0].set_box_aspect(1)
ax[0].set_title('t-SNE plot')

# UMAP
sc.tl.umap(adata, maxiter=1000)
sc.pl.umap(adata, ax=ax[1], show=False, color='leiden')
ax[1].set_box_aspect(1)
ax[1].set_title('UMAP plot')

# save figure
plt.tight_layout()
plt.savefig('tsne_umap.png', dpi=150)






























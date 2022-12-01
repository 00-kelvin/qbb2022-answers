#!/usr/bin/env python

import scanpy as sc
import matplotlib.pyplot as plt

# Read 10x dataset
adata = sc.read_10x_h5("neuron_10k_v3_filtered_feature_bc_matrix.h5")

# Make variable names (in this case the genes) unique
adata.var_names_make_unique()

############
## STEP 1 ##
############

### PCA BEFORE FILTERING ###

fig, ax = plt.subplots(nrows=2)
pca_data = sc.tl.pca(adata, copy=True)
sc.pl.pca(pca_data, na_color="blue", ax=ax[0], show=False)
ax[0].set_box_aspect(1)
ax[0].set_title('Before filtering')

### FILTERING ###

sc.pp.recipe_zheng17(adata)                            

### PCA AFTER FILTERING ###

pca_data_fltrd = sc.tl.pca(adata, copy=True)
sc.pl.pca(pca_data_fltrd, na_color="blue", ax=ax[1], show=False)
ax[1].set_box_aspect(1)
ax[1].set_title('After filtering')
plt.tight_layout()
plt.savefig("pca_plots.png", dpi=150)

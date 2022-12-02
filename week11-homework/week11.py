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
sc.pl.tsne(adata, ax=ax[0], color='leiden', legend_loc=None, show=False)
ax[0].set_box_aspect(1)
ax[0].set_title('t-SNE plot')

# UMAP
sc.tl.umap(adata, maxiter=1000)
sc.pl.umap(adata, ax=ax[1], color='leiden', show=False)
ax[1].set_box_aspect(1)
ax[1].set_title('UMAP plot')

plt.tight_layout()
plt.savefig('tsne_umap.png', dpi=150)

############
## STEP 3 ##
############

# t-test method
fig, ax = plt.subplots()
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
sc.pl.rank_genes_groups(adata, ax=ax, show=False)
ax.set_title('Distinguising genes using t-test')
plt.savefig('t-test_groups.png')

# logistic regression method
fig, ax = plt.subplots()
sc.tl.rank_genes_groups(adata, groupby='leiden', method='logreg')
sc.pl.rank_genes_groups(adata, ax=ax, show=False)
ax.set_title('Distinguising genes using logistic regression')
plt.savefig('logreg_groups.png')

############
## STEP 4 ##
############

# These are the marker genes I decided on...
marker_genes_dict = {
    'EC2': ['Tac2'],
    'OL': ['Olig1'],
    'PC': ['Rgs5', 'Cox4i2'],
    'MG': ['C1qc', 'Skap2', 'Tyrobp'],
    'EC3': ['Hbb-bt'],
    'FB1': ['Tshz1']
}

# Making some support t-SNE plots
genelist = ['Tac2', 'Cox4i2', 'Olig1', 'Rgs5', 'C1qc', 'Skap2',
    'Tyrobp', 'Hbb-bt', 'Tshz1']

for i in range(0, len(genelist)):
    sc.pl.tsne(adata, color=genelist[i], 
                save='_{}.png'.format(genelist[i]), show=False)

# Support dotplot
sc.pl.dotplot(adata, marker_genes_dict, 'leiden', 
                dendrogram=True, save=True, show=False)

# Create a dictionary to map cluster to annotation label
cluster2annotation = {
     '20': 'EC2',
     '21': 'OL',
     '25': 'PC',
     '26': 'MG',
     '16': 'EC3',
     '19': 'FB1'
}

# Add a new '.obs' column called 'cell type'
adata.obs['cell type'] = adata.obs['leiden'].map(cluster2annotation).astype('category')

# UMAP of my cell types 
sc.pl.umap(adata, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10, legend_fontoutline=2, 
           save='cell_types.png', show=False)





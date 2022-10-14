#!/usr/bin/env python

# import matplotlib and numpy
import matplotlib.pyplot as plt
import numpy as np

# generate pc coordinates from text file
pc_coord = np.genfromtxt("plink.eigenvec", 
								dtype = None,
								encoding = None,
								usecols = (2, 3),
								names = ["PC1", "PC2"]) 

# make my first scatter plot using appropriate columns from the pc_coord data
fig, ax = plt.subplots()
ax.scatter(pc_coord["PC1"], pc_coord["PC2"], label = "PC1 vs PC2")

# label my axes
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")

ax.set_title("Genetic relatedness between the cell lines")

# save figure
plt.savefig("pc1_vs_pc2.png", dpi=200)






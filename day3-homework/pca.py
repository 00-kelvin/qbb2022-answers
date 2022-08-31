#!/usr/bin/env python

# import the matplotlib and numpy
import matplotlib.pyplot as plt
import numpy as np

# generate pc coordinates from text file
pc_coord = np.genfromtxt("plink.eigenvec", 
								dtype = None,
								encoding = None,
								names = ["indiv", "fam", "PC1", "PC2", "PC3"]) 

# make my first scatter plot using appropriate columns from the pc_coord data
fig, ax = plt.subplots()
ax.scatter(pc_coord["PC1"], pc_coord["PC2"], label = "PC1 vs PC2")

# label my axes
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")

# save figure
plt.savefig("ex2_a.png")

# make my second scatter plot using appropriate columns from the pc_coord data
fig, ax = plt.subplots()
ax.scatter(pc_coord["PC1"], pc_coord["PC3"], label = "PC1 vs PC3")

# label my axes
ax.set_xlabel("PC1")
ax.set_ylabel("PC3")

# save figure
plt.savefig("ex2_b.png")


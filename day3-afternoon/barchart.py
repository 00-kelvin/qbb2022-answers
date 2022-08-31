#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

genes_per_chrom = np.genfromtxt("genes_per_chrom.txt", 
								dtype = None,			#let numpy guess
								encoding = None,
								names = ["gene_count", "chrom_num"])

fig, ax = plt.subplots()
ax.bar(genes_per_chrom["chrom_num"], genes_per_chrom["gene_count"])
plt.show()
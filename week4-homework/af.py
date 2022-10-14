#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

snp_af = np.genfromtxt("plink.frq", 
							dtype = float,
							encoding = None,
							skip_header = 1,
							usecols = (4))

fig, ax = plt.subplots()

ax.hist(snp_af, bins = 50)
ax.set_xlabel("Allele frequency")
ax.set_ylabel("Count")
ax.set_title("Allele frequency spectrum")

plt.savefig("af.png", dpi=200)
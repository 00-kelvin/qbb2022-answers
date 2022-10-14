#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

snp_af = np.genfromtxt("plink.frq", 
							dtype = float,
							encoding = None,
							skip_header = 1,
							usecols = (4))

print(snp_af[:20])

plt.hist(snp_af, bins = 50)

plt.savefig("af.png")
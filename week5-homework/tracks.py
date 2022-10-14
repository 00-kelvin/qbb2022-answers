#!/usr/bin/env python

import matplotlib.pyplot as plt
from bdg_loader import load_data

# load data

sox2 = load_data("callpeak-out/D2_Sox2_R1_treat_pileup_cropped.bdg")
klf4 = load_data("D2_Klf4_treat_cropped.bdg")
H3K27ac0 = load_data("D0_H3K27ac_treat_cropped.bdg")
H3K27ac2 = load_data("D2_H3K27ac_treat_cropped.bdg")

# make plots

titles = ["Sox2", "Klf4", "H3K27ac day 0", "H3K27ac day 2"]

fig, ax = plt.subplots(nrows = 4, figsize=(8, 8))

for i, data in enumerate([sox2, klf4, H3K27ac0, H3K27ac2]):
	ax[i].plot(data['X'], data['Y'])
	ax[i].fill_between(data['X'], data['Y'])
	ax[i].set_title(titles[i])

plt.tight_layout()

plt.savefig("tracks.pdf")
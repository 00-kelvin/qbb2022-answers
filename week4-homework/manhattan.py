#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_CB = pd.read_csv('CB1908_IC50_gwas_results.assoc.linear', 
	delim_whitespace=True)
df_GS = pd.read_csv('GS451_IC50_gwas_results.assoc.linear', 
	delim_whitespace=True)

# only care about actual genotypes not covariate p-values

df_CB = df_CB[df_CB['TEST'] == 'ADD']
df_GS = df_GS[df_GS['TEST'] == 'ADD']


titles = ['CB1908 IC50 GWAS Results', 'GS451 IC50 GWAS Results']

for i, df in enumerate([df_CB, df_GS]):
	df['minus_log10_pvalue'] = -np.log10(df.P)
	df.CHR = df.CHR.astype('category')
	df['index'] = range(len(df))
	df_grouped = df.groupby(('CHR'))
	df_sig = df[df['minus_log10_pvalue'] > 5]

	fig, ax = plt.subplots(figsize=(12, 6)) 
	colors = ['orange','gold','lightgreen', 'lightblue']

	x_labels = []
	x_ticks = []

	for j, (name, group) in enumerate(df_grouped):
	    group.plot(kind='scatter', x='index', y='minus_log10_pvalue',
			color=colors[j % len(colors)], ax=ax)
	    x_labels.append(name)
	    x_ticks.append((group['index'].iloc[0] + (group['index'].iloc[-1] - group['index'].iloc[0])/2))

	df_sig.plot(kind='scatter', x='index', y='minus_log10_pvalue', ax=ax)

	ax.set_title(titles[i])
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_labels, fontsize = 9)
	ax.set_xlabel('Chromosome')
	plt.tight_layout()

	plt.savefig(titles[i] + '.png', dpi = 100)


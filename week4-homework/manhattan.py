#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_CB = pd.read_csv('CB1908_IC50_gwas_results.assoc.linear', 
	delim_whitespace=True)
df_GS = pd.read_csv('GS451_IC50_gwas_results.assoc.linear', 
	delim_whitespace=True)

for df in [df_CB, df_GS]:
	df['minuslog10pvalue'] = -np.log10(df.P)
	df.CHR = df.CHR.astype('category')
	df['index'] = range(len(df))
	df_grouped = df.groupby(('CHR'))

	fig = plt.figure(figsize=(12, 6)) 
	ax = fig.add_subplot(111)
	colors = ['orange','gold','lightgreen', 'lightblue']

	for i, (name, group) in enumerate(df_grouped):
	    group.plot(kind='scatter', x='index', y='minuslog10pvalue',
	    	color=colors[i % len(colors)], ax=ax)

	plt.show()


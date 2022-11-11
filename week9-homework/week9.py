#!/usr/bin/env python

import numpy as np

input_arr = np.genfromtxt("dros_gene_expression.csv", 
 							delimiter=',', names=True, 
 							dtype=None, encoding='utf-8')

col_names = list(input_arr.dtype.names)[1:]
row_names = list(input_arr['t_name'])

fpkm_values = input_arr[col_names]

print(fpkm_values)
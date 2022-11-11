#!/usr/bin/env python

import numpy as np
import numpy.lib.recfunctions as rfn

############################
##                        ##
##  0a: read in the data  ##
##                        ##
############################

input_arr = np.genfromtxt("dros_gene_expression.csv", 
 							delimiter=',', names=True, 
 							dtype=None, encoding='utf-8')

# sample names only 
col_names = list(input_arr.dtype.names)[1:]
row_names = list(input_arr['t_name'])

# subset to exclude transcript names
fpkm_values = input_arr[col_names]

##############################
##                          ##
##  0b: process input data  ##
##                          ##
##############################

fpkm_values_2d = rfn.structured_to_unstructured(fpkm_values, dtype=float)

# indices of transcripts with median expression greater than 0
nonzero_ind = np.where(np.median(fpkm_values_2d, axis=1) > 0)[0]

# filter data and transcript names
fpkm_fltrd = fpkm_values_2d[nonzero_ind]
row_names_fltrd = [row_names[i] for i in nonzero_ind]

# log transform
fpkm_fltrd_log = np.log2(fpkm_fltrd + 0.1)


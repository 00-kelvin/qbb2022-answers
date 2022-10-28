#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns

mat_fname, out_fname = sys.argv[1:3]

mat_data = np.loadtxt(mat_fname, dtype=np.dtype([
    ('F1', int), ('F2', int), ('score', float)]))

start_bin = 54878
end_bin = 54951

fltrd = mat_data[np.where((mat_data['F1'] >= start_bin) & 
                            (mat_data['F2'] <= end_bin))]

j = 1000

for row in fltrd:
    row[2] = np.log(row[2])

    # this is my bootleg way of finding the minimum score...
    if row[2] < j:
        j = row[2]

# subtract out the lowest score and lowest bins
min_bin = fltrd[0][0]

# used this 'k' code to find the max scores to set the color map bar
# k = 0 
for row in fltrd:
    row[0] = row[0] - min_bin
    row[1] = row[1] - min_bin
    row[2] = row[2] - j

#     if row[2] > k:
#         k = row[2]

# print(k)

#gives the maximum bin after the min bin has been subtracted out
new_max = fltrd[len(fltrd) - 1][1]

mat = np.zeros((new_max + 1, new_max + 1))

# for each row in the filtered/log tfmd data file, set 
# the entry in the matrix where n = F1 and m = F2 to the score
for row in fltrd:
    mat[row[0], row[1]] = row[2]
    mat[row[1], row[0]] = row[2]

fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(5,6.25))
sns.heatmap(mat, ax=ax[0], vmin=0, vmax=10, cmap="magma_r", square=True,
                xticklabels=False, yticklabels=False, cbar_kws={"label": "Score"})

ax[0].set_title('dCTCF Interactions')
ax[0].set_xlabel('chr15:10400000-13400000')
ax[0].set_ylabel('chr15:10400000-13400000')

x = np.arange(5, new_max - 4)
y = []

for i in x:
    y.append(np.mean(mat[(i - 5):i, i:(i + 5)]))

ax[0].axis('off')
plt.margins(x=0)

ax[1].plot(x,y)
ax[1].set_xlim(0, new_max)
ax[1].set_title('Insulation score')
plt.subplots_adjust(left=0.15,
                bottom=0.1,
                right=1.0,
                top=1.0,
                wspace=0.4,
                hspace=0.0)



plt.tight_layout()
plt.savefig(out_fname, dpi = 200)

#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns


def main():
    # in1_fname should be ddCTCF
    # in2_fname should be dCTCF
    # bin_fname should be bed file with bin locations
    
    in1_fname, in2_fname, bin_fname, out_fname = sys.argv[1:5]
    data1 = np.loadtxt(in1_fname, dtype=np.dtype([
        ('F1', int), ('F2', int), ('score', float)]))
    data2 = np.loadtxt(in2_fname, dtype=np.dtype([
        ('F1', int), ('F2', int), ('score', float)]))
    frags = np.loadtxt(bin_fname, dtype=np.dtype([
        ('chr', 'S5'), ('start', int), ('end', int), ('bin', int)]))

    chrom = b'chr15'
    start = 11170245
    end = 12070245

    # print(data1)

    start_bin = frags['bin'][np.where((frags['chr'] == chrom) &
                                         (frags['start'] <= start) &
                                         (frags['end'] > start))[0][0]]
    end_bin = frags['bin'][np.where((frags['chr'] == chrom) &
                                       (frags['start'] <= end) &
                                       (frags['end'] > end))[0][0]] - 1

    d = {}

    titles = ['ddCTCF', 'dCTCF', 'dCTCF - ddCTCF']

    fig, ax = plt.subplots(ncols = 3, figsize = (14, 2.8))

    for i, file in enumerate([data1, data2]):

        # filter any out of range data
        fltrd = file[np.where((file['F1'] >= start_bin) & (file['F2'] <= end_bin))]

        j = 1000

        for row in fltrd:
            row[2] = np.log10(row[2])

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

            # if row[2] > k:
            #     k = row[2]

        # print(k)

        #gives the maximum bin after the min bin has been subtracted out
        new_max = fltrd[len(fltrd) - 1][1]

        mat = np.zeros((new_max + 1, new_max + 1))

        # for each row in the filtered/log tfmd/subtracted data file, set 
        # the entry in the matrix where n = F1 and m = F2 to the score
        for row in fltrd:
            mat[row[0], row[1]] = row[2]
            mat[row[1], row[0]] = row[2]

        d['mat{0}'.format(i + 1)] = mat

        # set up my heatmap
        sns.heatmap(mat, ax=ax[i], vmin=0, vmax=4, cmap="rocket_r", square=True,
                xticklabels=False, yticklabels=False, cbar_kws={"label": "Score"})

        ax[i].set_title(titles[i])
        ax[i].set_xlabel('chr15:11170245-12070245')
        ax[i].set_ylabel('chr15:11170245-12070245')


    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns

def remove_dd_bg(mat):
    N = mat.shape[0]
    mat2 = np.copy(mat)
    for i in range(N):
        bg = np.mean(mat[np.arange(i, N), np.arange(N - i)])
        mat2[np.arange(i, N), np.arange(N - i)] -= bg
        if i > 0:
            mat2[np.arange(N - i), np.arange(i, N)] -= bg
    return mat2

def smooth_matrix(mat):
    N = mat.shape[0]
    invalid = np.where(mat[1:-1, 1:-1] == 0)
    nmat = np.zeros((N - 2, N - 2), float)
    for i in range(3):
        for j in range(3):
            nmat += mat[i:(N - 2 + i), j:(N - 2 + j)]
    nmat /= 9
    nmat[invalid] = 0
    return nmat

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

        # for each row in the filtered/log tfmd/subtracted data file, set 
        # the entry in the matrix where n = F1 and m = F2 to the score
        for row in fltrd:
            mat[row[0], row[1]] = row[2]
            mat[row[1], row[0]] = row[2]

        

        # set up my heatmap
        sns.heatmap(mat, ax=ax[i], vmin=0, vmax=10, cmap="magma_r", square=True,
                xticklabels=False, yticklabels=False, cbar_kws={"label": "Score"})

        ax[i].set_title(titles[i])
        ax[i].set_xlabel('chr15:11170245-12070245')
        ax[i].set_ylabel('chr15:11170245-12070245')

        # remove distance dependent signal and smooth
        final_mat = smooth_matrix(remove_dd_bg(mat))

        # save each of my edited matrices in a dictionary for later use
        d['mat{0}'.format(i + 1)] = final_mat

    dif_mat = d['mat2'] - d['mat1']


    sns.heatmap(dif_mat, ax=ax[2], cmap="seismic", square=True,
                xticklabels=False, yticklabels=False, cbar_kws={"label": "Score"},
                norm=colors.CenteredNorm())

    ax[2].set_title(titles[2])
    ax[2].set_xlabel('chr15:11170245-12070245')
    ax[2].set_ylabel('chr15:11170245-12070245')

    plt.tight_layout()
    plt.savefig(out_fname, dpi = 200)

if __name__ == "__main__":
    main()
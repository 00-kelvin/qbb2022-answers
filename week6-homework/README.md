# Week 6 Homework

## Part 1: Getting, Exploring, and Commenting on the HiCPro analysis results

**What percentage of reads are valid interactions (duplicates do not count as valid)?**

* For the dCTCF data, 37.8% of read pairs are valid interaction pairs, and 92% of those "valid pairs" are "valid interactions" (not duplicates), which means ~34.8% of the reads are valid interactions.
* For the ddCTCF data, 36.6% of read pairs are valid interaction pairs, and 88% of those "valid pairs" are "valid interactions" (not duplicates), which means ~32.2% of the reads are valid interactions.

**What constitutes the majority of invalid 3C pairs?**

According to the Hi-C manual, the following are considered invalid pairs:

> * Dangling end, i.e. unligated fragments (both reads mapped on the same restriction fragment)
> * Self circles, i.e. fragments ligated on themselves (both reads mapped on the same restriction fragment in inverted orientation
> * Religation, i.e. ligation of juxtaposed fragments
> * Dumped pairs, i.e. any pairs that do not match the filtering criteria on inserts size, restriction fragments size or for which we were not able to reconstruct the ligation product.

I am not sure what the majority of the invalid pairs would be, but I would guess that the dumped pairs make up a large fraction, as it seems like any of a number of missed criteria could result in a pair being dumped.

## Part 2

```
./load_data.py analysis/hic_results/matrix/ddCTCF/iced/6400/ddCTCF_ontarget_6400_iced.matrix analysis/hic_results/matrix/dCTCF/iced/6400/dCTCF_ontarget_6400_iced.matrix analysis/hic_results/matrix/dCTCF/raw/6400/dCTCF_ontarget_6400_abs.bed subsampled.png

./load_data.py matrix/ddCTCF_full.6400.matrix matrix/dCTCF_full.6400.matrix matrix/6400_bins.bed full.png
```

1. Yes, I am able to see the highlighted difference!
2. Compared to the subsampled dataset, the full dataset had stronger interaction scores between most or all regions of the chromosome, and there appeared to be less extreme differences between the dCTCF and ddCTCF interaction profiles.
3. The highlighted signal indicates the position of the CTCF site.
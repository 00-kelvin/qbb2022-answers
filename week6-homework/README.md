# Week 6 Homework

## Part 1: Getting, Exploring, and Commenting on the HiCPro analysis results

**What percentage of reads are valid interactions (duplicates do not count as valid)?**

* For the dCTCF data, 37.8% of read pairs are valid interaction pairs, and 92% of those "valid pairs" are "valid interactions" (not duplicates), which means ~34.8% of the reads are valid interactions.
* For the ddCTCF data, 36.6% of read pairs are valid interaction pairs, and 88% of those "valid pairs" are "valid interactions" (not duplicates), which means ~32.2% of the reads are valid interactions.

**What constitutes the majority of invalid 3C pairs?**

According to the Hi-C manual, the following are considered invalid pairs:

> Dangling end, i.e. unligated fragments (both reads mapped on the same restriction fragment)
> Self circles, i.e. fragments ligated on themselves (both reads mapped on the same restriction fragment in inverted orientation
> Religation, i.e. ligation of juxtaposed fragments
> Dumped pairs, i.e. any pairs that do not match the filtering criteria on inserts size, restriction fragments size or for which we were not able to reconstruct the ligation product.
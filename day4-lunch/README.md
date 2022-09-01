#QBB2022 - Day 4 - Lunch Exercises Submission

**Exercise 1**

Portion of do_all.sh output (not script) that reports how many bp each feature covers:

```
--- Subsetting exons.chr21.bed.vcf
    + Covering 1107407 bp
--- Subsetting processed_pseudogene.chr21.bed.vcf
    + Covering 956640 bp
--- Subsetting protein_coding.chr21.bed.vcf
    + Covering 13780687 bp
```

One or more strategies to confirm that reproduced plots are the same as in the cache/ directory:

* After running the code, open the plots copied from the cache/ directory and the plots generated by running the code and compare visually
* Write a script that opens and reads both png files, and tests whether they are the same using the == operator, and returns True or prints "same" if they are

Three other gene_types present in the GENCODE .gtf that you find interesting and why:

* transcribed_unprocessed\_pseudogene: these are pseudogenes that have not been moved away from their parent gene, but are still transcribed, which means maybe they are used for something. interesting because generally pseudogenes are only transcribed if they are processed?
* miRNA: non-coding regulatory RNAs. maybe interesting to compare to the parts of the genome that they interact with 
* lncRNA: always interesting because we don't know what a lot of them do
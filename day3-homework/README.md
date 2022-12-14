# QBB2022 - Day 3 - Homework Exercises Submission

**Exercise 1**

Code to run PCA: 

```
plink --vcf ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz --pca 3
```

**Exercise 2**

My script for these plots is called "pca.py" and the plots themselves are "ex2_a-b.png"

There definitely is some structure among the points in both plots: both plots show a sort of crooked "L" shape; in both graphs, along PC1, the component showing the greatest variation, there is a large cluster towards the lower end and another cluster towards the positive end. In the PC1vPC2 graph, there appear to be two major clusters along PC2, and in the PC1vPC3 graph, there is a more uniform distribution with possibly 2 major clusters. This indicates that there are groups of the population that are more similar to each other, but since I don't yet have any information about the individuals, I can't connect that to any particular biological concept. Possibly the clusters could represent genetic similarity within geographical areas.

**Exercise 3**

To join the files: 

```
join -1 1 -2 1 <(sort /Users/cmdb/data/metadata_and_txt_files/integrated_call_samples.panel) <(sort plink.eigenvec) > joined.txt
```

What a doozy! My script for these plots is called "pca2.py" and the plots themselves are "ex3_a-c.png"
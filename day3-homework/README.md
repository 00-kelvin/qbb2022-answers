# QBB2022 - Day 3 - Homework Exercises Submission

**Exercise 1**

Code to run PCA: 

```
plink --vcf ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz --pca 3
```

**Exercise 2**

There definitely is some structure among the points in both plots: both plots show a sort of crooked "L" shape; in both graphs, along PC1, the component showing the greatest variation, there is a large cluster towards the lower end and another cluster towards the positive end. In the PC1vPC2 graph, there appear to be two major clusters along PC2, and in the PC1vPC3 graph, there is a more uniform distribution with possibly 2 major clusters. This indicates that there are groups of the population that are more similar to each other, but since I don't yet have any information about the individuals, I can't connect that to any particular biological concept. Possibly the clusters could represent genetic similarity within geographical areas.
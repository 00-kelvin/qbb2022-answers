 # QBB2022 - Day 1 - Lunch Exercises Submission

 1. I’m excited to learn how to create beautiful figures with Python.

 2. Mean number of exons per gene: 62

commands: 
```
$ wc exons.chr21.bed
$ wc genes.chr21.bed
$ echo $((13653/219)) 
```

To find the median: //
a) Find the number of exons in each gene using the start/stop positions in the exon file compared to start/stop positions in the gene file //
b) Store the exon counts in a new file //
c) sort -n (sort the list of counts in numerical order) //
d) use head 110 | tail 1 to print the 110th (middle) entry of the count list
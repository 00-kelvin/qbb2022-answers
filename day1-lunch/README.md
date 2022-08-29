 # QBB2022 - Day 1 - Lunch Exercises Submission

 1. I’m excited to learn how to create beautiful figures with Python.

 2. Mean number of exons per gene: 62

commands: 
```
$ wc exons.chr21.bed
$ wc genes.chr21.bed
$ echo $((13653/219)) 
```

To find the median:

* Find the number of exons in each gene using the start/stop positions in the exon file compared to start/stop positions in the gene file  
* Store the exon counts in a new file  
* sort -n (sort the list of counts in numerical order)  
* use head 110 | tail 1 to print the 110th (middle) entry of the count list  
3. Regions classified for each state: 
```
 305 chr21	1
 678 chr21	2
  79 chr21	3
 377 chr21	4
 808 chr21	5
 148 chr21	6
1050 chr21	7
 156 chr21	8
 654 chr21	9
  17 chr21	10
  17 chr21	11
  30 chr21	12
  62 chr21	13
 228 chr21	14
 992 chr21	15
```
command used: 
```
sort -nk 4 chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21.bed | uniq -f 3 -c | cut -f 1,4
```

To determine which state comprises largest fraction of genome:
* For each state, find the sum of the lengths of all the genes in that state
* Store sums as a column in a file with their corresponding states
* Sort numerically by sums
* use tail 1 to print the state with the largest sum
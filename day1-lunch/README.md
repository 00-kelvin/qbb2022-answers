 # QBB2022 - Day 1 - Lunch Exercises Submission

 1. I’m excited to learn how to create beautiful figures with Python.

 2. Mean number of exons per gene: 62

commands: 
```
wc exons.chr21.bed
wc genes.chr21.bed
echo $((13653/219)) 
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
* use tail -1 to print the state with the largest sum

4. Number of samples in each pop of AFR:

```
 123 HG01880	ACB	AFR
 112 NA19625	ASW	AFR
 173 HG02922	ESN	AFR
 180 HG02462	GWD	AFR
 122 NA19017	LWK	AFR
 128 HG03052	MSL	AFR
 206 NA18484	YRI	AFR
 ```
 
 commands used: 
 
 ```
 grep "AFR" integrated_call_samples.panel | sort -k 2 | cut -f -3 | uniq -cf 1
 ```
 
 For all 5 populations: 
 * Sort first by super population
 * Then sort by subpopulation as above
 * Use the uniq function to collapse the subpopulations and -c to get a count for how many samples are in each
 
 5. Explore SNP allele frequencies
 
 to create the HG00100 file:
 ```
 cut -f 1-9,13 random_snippet.vcf > HG00100.vcf
 ```

 counts for 0|0, 0|1, 1|0, and 1|1 values:
 ```
 9514 0|0
  127 0|1
  178 1|0
  181 1|1
  ```
 
 command to find the counts:
 ```
 sort -k 10 HG00100.vcf | cut -f 10 | uniq -c | head -17 | tail -4
 ```
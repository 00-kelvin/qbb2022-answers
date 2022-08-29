# QBB2022 - Day 1 - Homework Exercises Submission

**Exercise 1**

The error message we see is 'awk: illegal field $(), name "nuc"' for each of A, G, C, and T

To fix it, we need to set the identity of the variable with the -v flag within the awk command

Output from working script:

```
Considering  A
 354 C
1315 G
 358 T
Considering  C
 484 A
 384 G
2113 T
Considering  G
2041 A
 405 C
 485 T
Considering  T
 358 A
1317 C
 386 G
 ```
 
 It makes sense that A would be most commonly replaced with G (and vice versa) and C with T (and vice versa) since A and G are purines and thus similar sizes/shapes, while C and T are pyrimidines
 
 **Exercise 2**
 
 Promoters do not seem to be clearly and objectively defined: other things like transcription start sites (TSS) are.
 
 Chose state 2, "flanking active TSS" as best approximation of promoter.
 
 Created a new bed file that only has state 2 genes: 
 
 ```
 awk '{if ($4 == 2) {print}}' ~/data/bed_files/chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21.bed > genefile_state2.bed
 ```
 
 Then passed the new bedfile and the random_snippet.vcf through the bedtools intersect function, saving the output as intersect_out_ex2.vcf:
 
 ```
 bedtools intersect -a ~/data/vcf_files/random_snippet.vcf -b genefile_state2.bed > intersect_out_ex2.vcf
 ```
 
 We want the output as a vcf file rather than bed because we need information about each variant, not just number of variants.
 
 Finally, used the following command to print the alternates for all Cytosine reference alleles, sort them, collapse them to unique values, and count the frequency of each alternate: 
 
 ```
 awk '/^#/{next} {if ($4 == "C") {print $5}}' intersect_out_ex2.vcf | sort | uniq -c
 ```
 
 Resulting in: 
 ```
   7 A
   4 G
  24 T
 ```
 These results are very low numbers of variants in this region. One hypothesis might be that this area is highly conserved because of its importance in initiating transcription.
 
 **Exercise 3**
 
 Line 5 (the awk statement) creates a new bed file from the input VCF with 3 columns: the chromosome (chr21), the position minus 1, and the position of the variant. This is because the closest function requires bed file inputs, and bed files need a start and stop position (variants only have 1 position so pos-1, pos is a good approximation).
 
 Line 6 sorts the genes.bed file first by chromosome, then numerically by the starting position and saves it as a new file called 'genes.sorted.bed'. This is necessary because according to the bedtools documentation, "bedtools closest requires that all input files are presorted data by chromosome and then by start position".
 
 Line 7 feeds the new variants file and the sorted genes file into the bedtools closest function.
 
 The first error is Error: unable to open file or unable to determine types for file variants.bed. it suggests: - Please ensure that your file is TAB delimited (e.g., cat -t FILE).
 
 I changed the awk statement to make the output tab-delimited:
 
 ```
 awk '/^#/{next} BEGIN{OFS="\t"} {print $1,$2-1, $2}' $1 > variants.bed
 ```
 The second error was Error: Sorted input specified, but the file variants.bed has the following out of order record
 
 I added a line to sort the variants.bed file and save it as a new file, variants.sorted.bed which I then fed into the bedtools closest function.
 
 There are 10293 variants returned: 
 ```
 bash exercise3.sh ~/data/vcf_files/random_snippet.vcf | wc -l
    10293
 ```
 Used the following command to pull out the gene names, collapse to unique values, then count the number of unique genes, of which there are 200: 
 ```
 bash exercise3.sh ~/data/vcf_files/random_snippet.vcf | awk '{print $7}' | sort | uniq | wc -l 
      200
 ```
 On average, 10293/200 = ~51 variants are connected to each gene.
 
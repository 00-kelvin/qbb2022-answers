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
 These results are very low numbers of SNPs in this region. One hypothesis might be that this area is highly conserved because of its importance in initiating transcription.
 
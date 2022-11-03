# Week 8 Homework


**Command to call variants:**

```
for REGION in chr11:1900000-2800000 chr14:100700000-100990000 chr15:23600000-25900000 chr20:58800000-58912000; do medaka_variant -i methylation.bam -f hg38.fa -r ${REGION} -p -o ${REGION} -t 6; done
```

**Commands to mark reads with haplotype tags:**

```
whatshap haplotag -o chr11_marked_reads.bam --reference hg38.fa --regions chr11:1900000:2800000 --output-haplotag-list chr11_htlist.txt chr11:1900000-2800000/round_0_hap_mixed_phased.vcf.gz methylation.bam
```

```
whatshap haplotag -o chr14_marked_reads.bam --reference hg38.fa --regions chr14:100700000:100990000 --output-haplotag-list chr14_htlist.txt chr14:100700000-100990000/round_0_hap_mixed_phased.vcf.gz methylation.bam
```

```
whatshap haplotag -o chr15_marked_reads.bam --reference hg38.fa --regions chr15:23600000:25900000 --output-haplotag-list chr15_htlist.txt chr15:23600000-25900000/round_0_hap_mixed_phased.vcf.gz methylation.bam
```

```
whatshap haplotag -o chr20_marked_reads.bam --reference hg38.fa --regions chr20:58800000:58912000 --output-haplotag-list chr20_htlist.txt chr20:58800000-58912000/round_0_hap_mixed_phased.vcf.gz methylation.bam
```

**Command to split reads:**

```
for CHR in chr11 chr14 chr15 chr20; do whatshap split --output-h1 ${CHR}_marked_reads_h1.bam --output-h2 ${CHR}_marked_reads_h2.bam ${CHR}_marked_reads.bam ${CHR}_htlist.txt; done
```

**Commands to concatenate bam files:**

```
samtools cat -o h1_marked_reads.bam chr11_marked_reads_h1.bam chr14_marked_reads_h1.bam chr15_marked_reads_h1.bam chr20_marked_reads_h1.bam
```

```
samtools cat -o h2_marked_reads.bam chr11_marked_reads_h2.bam chr14_marked_reads_h2.bam chr15_marked_reads_h2.bam chr20_marked_reads_h2.bam
```

**Do you expect each region in H1 or H2 to correspond to the same parent of origin (i.e. the same haplotype)? Explain your reasoning.**

Not necessarily: since I ran whatshap haplotag separately on each of the 4 regions, it could have assigned H1 as the mother and H2 as the father in one region, but the opposite in another region. (It's possible that whatshap has a standardized way of doing things such that this doesn't happen, but I can't find any evidence of that it its documentation.)
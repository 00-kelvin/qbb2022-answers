# Week 8 Homework


**Command to call variants:**
```for REGION in chr11:1900000-2800000 chr14:100700000-100990000 chr15:23600000-25900000 chr20:58800000-58912000; do medaka_variant -i methylation.bam -f hg38.fa -r ${REGION} -p -o ${REGION} -t 6; done```

**Commands to mark reads with haplotype tags:**
```whatshap haplotag -o chr11_marked_reads.bam --reference hg38.fa --regions chr11:1900000:2800000 --output-haplotag-list chr11_htlist.txt chr11:1900000-2800000/round_0_hap_mixed_phased.vcf.gz methylation.bam```

```whatshap haplotag -o chr14_marked_reads.bam --reference hg38.fa --regions chr14:100700000:100990000 --output-haplotag-list chr14_htlist.txt chr14:100700000-100990000/round_0_hap_mixed_phased.vcf.gz methylation.bam```

```whatshap haplotag -o chr15_marked_reads.bam --reference hg38.fa --regions chr15:23600000:25900000 --output-haplotag-list chr15_htlist.txt chr15:23600000-25900000/round_0_hap_mixed_phased.vcf.gz methylation.bam```

```whatshap haplotag -o chr20_marked_reads.bam --reference hg38.fa --regions chr20:58800000:58912000 --output-haplotag-list chr20_htlist.txt chr20:58800000-58912000/round_0_hap_mixed_phased.vcf.gz methylation.bam```
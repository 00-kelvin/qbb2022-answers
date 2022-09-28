#!/bin/bash

# # Step 1: index the reference genome
# bwa index sacCer3.fa 

# # Step 2: align the reads with the reference genome
# for SAMPLE in A01_09 A01_11 A01_23 A01_24 A01_27 A01_31 A01_35 A01_39 A01_62 A01_63
# do
# 	bwa mem \
# 	-t 4 \
# 	-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
# 	-o ${SAMPLE}.sam \
# 	sacCer3.fa ${SAMPLE}.fastq
# done
	
# # Step 3a: create sorted bam files
# for SAMPLE in A01_09 A01_11 A01_23 A01_24 A01_27 A01_31 A01_35 A01_39 A01_62 A01_63
# do
# 	samtools sort -@4 -O bam -o ${SAMPLE}.bam ${SAMPLE}.sam
# done
	
# # Step 3b: create an index for each sorted bam file
# for SAMPLE in A01_09 A01_11 A01_23 A01_24 A01_27 A01_31 A01_35 A01_39 A01_62 A01_63
# do
# 	samtools index ${SAMPLE}.bam
# done

# Step 4-5: variant calling with freebayes, filter by genotype quality
freebayes -f sacCer3.fa \
	-p 1 \
	--genotype-qualities \
	*.bam \
	| vcffilter -f "QUAL > 20" \
	| vcfallelicprimitives -k -g \
	> sacCer3.vcf

# Step 6: 

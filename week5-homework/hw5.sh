#!/bin/bash

for FILE in D2_Sox2_R1 D2_Sox2_R1_input D2_Sox2_R2 D2_Sox2_R2_input 
do 
	samtools view \
	-q 10 \
	-o ${FILE}_filtered.bam \
	${FILE}.bam
done
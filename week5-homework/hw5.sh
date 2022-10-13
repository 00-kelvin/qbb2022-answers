#!/bin/bash

for READ in D2_Sox2_R1 D2_Sox2_R1_input D2_Sox2_R2 D2_Sox2_R2_input 
do 
	samtools view \
	-q 10 \
	-o ${READ}_fltrd.bam \
	${READ}.bam
done

for REPLICATE in D2_Sox2_R1 D2_Sox2_R2
do
	macs2 callpeak \
	-t ${REPLICATE}_fltrd.bam \
	-c ${REPLICATE}_input_fltrd.bam \
	-g 94987271 \
	-B 
done

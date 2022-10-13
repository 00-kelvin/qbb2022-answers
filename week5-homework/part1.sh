#!/bin/bash

# Step 1: Filtering reads

for READ in D2_Sox2_R1 D2_Sox2_R1_input D2_Sox2_R2 D2_Sox2_R2_input 
do 
	samtools view \
	-q 10 \
	-b \
	-o ${READ}_fltrd.bam \
	${READ}.bam
done

# Step 2: Calling peaks

for REPLICATE in D2_Sox2_R1 D2_Sox2_R2
do
	macs2 callpeak \
	-t ${REPLICATE}_fltrd.bam \
	-c ${REPLICATE}_input_fltrd.bam \
	-g 94987271 \
	-B \
	-n ${REPLICATE} \
	--outdir ./callpeak-out/
done

# Step 3: Intersecting peaks

bedtools intersect -wa \
-a ./callpeak-out/D2_Sox2_R1_peaks.narrowPeak \
-b ./callpeak-out/D2_Sox2_R2_peaks.narrowPeak \
> D2_Sox2_peaks.bed

# Step 4: Colocalization of Sox2 and Klf4

echo "Number of Sox2 peaks:" > peak_counts.txt
wc -l < D2_Sox2_peaks.bed >> peak_counts.txt

echo "Number of Klf4 peaks:" >> peak_counts.txt
wc -l < D2_Klf4_peaks.bed >> peak_counts.txt

echo "Number of overlapping peaks:" >> peak_counts.txt
bedtools intersect -wa -a D2_Klf4_peaks.bed -b D2_Sox2_peaks.bed | wc -l \
>> peak_counts.txt

echo "Percentage of overlapping peaks:" >> peak_counts.txt
awk "BEGIN {print (41/60)*100}" >> peak_counts.txt

# Step 5: Plot

for BDG in D2_Klf4_treat D0_H3K27ac_treat D2_H3K27ac_treat callpeak-out/D2_Sox2_R1_treat_pileup
do
	python scale_bdg.py ${BDG}.bdg ${BDG}_scaled.bdg
	awk '{ if ($2 < 35507055 && $3 > 35502055) print $0 }' ${BDG}_scaled.bdg > ${BDG}_cropped.bdg
done
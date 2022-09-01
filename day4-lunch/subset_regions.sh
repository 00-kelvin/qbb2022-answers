#!/bin/bash

GTF=$1
CHR=chr21

if [ ! -f $CHR.gtf ]
then
    echo "--- Creating $CHR.gtf"
    grep -w $CHR $GTF > $CHR.gtf
fi

# add lncRNA
for TYPE in protein_coding processed_pseudogene lncRNA
do
    echo "--- Creating $TYPE.$CHR.bed"
    grep $TYPE $CHR.gtf | awk 'BEGIN {OFS="\t"}{if ($3=="gene"){print $1, $4-1, $5}}' > $TYPE.$CHR.bed
done

echo "--- Creating exons.$CHR.bed"
grep "protein_coding" $CHR.gtf | awk 'BEGIN {OFS="\t"}{if ($3=="exon"){print $1, $4-1, $5}}' > exons.$CHR.bed


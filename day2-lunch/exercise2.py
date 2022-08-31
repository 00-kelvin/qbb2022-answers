#!/usr/bin/env python3
#usage: exercise2.py bedfile.bed

#Import BED parser function
from bed_parser import parse_bed

#Load a BED file
import sys
bed = parse_bed(sys.argv[1])

#Find the median number of exons for genes in the BED file 

#Initialize empty list for sorting
blockCountList = []

#Add block counts from each line to empty list
for i in range(len(bed)):
	blockCountList.append(bed[i][9])

#Sort blockCountList
blockCountList.sort()

#Find the median if the length of the bed file is even
if len(blockCountList) % 2 == 0:
	pos1 = int(len(blockCountList) / 2)
	pos2 = pos1 + 1
	median = (blockCountList[pos1] + blockCountList[pos2]) / 2

#if there are an odd number of records
else: 
	median = blockCountList[int(len(blockCountList) / 2 + 0.5)]

#print median number
print(f"The median number of exons is {median}")


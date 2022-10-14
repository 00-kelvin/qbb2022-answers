# Week 5 Homework

## Part 1

Command line commands are in the part1.sh file.

Found effective genome size from: https://genome-test.gi.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A1%2D94987271&hgsid=396441278_yk94l0t9HFqDzc04ZhxvPmFKFNGD


## Part 2

Commands: 

```
sort -k 5,5rn D2_Sox2_peaks.bed | head -300 | awk '{ printf "%s:%i-%i\n", $1, $2, $3 }' > D2_Sox2_peaks_srtd.bed

samtools faidx mm10.fa -r D2_Sox2_peaks_srtd.bed -o peak_sequences.fa

```

then, after activating meme environment: 

`meme-chip -maxw 7 peak_sequences.fa`

## Part 3


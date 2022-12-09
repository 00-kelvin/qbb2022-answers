# Week 12 Homework

## Step 1.

### Step 1B. 

I used the provided code, saved as ```parse_kraken.py```

```
for KRAKEN in 83 86 88 89 90 93 94 97; do ./parse_kraken.py ./metagenomics_data/step0_givendata/KRAKEN/SRR4921${KRAKEN}.kraken SRR4921${KRAKEN}; done
```

### Step 1C. 

Rank ```ktImportText``` on the samples

```
ktImportText -q *_krona.txt
```

**QUESTION 1**

Throughout the first week, _Enterococcus faecalis_ is the most prevalent group of bacteria, with the other major groups being _Bacillales_ (mostly including different _staphyloccocus_ strains) and _Actinobacteria_. On the first day the amount of _Bacillales_ and _Actinobacteria_ are approximately the same (~18% of total) with the rest being _Enterococcus faecalis_. Then the _Actinobacteria_ group collapses entirely and _Bacillales_ decreases, with _Enterococcus faecalis_ taking up 92% of the samples. Over the rest of the week, _Bacillales_ gradually increases again, and in the last 2 days, _Actinobacteria_ reappears and then rapidly increases to 27% on the last day.

## Step 2.

**QUESTION 2**
??

### Step 2A.

```
bwa index metagenomics_data/step0_givendata/assembly.fasta
```

### Step 2B.

From my step0_givendata directory (in which I added a ```sorted_bam_files``` directory): 

```
for SAMPLE in 83 86 88 89 90 93 94 97; do bwa mem -t 4 assembly.fasta ./READS/SRR4921${SAMPLE}_1.fastq ./READS/SRR4921${SAMPLE}_2.fastq | samtools sort -o ./sorted_bam_files/SRR4921${SAMPLE}.bam; done
```

### Step 2D.

```
jgi_summarize_bam_contig_depths --outputDepth depth.txt ./sorted_bam_files/*.bam

metabat2 -i assembly.fasta -a depth.txt -o bins_dir/bin

```

**QUESTION 3A**
I got 6 bins

**QUESTION 3B**
There are 197 contigs in the 6 bin files: 

```
grep '>' bins_dir/*.fa | wc -l
     197
```

and 4103 in the whole assembly: 

```
grep '>' assembly.fasta | wc -l
    4103
```	

So these represent only about 4.8% of contigs in the whole the assembly.

In terms of length of the whole assembly, which is 38708237 bases:
```
grep -v ">" assembly.fasta | wc -m
 38708237
```
the proportion is much higher, the bins take up 17365118 bases:

```
grep -v ">" bins_dir/bin.*.fa | wc -m
 17365118
```

which is 44.9%.

**QUESTION 3C**

I found the lengths of all the contigs in each of the 6 bins:
```
for BIN in bins_dir/*; do grep -v ">" ${BIN} | wc -m; done
 2750133
 2289418
 1683638
 1248387
 2525062
 2910568
 ```

Prokaryotic genomes tend to be 0.6 to 8.0 Mb; these values are all within that range, so they seem reasonable.

**QUESTION 3D**
?? 

## Step 3. 

Made lists of contigs for each bin (have to cut out the ">" so it can be used to ```grep``` in the assembly.kraken):

```
for BIN in 1 2 3 4 5 6; do grep ">" bins_dir/bin.${BIN}.fa | cut -c 2- > ./contig_lists/bin${BIN}.txt; done
```
then search for those contigs in the assembly.kraken and save all the matches to a file, then collapse to unique matches for each (with counts):

```
for BIN in 1 2 3 4 5 6; do grep -f contig_lists/bin${BIN}.txt ./KRAKEN/assembly.kraken | cut -f 2 | sort | uniq -c > ./taxonomies/bin${BIN}.txt; done
```


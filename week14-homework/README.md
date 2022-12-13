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

We could use amount of coverage per contig, codon usage, GC content/nucleotide frequency, and relative abundance over a time series to group them.

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

Once we have a guess of what species each bin corresponds to, we could align the bins to existing genomes for those species and compare to figure out what contaminations there are and what parts are missing.

## Step 3. 

Made lists of contigs for each bin (have to cut out the ">" so it can be used to ```grep``` in the assembly.kraken):

```
for BIN in 1 2 3 4 5 6; do grep ">" bins_dir/bin.${BIN}.fa | cut -c 2- > ./contig_lists/bin${BIN}.txt; done
```
then search for those contigs in the assembly.kraken and save all the matches to a file, then collapse to unique matches for each (with counts):

```
for BIN in 1 2 3 4 5 6; do grep -f contig_lists/bin${BIN}.txt ./KRAKEN/assembly.kraken | cut -f 2 | sort | uniq -c > ./taxonomies/bin${BIN}.txt; done
```

From these results I chose the most frequent hit as my prediction (used ```wc -l``` to find total number of contigs in each bin for the proportions reported below).

**QUESTION 4A**

Predictions for each bin:

* **Bin 1:** 49/55 root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus ST72;Staphylococcus aureus subsp. aureus CN1
	* All the contigs (55/55) fall under the same hierarchies up to _Staphylococcus aureus subsp. aureus_, so we can say Bin 1 is _Staphylococcus aureus subsp. aureus_ with high confidence
* **Bin 2:** (51/78) root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus epidermidis;Staphylococcus epidermidis RP62A
	* Nearly all contigs (77/78) == _Staphylococcus epidermidis_, the only other 1 is _Staphylococcus aureus_, so we can predict with high confidence that Bin 2 corresponds to _Staphylococcus epidermidis_ and with even higher confidence that it it is part of the _Staphylococcus_ genus
* **Bin 3:** (3/8) root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Tissierellia;Tissierellales;Peptoniphilaceae;Anaerococcus;Anaerococcus prevotii;Anaerococcus prevotii DSM 20548
	* Since there are only 8 contigs in this bin, and they are somewhat spread out, the only prediction we can make with high confidence is that Bin 3 corresponds to a species in the _Firmicutes_ phylum
* **Bin 4:** (24/37) root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus haemolyticus;Staphylococcus haemolyticus JCSC1435
	* The rest of the contigs that do not belong to this species are other _Staphylococcus_ bacteria, but not _Staphylococcus haemolyticus_
* **Bin 5:** (13/13) root;cellular organisms;Bacteria;Terrabacteria group;Actinobacteria;Actinobacteria;Propionibacteriales;Propionibacteriaceae;Cutibacterium;Cutibacterium avidum;Cutibacterium avidum 44067
	* All contigs agree!
* **Bin 6:** all 6 belong to species _Enterococcus faecalis_, two strains have 2 hits each: OG1RF and V583.

**QUESTION 4B**

We could more robustly infer the taxonomy by calculating some measure of significance of the match, such as an e-value. 


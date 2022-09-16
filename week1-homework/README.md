# Week 1 Homework

## Exercise 1

### Question 1: Coverage Simulator

#### Question 1.1

1 Mbp genome x 5x coverage = 5 Mbp of data

5 Mbp / 100 bp reads = 50,000 reads

For 15x coverage, 150,000 reads

#### Question 1.3

In the 5x coverage simulation, 6953 bp were not sequenced. This is comparable to the Poisson expectations, which I calculated in my script to be 6737.946999085467.

#### Question 1.4

In the 15x coverage simulation, 3 bp were not sequenced. The Poisson expectation was that only 0.3059023205018258 bp would not be sequenced, so 3 is quite close to this expectation of basically no missed bp.

### Question 2: De novo assembly

The bash script I wrote to find these answers is saved in this directory under `/asm/asm/exercise2.sh`; it worked when saved in the same directory as the fasta files (`~/qbb2022-answers/week1-homework/asm/asm`) and produced the following output: 

```
Question 2.1:
4
Question 2.2:
234467
Question 2.3:
105830
```

#### Question 2.1

4 contigs were produced.

#### Question 2.2

The total length of the contigs is 234467

#### Question 2.3

The largest contig is 105830

#### Question 2.4

Since the total length is 234467, 50% of the genome is 117233.5 bp. The contigs ordered from largest to smallest have lengths: 
```
105830
47860
41351
39426
```
The first contig does not contain 50% of the genome, but the first 2 do (together they contain 153690 bp). Therefore, the N50 value is 47860.

### Question 3: Whole Genome Alignment

#### Question 3.1

I ran: `$ dnadiff ./asm/ref.fa ./asm/asm/contigs.fasta`

According to the out.report, the average identity of my assembly compared to the reference was 100.00

#### Question 3.2

The output of

`$ nucmer ./asm/ref.fa ./asm/asm/contigs.fasta`
`$ show-coords out.delta`

was

```
    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
=====================================================================================
  127965   233794  |        1   105830  |   105830   105830  |    99.99  | Halomonas	NODE_1_length_105830_cov_20.649193
   40651    88510  |        1    47860  |    47860    47860  |   100.00  | Halomonas	NODE_2_length_47860_cov_20.367392
       3    26789  |        1    26787  |    26787    26787  |   100.00  | Halomonas	NODE_3_length_41351_cov_20.528098
   26790    40641  |    27500    41351  |    13852    13852  |   100.00  | Halomonas	NODE_3_length_41351_cov_20.528098
   88532   127957  |        1    39426  |    39426    39426  |   100.00  | Halomonas	NODE_4_length_39426_cov_20.336388
```

It appears that the length of the longest alignment was 105830.

#### Question 3.3

From the out.report, there is 1 insertion of length 712 and 5 deletions of total length 51.

### Question 4: Decoding the insertion

#### Question 4.1

The insertion begins at position 26788 in Node 3 of the assembly. This corresponds to position 26790 in the reference.

#### Question 4.2

The insertion is 712 bp long.

#### Question 4.3

`samtools faidx contigs.fasta NODE_3_length_41351_cov_20.528098:26788-27500` gives the following result: 

>NODE_3_length_41351_cov_20.528098:26788-27500
CGCCCATGCGTAGGGGCTTCTTTAATTACTTGATTGACGCATGCCCCTCGTTCTACATGT
CTAGCTTCGTAACTGCCCCGATTTATACAGGAGCATATGCGTTTCGTAGTGCCGGGAATG
CATACCAAAGGGCTCACGGCGGGTACGCCACAATGGCTCAAGTCGAAAATGAATCGAAGA
CAACAAGGAATACCGTACCCAATTACTCAAGGACCTCATACACCATCCCATGCTACTTAT
CTACAGACATACACGCCAGCACCCAGCAACCAAAGCACACCGACGATAAGACTACAATCG
CGATAAGCACAACTTACATTAGGAGGCCCGGCAAATCTTGACGGCGTTAAGTCCGACACG
AATACCCCCCGACAAAAGCCTCGTATTCCGAGAGTACGAGAGTGCACAAAGCACCAAGGC
GGGGCTTCGGTACATCCACCAGTAGTCCCGTCGTGGCGGATTTTCGTCGCGGATGATCCG
AGGATTTCCTGCCTTGCCGAACACCTTACGTCATTCGGGGATGTCATAAAGCCAAACTTA
GGCAAGTAGAAGATGGAGCACGGTCTAAAGGATTAAAGTCCTCGAATAACAAAGGACTGG
AGTGCCTCAGGCATCTCTGCCGATCTGATTGCAAGAAAAAATGACAATATTAGTAAATTA
GCCTATGAATAGCGGCTTTAAGTTAATGCCGAGGTCAATATTGACATCGGTAG

#### Question 4.4

I saved the above sequence to a fasta file: `samtools faidx contigs.fasta NODE_3_length_41351_cov_20.528098:26788-27500 > insertion.fasta` then ran `python ../dna-decode.py --decode --input insertion.fasta`

The decoded message :  Congratulations to the 2021 CMDB @ JHU class!  Keep on looking for little green aliens...



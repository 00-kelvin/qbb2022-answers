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

The bash script I wrote to find these answers is saved as exercise2.sh; it works when saved in the same directory as the fasta files (`~/qbb2022-answers/week1-homework/asm/asm`) and produces the following output: 

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



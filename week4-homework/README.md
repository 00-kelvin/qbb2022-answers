## Week 4 Homework

### Question 2

Command to generate PC values: `plink --vcf ./gwas_data/genotypes.vcf --pca 10`

### Question 3

To generate allele frequencies: `plink --vcf ./gwas_data/genotypes.vcf --freq`

### Question 4

Quantitative association test: `plink --vcf ./gwas_data/genotypes.vcf --linear --pheno ./gwas_data/CB1908_IC50.txt --covar plink.eigenvec --allow-no-sex --out CB1908_IC50_gwas_results` and `plink --vcf ./gwas_data/genotypes.vcf --linear --pheno ./gwas_data/GS451_IC50.txt --covar plink.eigenvec --allow-no-sex --out GS451_IC50_gwas_results`

### Question 5

For the manhattan plots, I drew a lot of inspiration from here: https://www.python-graph-gallery.com/manhattan-plot-with-matplotlib 

I rewrote all the code myself and made sure I know what's going on in it but let me know if it's too close.

### Question 7

**Top SNP in the CB1908 dataset: rs10876043**
* According to the UCSC Genome Browser, this SNP may be affecting the gene DIP2B which codes for Disco-interacting protein 2 homolog B. This protein "contains a binding site for the transcriptional regulator DNA methyltransferase 1 associated protein 1 as well as AMP-binding sites" which suggests that it may be involved in DNA methylation
* If this protein is mutated, it is possible that genes that should be methylated and thus repressed will not be. The expression of these usually repressed genes could lead to lymphocytopenia upon interaction with the drug.
* Code to locate this SNP is in the effect_size.py script


**Top SNP in the GS451 dataset: rs7257475**
* This SNP appears to be in the gene ZNF826. This is a gene for a zinc finger protein, but according to the genome browser, it is a pseudogene (non-coding RNA). 
* I am not sure how this could contribute to the disease phenotype. Perhaps the SNP somehow makes this gene active again, and it contributes to the lymphocytopenia phenotype. 
* Code to located this SNP is in the manhattan.py script


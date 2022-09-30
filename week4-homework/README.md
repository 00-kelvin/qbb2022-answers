## Week 4 Homework

Command to generate PC values: `plink --vcf ./gwas_data/genotypes.vcf --pca 10`

To generate allele frequencies: `plink --vcf ./gwas_data/genotypes.vcf --freq`

Quantitative association test: `plink --vcf ./gwas_data/genotypes.vcf --linear --pheno ./gwas_data/CB1908_IC50.txt --covar plink.eigenvec --allow-no-sex` and `plink --vcf ./gwas_data/genotypes.vcf --linear --pheno ./gwas_data/GS451_IC50.txt --covar plink.eigenvec --allow-no-sex`


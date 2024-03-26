#!/bin/bash
cd ~/Bioinformatics/macaco/all_cynos

# Phenotype table as input and aggregated vcf
pheno_table="all_cynos_pheno.txt"

# file must be bgzipped
VCF="raw_vcf.vcf.gz"

# index using tabix
bcftools index -t $VCF

# Filter VCF down to the samples of interest only
# We are using all samples in this VCF so skip

# Filter VCF down to chr1-20, X, Y, MT
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,X,Y $VCF -o vcf_qc1.vcf.gz

# PLINK --vcf convert to --bfile .bed
plink --vcf vcf_qc1.vcf.gz --make-bed --out cynos

# PLINK filter out SNP's below a represented frequency of 0.8
plink --bfile cynos --maf 0.8 --make-bed --out cynos_maf

# Assign phenotypes with PLINK
plink --bfile cynos_maf --pheno $pheno_table --make-bed --out cynos_w_pheno --allow-no-sex

# QC step 1: --mind 0.1, missing genotype data which means exclude samples with more than 10% missing genotypes overall
plink --bfile cynos_w_pheno --mind 0.1 --make-bed --out cynos_qc1 --allow-no-sex

# PLINK QC genotype reads
plink --bfile cynos_qc1 --geno 0.1 --make-bed --out cynos_qc2 --allow-no-sex

# PLINK QC HWE
plink --bfile cynos_qc2 --hwe 1e-6 --make-bed --out cynos_qc3

# Assoc is used as our disease phenotype is either healthy or sick (binary)
# PLINK allele frequency test
plink --bfile cynos_qc3 --assoc --allow-no-sex --out cynos_assoc

#* Make sure to --allow-no-sex
# To check p-values, cat cynos_assoc.assoc | awk 'BEGIN {OFS="\t"}; {print $9}' | sort | uniq

# Run python program on this code
cd ..

python3 bin/manhattan_plot.py -g MMul_10/ncbi_dataset/data/GCF_003339765.1/genomic.gff -a all_cynos/cynos_assoc.assoc
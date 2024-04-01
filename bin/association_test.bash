#!/bin/bash

#########################################
############ VARIABLE NAMING ############
#########################################
# Running directory
curr_dir=$(pwd)

#Folder name of your output where you want to have all the files created and have placed the four required files, gff, acc_to_chr, pheno.txt, and .vcf
folder_name="all_cynos"

#path to where your files are located
folder_path=${curr_dir}/${folder_name}

cd $folder_path

# Phenotype table as input and aggregated vcf
pheno_table=${folder_path}/pheno.txt

# GFF file path
GFF=${folder_path}/MMul_10/ncbi_dataset/data/GCF_003339765.1/genomic.gff

# file must be bgzipped
VCF=${folder_path}/raw_vcf.vcf.gz

# Accession No to Chr csv
acc_to_chr=${folder_path}/accession_no_to_Chrom.csv

#############################################
############ BASH SCRIPT SECTION ############
#############################################

# index using tabix
bcftools index -t $VCF

# Filter VCF down to the samples of interest only
# We are using all samples in this VCF so skip

# Filter VCF down to chr1-20, X, Y, MT
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,X,Y $VCF -o plink_pre.vcf.gz

# PLINK --vcf convert to --bfile .bed
plink --vcf plink_pre.vcf.gz --make-bed --out plink_file

# PLINK filter out SNP's below a represented frequency of 0.01 (default) and this really improved things
plink --bfile plink_file --maf 0.01 --make-bed --out plink_qc1

# Assign phenotypes with PLINK
plink --bfile plink_qc1 --pheno $pheno_table --make-bed --out plink_pheno_qc1 --allow-no-sex

# QC step 1: --mind 0.1, missing genotype data which means exclude samples with more than 10% missing genotypes overall
plink --bfile plink_pheno_qc1 --mind 0.1 --make-bed --out plink_qc2 --allow-no-sex

# PLINK QC genotype reads
plink --bfile plink_qc2 --geno 0.1 --make-bed --out plink_qc3 --allow-no-sex

# PLINK QC HWE
plink --bfile plink_qc3 --hwe 1e-6 --make-bed --out plink_qc4

# Assoc is used as our disease phenotype is either healthy or sick (binary)
# PLINK allele frequency test
plink --bfile plink_qc4 --assoc --allow-no-sex --out plink_final

#* Make sure to --allow-no-sex
# To check p-values, cat cynos_assoc.assoc | awk 'BEGIN {OFS="\t"}; {print $9}' | sort | uniq

# Run python program
python3 ${curr_dir}/bin/manhattan_plot.py -ng $GFF -a $folder_path/plink_final.assoc -c $acc_to_chr -ng ${GFF} -og ${folder_path}/corrected.gff -v ${VCF} -o ${folder_path}/macic_results
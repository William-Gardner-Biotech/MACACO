```
                    ███╗   ███╗  █████╗   ██████╗ ██╗  ██████╗
                    ████╗ ████║ ██╔══██╗ ██╔════╝ ██║ ██╔════╝
                    ██╔████╔██║ ███████║ ██║      ██║ ██║     
                    ██║╚██╔╝██║ ██╔══██║ ██║      ██║ ██║     
                    ██║ ╚═╝ ██║ ██║  ██║ ╚██████╗ ██║ ╚██████╗
                    ╚═╝     ╚═╝ ╚═╝  ╚═╝  ╚═════╝ ╚═╝  ╚═════╝                  
```
![MACIC_image](https://github.com/William-Gardner-Biotech/MACIC/assets/99355149/0a37dfef-baa0-4af9-b456-e6c2b4f20f80)


Multiple Allelic Comparison to Identify Candidates (MACIC) is a Bioinformatic workflow designed with the express goal of running a GWAS (Genome Wide Association Study) on large VCF files from the same species to identify statistically significant SNPs. MACIC was built using two great open source tools, [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/) and [bcftools](https://github.com/samtools/bcftools), and uses Python to generate plots and assemble reports.



This program requires 4 input files to run:
1. A VCF file containing your samples of interest (Must be bgzipped).
2. A PLINK phenotype .txt file to annotate sample phenotypes.
  - Phenotype file is a tab separated .txt file
    
         Animal_ID  Animal_ID_with_family  Phenotype
         Animal_ID_1  Animal_ID_1(same)  1

  - Phenotype values: 1 = control, 2 = case, 0 = missing phenotype

3. The genbank annotation file (.gff) corresponding to the species.
4. A handmade .csv file that lists the Accession numbers and chromosomes contained within the .gff

This project is still a work in progress and is very rough around the edges but it is working as of 29Mar2024. 

Some people have asked, "Is it Macaque or Macic?" To which I have only one answer, it is Magic!

Assumptions made by this program:
1. Your input vcf file has standardized chromosome numbers i.e. chr1, 1, 01 and your accession_to_chromosome.txt file corresponds with this numbering convention.
2. Sample's sex is often hard to deduce using PLINK and unless known by other means (study notes) should be ignored. If you want to include sex as a confounder I would suggest referring to the plink documentation and adding sex labels using their built in method, then removing the "--allow-no-sex" flag from the bash script.

import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Convert the chromosome numbers into one long number set 1-3GbpMMul_10/ncbi_dataset/data/GCF_003339765.1/genomic.gff

# both functions lifted from xtract_genes_from_gff_by_gene_name.py

def build_acc_chrom_dict():
    with open('accession_no_to_Chrom.csv', 'r') as table:
        chrom_dict = {}
        for row in table:
            r_ow = row.split(',')
            # [0] is accession number [1] is chrom
            #print(f"ROW: {r_ow[1]}")
            chrom_dict[r_ow[0]] = r_ow[1].strip()
    return(chrom_dict)

# The Mmul_10 ref genome uses NC_ Accession numbers to name the chromosomes instead of numbering or chr1 etc. This command will replace those acc numbers with their respective chrom number for genes
def replace_acc_with_chrom(chrom_dict: dict):
    if os.path.exists(os.path.join(os.getcwd(), "Mmul10_chrom_converted.gff")):
        return print(f'Running redundant function, Converted GFF already exists {os.path.join(os.getcwd(), "Mmul10_chrom_converted.gff")}')
    else:
        with open("MMul_10/ncbi_dataset/data/GCF_003339765.1/genomic.gff", 'r') as filtered_gff, open("Mmul10_chrom_converted.gff", 'w') as out_bedfile:
            pattern = r"Name=([^;]+)"
            for entry in filtered_gff:
                # skip headers
                if entry.startswith('#'):
                    continue
                if "ID=gene" in entry:
                    ent_ry = entry.split('\t')
                    #[0]=accession,[3]=start,[4]=end,[5]=score,[6]=strand,[8]=description
                    #print(ent_ry)
                    try:
                        bed_line = f"{chrom_dict[ent_ry[0]]}\t{ent_ry[3]}\t{ent_ry[4]}\t{re.search(pattern, ent_ry[8]).group(1)}\t{ent_ry[5]}\t{ent_ry[6]}\n"
                    except: pass
                    out_bedfile.write(bed_line)

def chr_to_genomic_coords(genomic_gff):
    # Now we need to make a conversion that will build a list of all chromosome positions and add the 
    with open(genomic_gff, 'r') as gff:
        chroms = list(build_acc_chrom_dict().keys())
        # turn chroms into a list that can be compared against by "in" method
        chroms = [str(chrom) for chrom in chroms]
        #print(chroms)
        # A dictionary that matches each chromosome to its relative genomic start position
        genomic_chrom_start = {}
        genome_pos = 1 # GFF is indexed at 1
        for line in gff:
            # remove header lines
            if line.startswith('#'):
                continue
            if "ID=gene" in line:
                #print(line)
                pass
            # working to isolate the chromosomes
            elif "region" in line.split('\t')[2] and line.split('\t')[0] in chroms:
                #print(line)
                pattern = r"\=\w*\;"
                chromosome = re.search(pattern, line).group(0)
                if chromosome:
                    chromosome = chromosome.strip("=")
                    chromosome = chromosome.strip(";")
                    # print(chromosome)
                    # 5th line has chr end pos which I'll add to genome_pos after dict addition
                    chr_end = (line.split('\t')[4]) 
                    # key is chr number and value is its relative start position
                    genomic_chrom_start[chromosome] = genome_pos
                    genome_pos += int(chr_end)
                    print(f"Chr {chromosome}, start: {genomic_chrom_start[chromosome]}")
    return genomic_chrom_start

def gather_dots(genomic_chr_dict, association_file):
    """
    Function that will take every point on the .assoc file and simplify it to a pandas
    Dataframe object with the genomic position column and p score column
    """
    print(genomic_chr_dict.keys())

    # loop through the association file
    # CHR 23 is listed as the number for CHR X
    with open(association_file, 'r') as assoc:
        # iterate through and build lists of all 3 (plus POSition) to build into a pandas df
        CHRs = []
        SNPs = []
        POSs = []
        Ps = []
        for line in assoc:
            # 0 = CHR, 1 = SNP, 7 = P-value 
            CHR, SNP, _, _, _, _, _, _, P, _ = line.split()
            if len(SNP.split(':')) != 4:
                print(f'Excluding this line: {line}')
                continue
            else:
                if P == 'NA':
                    continue
                if CHR == '23':
                    CHR = 'X'
                CHRs.append(CHR)
                SNPs.append(SNP)
                Ps.append(P)
                # Split SNP by colon and retrieve SNP pos relative to the CHR it is on
                POSs.append(SNP.split(':')[1])

    fd = {'CHR': CHRs, 'SNP': SNPs, 'POS': POSs, 'p-value':Ps}
    df = pd.DataFrame(data=fd)

    print(df['CHR'].unique())

    # use CHR value as key to dict and add this to POS from the SNP
    def calculate_genomic_position(row):
        return genomic_chr_dict[row['CHR']] + int(row['POS']) - 1

    # Apply the function to create the new column
    df['Genomic_Position'] = df.apply(calculate_genomic_position, axis=1)

    # Add -log10 scale to make the smaller p-values stand out graphically
    df['neg_log10_p'] = df['p-value'].apply(lambda x: -np.log10(float(x)))

    return df

def generate_plot(df):
    '''
    Function that will generate the manhattan plot and add coloring based upon chromosome
    '''
    # Grouping SNPs by CHR will allow for coloring
    chr_groups = df.groupby('CHR')

    plt.figure(figsize=(20,12))

    chr_names = []
    chr_positions = []
    # .groupby() method will create an iterable of the unique group name along with the subset df of all data found in that group
    for chr_name, chr_df in chr_groups:
        plt.scatter(chr_df['Genomic_Position'], chr_df['neg_log10_p'], label=chr_name)

        chr_names.append(chr_name)
        chr_positions.append(chr_df['Genomic_Position'].mean()) # finding avg of the group to add the middle x-tick

    # Bonferroni is way too high, let's exclude all p-values == 1
    bad_observations = ((df['p-value']) == '1.0').sum()
    print(bad_observations)

    # Calculate the Bonferroni p-value to create threshold for p-value over false positives
    Bonferroni = 0.05 / (df.shape[0] - bad_observations)

    print(f'Bonferroni corrected p-value: {Bonferroni}')
    # Add a horizontal dotted line at the threshold -log10 value
    plt.axhline(y=(-np.log10(Bonferroni)), color='gray', linestyle=':', linewidth=1)
    plt.axhline(y=(-np.log10(0.05)), color='red', linestyle=':', linewidth=1)


    plt.xlabel('Chromosome')
    plt.ylabel('-log10 p-values')
    # Set x-axis ticks and labels
    #plt.xticks(chr_positions, chr_names)
    plt.show()




def main():
    # Chrom gene conversion done
    #x = build_acc_chrom_dict()
    #replace_acc_with_chrom(x)
    parser = argparse.ArgumentParser(
        description="Parser that will handle input files for the mahattan plot generation",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-g', '--gff', type=str, default='MMul_10/ncbi_dataset/data/GCF_003339765.1/genomic.gff', required=False,
                            help='Mmul10 gff to extract chromosome regions from and gene regions')

    parser.add_argument('-a', '--assoc', type=str, default='/Users/willgardner/Bioinformatics/macaco/all_cynos/cynos_assoc.assoc', required=False,
                            help='Association statistic file generated by the macaco nextflow script')
    args = parser.parse_args()

    # Build the relational dict so we can begin converting points for the manhattan plot
    genome_chr_coords = chr_to_genomic_coords(args.gff)

    plot_df = gather_dots(genome_chr_coords, args.assoc)

    generate_plot(plot_df)

main()
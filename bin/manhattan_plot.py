import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from accession_to_chr import build_acc_chrom_dict, replace_acc_with_chrom
from statsmodels.stats.multitest import multipletests

def chr_to_genomic_coords(genomic_gff, acc_to_chr):
    # Now we need to make a conversion that will build a list of all chromosome positions and add the 
    with open(genomic_gff, 'r') as gff:
        chroms = list(build_acc_chrom_dict(acc_to_chr).keys())
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

    df['p-value'] = df['p-value'].astype(float)

    # Adjust the p-values using Hol-Bonferroni
    def adjust_p_values(p_values, chosen_method='holm'):
        # function will return a tuple 0 = reject, 1 = adjusted p-value, 2 = Sidak adjusted, 3 = Bonferroni adjusted
        return multipletests(p_values, method = chosen_method, is_sorted = False, alpha = 0.05)
    
    df['reject'] = adjust_p_values(df['p-value'])[0]
    df['adj_p'] = adjust_p_values(df['p-value'])[1]

    #pd.set_option('display.max_rows', None)
    #print(df.head(500))

    #print('Adjusted')
    #print(df['adj_p'].unique())

    adj_df = (df[df['reject'] == True])
    adj_df['neg_log10_p'] = df['adj_p'].apply(lambda x: -np.log10(float(x)))

    return (df, adj_df)


def generate_plot(df, title = 'Manhattan plot (All p-values)', holm_bon = False):
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


    # Only do this section for first graph and skip on second
    if holm_bon == False:
        # Any obvious p-value ~> 0.75 should be discarded for bonferroni
        # df['float_p'] = df['p-value'].astype(float)
        obvious_obs = (df['p-value'] >= 0.5).sum()
        print(f"OBV OBS: {obvious_obs}")

        # Calculate the Bonferroni p-value to create threshold for p-value over false positives
        Bonferroni = 0.05 / (df.shape[0])

        print(f'Bonferroni corrected p-value: {Bonferroni}\nBonferroni minus Obvious OBV: {0.05 / (df.shape[0]- obvious_obs)}')
        # Add a horizontal dotted line at the threshold -log10 value
        plt.axhline(y=(-np.log10(Bonferroni)), color='gray', linestyle=':', linewidth=1)
        plt.axhline(y=(-np.log10(0.05)), color='red', linestyle=':', linewidth=1)

    plt.xlabel('Chromosome')
    plt.ylabel('-log10 p-values')
    # Set x-axis ticks and labels TODO return x axis to chr numbers
    #plt.xticks(chr_positions, chr_names)
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, labels, title='CHROMOSOME', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.title(title)

    plt.show()

    plt.clf()

    for chr_name, chr_df in chr_groups:
        plt.scatter(chr_df['Genomic_Position'], chr_df['neg_log10_p'], label=chr_name)

        chr_names.append(chr_name)
        chr_positions.append(chr_df['Genomic_Position'].mean())
    
    # TODO fixed 26Mar24 MAF = 0.05 to remove the streaks caused by low maf inclusion, left in as a way to visualize the SNP counts
    #pd.set_option('display.max_rows', None)
    #print(f'Value Counts:\n {df["p-value"].value_counts()}')
    #print('LOG 10 P == 9.338377:')
    #df['neg_log10_p'] = df['neg_log10_p'].astype(float)
    #print(df[(df['neg_log10_p']>=9.3) & (df['neg_log10_p'] <= 9.4)])
        
    def generate_report(df, alpha, gff_path:str):
        """
        This function will generate a report on all of the findings including, list of all SNP's identified as significant using Bonferroni,
        Holm-Bonferroni and list their corresponding gene regions, if any. It will also provide metadata statistics on the overall process run.
        Finally it will also include the SNP effect. It will all be gathered into one large pandas dataframe.
        """


        # First step is to find where the significant SNP's intersect gene regions






def main():
    # Chrom gene conversion done
    #x = build_acc_chrom_dict()
    #replace_acc_with_chrom(x)
    parser = argparse.ArgumentParser(
        description="Parser that will handle input files for the mahattan plot generation",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-ng', '--ncbi_gff', type=str, default='MMul_10/ncbi_dataset/data/GCF_003339765.1/genomic.gff', required=False,
                            help='Mmul10 gff to extract chromosome regions from and gene regions, this has all chromosomes as accession numbers not ')
    
    parser.add_argument('-og', '--out_gff', type=str, default='Mmul10_chrom_converted.gff', required=False,
                            help='Genomic Feature annotation file with corrected chromosome names.') 
    
    parser.add_argument('-c', '--acc_to_chr', type=str, default="accession_no_to_Chrom.csv", required=False,
                            help='This custom file is required for the pipeline to run. It is not generated by this pipeline and instead expected to be built from GenBank.\n\
                                It is a .csv file with two columns, ACC_NO and CHROM in that order. The program relies on them to be ordered this way to build the new gff file.\n\
                                This is necessary to annotate the SNPs and find the significant gene regions.')

    parser.add_argument('-a', '--assoc', type=str, default='/Users/willgardner/Bioinformatics/macaco/all_cynos/cynos_assoc.assoc', required=False,
                            help='Association statistic file generated by the macaco nextflow script')
    parser.add_argument('-p', '--alpha', type=float, default=0.05, required=False,
                            help='Desired alpha level used for the statistical hypthesis testing.\n\
                                DEFAULT: 0.05')
    args = parser.parse_args()

    # Build the converted gff file
    x = build_acc_chrom_dict(args.acc_to_chr)
    replace_acc_with_chrom(x, args.ncbi_gff, args.out_gff) # args.out_gff will be used later for opening with generate report

    # Build the relational dict so we can begin converting points for the manhattan plot
    genome_chr_coords = chr_to_genomic_coords(args.ncbi_gff, args.acc_to_chr)

    plot_df, adj_df = gather_dots(genome_chr_coords, args.assoc)

    # Default df plot
    generate_plot(plot_df)

    # Holm-Bonferroni Adj manhattan
    generate_plot(adj_df, 'Manhattan plot (Holm-Bon significant SNPs)', True)

main()
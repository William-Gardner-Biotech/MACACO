import os
import re

def build_acc_chrom_dict(acc_to_chrom):
    with open(acc_to_chrom, 'r') as table:
        chrom_dict = {}
        for row in table:
            r_ow = row.split(',')
            # [0] is accession number [1] is chrom
            #print(f"ROW: {r_ow[1]}")
            chrom_dict[r_ow[0]] = r_ow[1].strip()
    return(chrom_dict)

# The Mmul_10 ref genome uses NC_ Accession numbers to name the chromosomes instead of numbering or chr1 etc. This command will replace those acc numbers with their respective chrom number for genes
def replace_acc_with_chrom(chrom_dict: dict, in_gff = "MMul_10/ncbi_dataset/data/GCF_003339765.1/genomic.gff", out_gff = "Mmul10_chrom_converted.gff"):
    if os.path.exists(os.path.join(os.getcwd(), "Mmul10_chrom_converted.gff")):
        return print(f'Running redundant function, Converted GFF already exists {os.path.join(os.getcwd(), out_gff)}')
    else:
        with open(in_gff, 'r') as filtered_gff, open(out_gff, 'w') as out_bedfile:
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
                        # Easy to miss but in this line we use the chrom_dict to index and replace the chromosome accession with 1-23,X,Y
                        bed_line = f"{chrom_dict[ent_ry[0]]}\t{ent_ry[3]}\t{ent_ry[4]}\t{re.search(pattern, ent_ry[8]).group(1)}\t{ent_ry[5]}\t{ent_ry[6]}\n"
                    except: 
                        pass
                    out_bedfile.write(bed_line)

def main():
    """
    This function will take the CSV file that lists the Accession numbers and their corresponding chr and impute all parts of the gbff file with chr.
    NOTE: the accession_to_Chrom.csv file was handmade for Mmul_10 reference genome using genbank.
    """
    x = build_acc_chrom_dict('accession_no_to_Chrom.csv')
    replace_acc_with_chrom(x)

if __name__ == "__main__":
    main()
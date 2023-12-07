import argparse


# ! Args:
parser = argparse.ArgumentParser(description="Process input gene table, subset the gff file. 5 files will be generated")

### ~ Add the arguments
parser.add_argument("--target_table", type=str, required=True, help="Path to the gene information table.")
parser.add_argument("--gff", type=str, required=True, help="Path to the GFF3 file.")
parser.add_argument("--output_path", type=str, required=True, help="Path to the output file")
parser.add_argument("--output_prefix", type=str, required=True, help="Prefix of the output")

### ~ Parse the arguments
args = parser.parse_args()

# ! ~ Target output
### ~  5 subsets of interval list for
#### ~ exon: 
#### ~ intron:
#### ~ promoter:
#### ~ utr
#### ~ ncrna:


# ! Read the input file and transfer it to a
# ! dictionary for fast query
target_gene_dict = {}
with open(args.target_table) as f:
    next(f)
    for line in f:
        # @ gene, alt_name, intron, others
        line = line.strip().split("\t")
        dict_item = ["","",""]
        if len(line) > 1 :
            for i in range(len(line[1:])):
                dict_item[i] = line[i + 1]

        target_gene_dict[line[0]] = dict_item

# print(target_gene_dict)

# ! Prepare all the output files:
output_path = args.output_path
if output_path[-1] != "/":
    output_path += "/"

exon_out = open(f'{output_path}{args.output_prefix}.exon.interval_list.raw', "w+")
intron_out = open(f'{output_path}{args.output_prefix}.intron.interval_list.raw', "w+")
promoter_out = open(f'{output_path}{args.output_prefix}.promoter.interval_list.raw', "w+")
three_utr_out = open(f'{output_path}{args.output_prefix}.3utr.interval_list.raw', "w+")
ncrna_out = open(f'{output_path}{args.output_prefix}.ncrna.interval_list.raw', "w+")

# ! Open the 
with open(args.gff) as f:
    for line in f:
        if(line[0]) == "#":
            continue
        # @ chrm, _ , _ , start, end, _ , strand, _ , info
        line = line.split("\t")

        # @ Find the target gene 
        for current_gene in target_gene_dict.keys():
            if f';gene_name={current_gene};' not in line[8] and f';gene_name={target_gene_dict[current_gene][0]};' not in line[8]:
                continue
            else:
                # @ Whether non-extron/intron? :
                if target_gene_dict[current_gene][2] == "Promoter":
                    if line[2] == "five_prime_UTR":
                        promoter_out.write(f'{line[0]}\t{line[3]}\t{int(line[4]) - 1}\n')
                elif target_gene_dict[current_gene][2] == "3UTR":
                    if line[2] == "three_prime_UTR":
                        three_utr_out.write(f'{line[0]}\t{line[3]}\t{int(line[4]) - 1}\n')
                elif target_gene_dict[current_gene][2] == "ncRNA":
                    if line[2] == "gene" and f';gene_type=lncRNA;' in line[8]:
                        ncrna_out.write(f'{line[0]}\t{line[3]}\t{int(line[4]) - 1}\n')
                
                # @ Extron or intron?
                if target_gene_dict[current_gene][1] != "":
                    # * Extron:
                    if line[2] == "exon":
                        exon_out.write(f'{line[0]}\t{line[3]}\t{int(line[4]) - 1}\n')
                else:
                    # * Intron
                    if line[2] == "exon":
                        extend_info = line[8].strip().split(";")
                        extend_info = line[8].strip().split(";")
                        intron_out.write(f'{line[0]}\t{line[3]}\t{int(line[4]) - 1}\t{extend_info[5]}\t{extend_info[7]}\t{extend_info[8]}\n')
        